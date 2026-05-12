#include "catalog.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// --- Phase 3f.4: mmap-based catalog load ---
//
// We mmap the pair, k-vector, and pattern files PROT_READ MAP_PRIVATE,
// madvise(MADV_RANDOM) (all three are random-access binary-search targets),
// and mlock (best-effort) to pin them in physical memory so the first solve
// has the same latency as steady-state. mlock implicitly prefaults; on the
// mlock-failed fallback path we walk every page explicitly so the first
// solve doesn't pay a major-fault avalanche on the cold cache.
//
// The on-disk struct layouts must match the in-memory layouts byte-for-byte.
// Both target architectures (x86-64, arm64) are little-endian; the relevant
// types (double, int32, int64) have well-known sizes. These static_asserts
// fail the build loudly if anyone reorders fields or adds padding.
static_assert(sizeof(CatalogPair) == 16,
              "CatalogPair must be exactly 16 bytes (double + 2*int32) to "
              "match the on-disk layout written by tools/generate_catalog.py");
static_assert(sizeof(StarPattern) == 24,
              "StarPattern must be exactly 24 bytes (uint64 + 4*int32) to "
              "match the on-disk layout written by tools/generate_catalog.py");
static_assert(sizeof(int32_t) == sizeof(int),
              "This file assumes int is 32-bit. The on-disk layouts and the "
              "ifstream/mmap loaders both write/read raw int into 4-byte "
              "fields; a 16-bit or 64-bit int target would corrupt them.");

namespace {
constexpr size_t kPageSize = 4096; // assumed page granularity for prefault

// Touch every page of [ptr, ptr+size) so the kernel populates the page
// table eagerly. Without this, the first access pays a major page fault.
// Reading into a volatile sink prevents the optimizer from eliding the load.
void prefault(const uint8_t *ptr, size_t size) {
  volatile uint8_t sink = 0;
  for (size_t i = 0; i < size; i += kPageSize)
    sink = static_cast<uint8_t>(sink ^ ptr[i]);
  (void)sink;
}
} // namespace

std::pair<const uint8_t *, size_t>
StarDatabase::mmap_file(const std::string &path) {
  return mmap_file_ex(path, /*pin_and_prefault=*/true);
}

std::pair<const uint8_t *, size_t>
StarDatabase::mmap_file_ex(const std::string &path, bool pin_and_prefault) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw std::runtime_error("Could not open " + path + ": " +
                             std::strerror(errno));
  struct stat st {};
  if (::fstat(fd, &st) != 0) {
    int err = errno;
    ::close(fd);
    throw std::runtime_error("fstat failed on " + path + ": " +
                             std::strerror(err));
  }
  if (st.st_size <= 0) {
    ::close(fd);
    throw std::runtime_error("Empty or unreadable file: " + path);
  }
  size_t size = static_cast<size_t>(st.st_size);
  void *ptr = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
  // The fd can be closed immediately after mmap; the mapping holds its own
  // reference.
  ::close(fd);
  if (ptr == MAP_FAILED)
    throw std::runtime_error("mmap failed on " + path + ": " +
                             std::strerror(errno));
  // We binary-search the mapped region, so access is random, not sequential.
  // Advise the kernel accordingly so it doesn't waste read-ahead bandwidth.
  if (::madvise(ptr, size, MADV_RANDOM) != 0) {
    // Non-fatal: just a hint.
    std::cerr << "[mmap] madvise(MADV_RANDOM) failed on " << path << ": "
              << std::strerror(errno) << " (continuing)\n";
  }
  if (pin_and_prefault) {
    // Pin in physical memory so the first solve doesn't pay a page-in cost.
    // mlock may fail on Pi 4 without `ulimit -l unlimited` / a raised
    // RLIMIT_MEMLOCK; we treat that as a soft warning, not an error. The
    // mapping is still valid; pages just may be evicted under pressure.
    //
    // mlock is documented to populate the resident set ("All pages which
    // contain a part of the specified address range are guaranteed to be
    // resident in RAM when the call returns successfully"), so a successful
    // mlock implicitly prefaults. We still walk the pages on the mlock-failed
    // path so the first solve doesn't pay a major-fault avalanche on the cold
    // cache.
    bool mlocked = (::mlock(ptr, size) == 0);
    if (!mlocked) {
      std::cerr << "[mmap] mlock failed on " << path << ": "
                << std::strerror(errno)
                << " (RLIMIT_MEMLOCK too low?); continuing without pin\n";
      prefault(static_cast<const uint8_t *>(ptr), size);
    }
  }
  // pin_and_prefault==false: we deliberately leave pages unfaulted. Caller
  // owns the cost/latency tradeoff (used by the 198 MB partner index whose
  // access pattern is sparse — each find_partners() call touches one tiny
  // entry block, so prefaulting all 198 MB just to save a few hundred µs
  // of fault cost on first identify is a bad cold-start tradeoff).
  mappings_.push_back({ptr, size});
  return {static_cast<const uint8_t *>(ptr), size};
}

StarDatabase::~StarDatabase() {
  // Tear down in reverse insertion order purely for symmetry; munmap doesn't
  // care about ordering.
  for (auto it = mappings_.rbegin(); it != mappings_.rend(); ++it) {
    if (it->ptr && it->size > 0) {
      // munlock is best-effort: if we never successfully mlock'd (Pi 4
      // without RLIMIT_MEMLOCK), this returns EINVAL/EPERM. Either way the
      // subsequent munmap will release the mapping.
      (void)::munlock(it->ptr, it->size);
      (void)::munmap(it->ptr, it->size);
    }
  }
  mappings_.clear();
}

StarDatabase::StarDatabase(const std::string &star_file,
                           const std::string &pair_file) {
  // --- Stars: still ifstream. Tiny (~435 KB), feeds an unordered_map. ---
  // mmap would not help here: the destination is a hash table that must be
  // built in heap memory anyway, so any savings on the read are dwarfed by
  // the map insertions.
  std::ifstream fs_stars(star_file, std::ios::binary);
  if (!fs_stars)
    throw std::runtime_error("Could not open star file");

  int num_stars;
  fs_stars.read(reinterpret_cast<char *>(&num_stars), sizeof(int));

  for (int i = 0; i < num_stars; ++i) {
    int id;
    double x, y, z;
    fs_stars.read(reinterpret_cast<char *>(&id), sizeof(int));
    fs_stars.read(reinterpret_cast<char *>(&x), sizeof(double));
    fs_stars.read(reinterpret_cast<char *>(&y), sizeof(double));
    fs_stars.read(reinterpret_cast<char *>(&z), sizeof(double));
    star_map[id] = {id, x, y, z};
  }

  // --- Pairs: mmap. This is the 99 MB file that previously dominated
  // catalog_load (~360 ms of ifstream read+memcpy on M-series; much worse on
  // Pi 4 / SD card). ---
  auto [pair_bytes, pair_size] = mmap_file(pair_file);
  // Header: int32 num_pairs.
  if (pair_size < sizeof(int32_t))
    throw std::runtime_error("Pair file truncated (no header): " + pair_file);
  int32_t num_pairs = 0;
  std::memcpy(&num_pairs, pair_bytes, sizeof(int32_t));
  if (num_pairs < 0)
    throw std::runtime_error("Pair file claims negative num_pairs");
  size_t body_size = sizeof(int32_t) + static_cast<size_t>(num_pairs) *
                                           sizeof(CatalogPair);
  if (pair_size < body_size)
    throw std::runtime_error("Pair file truncated (header says " +
                             std::to_string(num_pairs) +
                             " pairs, file too short)");
  // CatalogPair is 16 bytes; the file header is 4 bytes, so the array starts
  // at byte 4. That's 4-byte aligned but the natural alignment of double is
  // 8 bytes. In practice both x86-64 and arm64 tolerate unaligned 8-byte
  // loads, and the generator writes the file with this exact 4-byte offset.
  // We reinterpret directly — the static_assert above guards struct size,
  // and on these targets the unaligned-double access is well-defined at the
  // hardware level (though strictly UB per the C++ standard; if a strict
  // architecture is added later, switch to memcpy-per-element).
  pairs_ptr_ =
      reinterpret_cast<const CatalogPair *>(pair_bytes + sizeof(int32_t));
  pairs_count_ = static_cast<size_t>(num_pairs);

  // --- Per-star partner index (Phase 3f.5): mmap if catalog_partners.bin
  // exists, otherwise fall back to the old in-memory build for old `data/`
  // directories. ---
  static_assert(sizeof(PartnerEntry) == 16,
                "PartnerEntry must be exactly 16 bytes (double + int32 + "
                "int32 pad) to match the on-disk layout written by "
                "tools/generate_catalog.py");
  std::string partners_file;
  {
    auto slash = pair_file.find_last_of("/\\");
    if (slash == std::string::npos)
      partners_file = "catalog_partners.bin";
    else
      partners_file = pair_file.substr(0, slash + 1) + "catalog_partners.bin";
  }
  bool partners_loaded = false;
  if (::access(partners_file.c_str(), R_OK) == 0) {
    try {
      // Use the no-prefault variant: the partners file is ~200 MB (~50k
      // pages) but find_partners only ever touches a tiny entry block per
      // call. Prefaulting all of it would add ~30 ms to catalog_load to save
      // a few hundred microseconds spread across the first solve — a bad
      // trade. We explicitly touch the header + directory below (which IS
      // walked at startup) so steady-state determinism on those is preserved.
      auto [pbytes, psize] = mmap_file_ex(partners_file,
                                          /*pin_and_prefault=*/false);
      // Header: int32 magic, int32 num_stars.
      constexpr int32_t kPartnersMagic = 0x50415254; // 'PART'
      const size_t hdr_size = 2 * sizeof(int32_t);
      if (psize < hdr_size)
        throw std::runtime_error("Partners file truncated (no header)");
      int32_t magic = 0, num_partner_stars = 0;
      std::memcpy(&magic, pbytes, sizeof(int32_t));
      std::memcpy(&num_partner_stars, pbytes + sizeof(int32_t),
                  sizeof(int32_t));
      if (magic != kPartnersMagic)
        throw std::runtime_error("Partners file bad magic");
      if (num_partner_stars < 0)
        throw std::runtime_error("Partners file negative num_stars");

      // Directory: num_partner_stars × (int32 hip, int32 count, int64 offset).
      constexpr size_t kDirEntrySize =
          sizeof(int32_t) + sizeof(int32_t) + sizeof(int64_t);
      const size_t dir_size =
          static_cast<size_t>(num_partner_stars) * kDirEntrySize;
      if (psize < hdr_size + dir_size)
        throw std::runtime_error("Partners file truncated (directory)");

      partners_index_.reserve(static_cast<size_t>(num_partner_stars));
      const uint8_t *dir_ptr = pbytes + hdr_size;
      for (int32_t i = 0; i < num_partner_stars; ++i) {
        int32_t hip = 0, count = 0;
        int64_t offset = 0;
        const uint8_t *rec = dir_ptr + i * kDirEntrySize;
        std::memcpy(&hip, rec, sizeof(int32_t));
        std::memcpy(&count, rec + sizeof(int32_t), sizeof(int32_t));
        std::memcpy(&offset, rec + 2 * sizeof(int32_t), sizeof(int64_t));
        if (count < 0 || offset < 0)
          throw std::runtime_error("Partners file negative count/offset");
        const size_t entry_bytes =
            static_cast<size_t>(count) * sizeof(PartnerEntry);
        if (static_cast<size_t>(offset) + entry_bytes > psize)
          throw std::runtime_error(
              "Partners file truncated (entry block past EOF for hip " +
              std::to_string(hip) + ")");
        const PartnerEntry *eptr = reinterpret_cast<const PartnerEntry *>(
            pbytes + static_cast<size_t>(offset));
        partners_index_.emplace(
            hip, std::make_pair(eptr, static_cast<size_t>(count)));
      }
      partners_loaded = true;
      std::cout << "Loaded per-star partner index: " << num_partner_stars
                << " stars (mmap'd " << partners_file << ")\n";
    } catch (const std::exception &e) {
      std::cerr << "Warning: failed to mmap partners file " << partners_file
                << ": " << e.what()
                << "; falling back to in-memory build.\n";
      partners_index_.clear();
      partners_loaded = false;
    }
  } else {
    std::cerr << "Warning: partners file " << partners_file
              << " not found; building per-star index in memory "
                 "(adds ~225ms to catalog_load). Regenerate the catalog with "
                 "tools/generate_catalog.py to eliminate this cost.\n";
  }
  if (!partners_loaded) {
    // Legacy build: bucket pairs by HIP, sort each bucket ascending by
    // cos_val. Same semantics as the mmap path; we just own the storage.
    partners_fallback_.reserve(num_stars);
    for (size_t i = 0; i < pairs_count_; ++i) {
      const CatalogPair &p = pairs_ptr_[i];
      partners_fallback_[p.id1].push_back({p.cos_val, p.id2, 0});
      partners_fallback_[p.id2].push_back({p.cos_val, p.id1, 0});
    }
    for (auto &kv : partners_fallback_) {
      std::sort(kv.second.begin(), kv.second.end(),
                [](const PartnerEntry &a, const PartnerEntry &b) {
                  if (a.cos_val != b.cos_val)
                    return a.cos_val < b.cos_val;
                  return a.partner_hip < b.partner_hip;
                });
      partners_index_.emplace(
          kv.first, std::make_pair(kv.second.data(), kv.second.size()));
    }
  }

  // --- Load Mortari k-vector index (optional, via mmap) ---
  // The index file lives alongside the pair file. If it is missing, the
  // database still works; find_pairs_kvec falls back to the binary-search
  // implementation.
  std::string kvec_file;
  {
    auto slash = pair_file.find_last_of("/\\");
    if (slash == std::string::npos)
      kvec_file = "catalog_kvec.bin";
    else
      kvec_file = pair_file.substr(0, slash + 1) + "catalog_kvec.bin";
  }
  // Probe for existence before going through mmap_file (which throws on
  // missing). The optional nature of kvec is preserved.
  if (::access(kvec_file.c_str(), R_OK) == 0) {
    auto [kvec_bytes, kvec_size] = mmap_file(kvec_file);
    // Header layout (must match tools/generate_catalog.py):
    //   int32 M, double y_min, double y_max, double dq, then (M+1) int32.
    const size_t hdr_size = sizeof(int32_t) + 3 * sizeof(double);
    if (kvec_size >= hdr_size) {
      int32_t M = 0;
      double y_min = 0, y_max = 0, dq = 0;
      std::memcpy(&M, kvec_bytes, sizeof(int32_t));
      std::memcpy(&y_min, kvec_bytes + sizeof(int32_t), sizeof(double));
      std::memcpy(&y_max, kvec_bytes + sizeof(int32_t) + sizeof(double),
                  sizeof(double));
      std::memcpy(&dq, kvec_bytes + sizeof(int32_t) + 2 * sizeof(double),
                  sizeof(double));
      const size_t body_bytes =
          hdr_size + static_cast<size_t>(M + 1) * sizeof(int32_t);
      if (M > 0 && kvec_size >= body_bytes) {
        kvec_K_ptr_ =
            reinterpret_cast<const int32_t *>(kvec_bytes + hdr_size);
        kvec_K_count_ = static_cast<size_t>(M) + 1;
        kvec_y_min = y_min;
        kvec_y_max = y_max;
        kvec_dq = dq;
        kvec_M = M;
        std::cout << "Loaded k-vector index: " << kvec_M
                  << " bins, dq=" << kvec_dq << "\n";
      } else {
        std::cerr << "Warning: k-vector file " << kvec_file
                  << " truncated or invalid M; falling back to binary "
                     "search.\n";
      }
    } else {
      std::cerr << "Warning: k-vector file " << kvec_file
                << " has invalid header; falling back to binary search.\n";
    }
  } else {
    std::cerr << "Warning: k-vector file " << kvec_file
              << " not found; find_pairs_kvec will fall back to find_pairs.\n";
  }

  std::cout << "Loaded database: " << num_stars << " stars, " << num_pairs
            << " pairs.\n";
}

std::vector<int> StarDatabase::find_partners(int hip, double cos_target,
                                             double cos_tolerance) const {
  std::vector<int> result;
  auto it = partners_index_.find(hip);
  if (it == partners_index_.end())
    return result;
  const PartnerEntry *begin = it->second.first;
  const PartnerEntry *end = begin + it->second.second;
  // Sorted ascending by cos_val. Find range [cos_target - tol,
  // cos_target + tol]. Compare against cos_val only — exactly matches the
  // previous std::pair<double, int> behavior when we used the smallest /
  // largest int sentinel for the second field.
  const double cos_low = cos_target - cos_tolerance;
  const double cos_high = cos_target + cos_tolerance;
  auto lo = std::lower_bound(
      begin, end, cos_low,
      [](const PartnerEntry &e, double v) { return e.cos_val < v; });
  auto hi = std::upper_bound(
      lo, end, cos_high,
      [](double v, const PartnerEntry &e) { return v < e.cos_val; });
  result.reserve(static_cast<size_t>(hi - lo));
  for (auto p = lo; p != hi; ++p)
    result.push_back(p->partner_hip);
  return result;
}

std::vector<CatalogPair> StarDatabase::find_pairs(double cos_target,
                                                  double cos_tolerance) const {
  // Array is sorted descending by cos_val (meaning closest pairs first)
  // We want cos_val in [cos_target - cos_tolerance, cos_target + cos_tolerance]

  // upper_bound/lower_bound requires binary search logic adapted for descending
  // order
  auto comp = [](const CatalogPair &p, double val) {
    return p.cos_val > val; // descending
  };

  const CatalogPair *begin = pairs_ptr_;
  const CatalogPair *end = pairs_ptr_ + pairs_count_;
  auto it_begin =
      std::lower_bound(begin, end, cos_target + cos_tolerance, comp);
  auto it_end = std::lower_bound(begin, end, cos_target - cos_tolerance, comp);

  return std::vector<CatalogPair>(it_begin, it_end);
}

std::vector<CatalogPair>
StarDatabase::find_pairs_kvec(double cos_target,
                              double cos_tolerance) const {
  // Fall back to binary search when the k-vector index is unavailable.
  if (kvec_K_ptr_ == nullptr || kvec_M <= 0 || kvec_dq <= 0.0 ||
      pairs_count_ == 0)
    return find_pairs(cos_target, cos_tolerance);

  const double cos_low = cos_target - cos_tolerance;
  const double cos_high = cos_target + cos_tolerance;

  // Bin lookup: i_low = floor((cos_low - y_min)/dq), i_high = ceil((cos_high
  // - y_min)/dq), both clamped to [0, M].
  auto clamp_bin = [this](double f) -> int {
    if (!(f == f)) // NaN guard
      return 0;
    if (f < 0.0)
      return 0;
    if (f > static_cast<double>(kvec_M))
      return kvec_M;
    return static_cast<int>(f);
  };

  const double f_low = (cos_low - kvec_y_min) / kvec_dq;
  const double f_high = (cos_high - kvec_y_min) / kvec_dq;
  const int i_low = clamp_bin(std::floor(f_low));
  const int i_high = clamp_bin(std::ceil(f_high));

  // kvec_K[i] is the largest index j (in the descending pairs array) such
  // that pairs[j].cos_val >= y_min + i*dq. Pairs are descending in cos_val:
  //   - All pairs at index <= K[i_high] have cos_val >= y_min + i_high*dq.
  //     Because i_high = ceil((cos_high - y_min)/dq), we have
  //     y_min + i_high*dq >= cos_high. So pairs at index <= K[i_high] could
  //     have cos_val > cos_high (out of range high); the *first* candidate
  //     that might satisfy cos_val <= cos_high lies at index K[i_high] (we
  //     include it for safety) or later.
  //   - All pairs at index > K[i_low] have cos_val < y_min + i_low*dq.
  //     Because i_low = floor((cos_low - y_min)/dq), we have
  //     y_min + i_low*dq <= cos_low. So pairs beyond K[i_low] are below
  //     cos_low (out of range low). The last candidate is at index K[i_low].
  //
  // Inclusive candidate range: [start, end] = [K[i_high], K[i_low]].
  // Since dq is chosen << cos_tolerance, the slack at each edge is at most
  // one or two entries. We trim with a tiny linear walk so the final return
  // can be a single std::vector range copy (memcpy), matching find_pairs'
  // throughput.
  int start = kvec_K_ptr_[i_high];
  int end = kvec_K_ptr_[i_low];
  if (start < 0)
    start = 0;
  const int P = static_cast<int>(pairs_count_);
  if (end >= P)
    end = P - 1;
  if (end < 0 || start > end || start >= P)
    return {};

  // Trim the high-cos slack from the front: advance `start` while
  // pairs[start].cos_val > cos_high. Bounded by ceil(dq * (M+1) / dq) = O(1)
  // entries in practice because dq << tolerance.
  while (start <= end && pairs_ptr_[start].cos_val > cos_high)
    ++start;
  // Trim the low-cos slack from the back: retreat `end` while
  // pairs[end].cos_val < cos_low.
  while (end >= start && pairs_ptr_[end].cos_val < cos_low)
    --end;
  if (start > end)
    return {};

  return std::vector<CatalogPair>(pairs_ptr_ + start, pairs_ptr_ + end + 1);
}

CatalogStar StarDatabase::get_star(int hip_id) const {
  auto it = star_map.find(hip_id);
  if (it != star_map.end()) {
    return it->second;
  }
  throw std::runtime_error("Star ID not found");
}

// --- Phase 3e.2: pattern-hash catalog loader + query ---

namespace {
// Must match tools/generate_catalog.py: PATTERN_MAGIC = 0x50415431 ('PAT1').
constexpr int32_t kPatternMagic = 0x50415431;
constexpr int kQuantBits = 10;
constexpr uint64_t kQuantMask = (1ULL << kQuantBits) - 1ULL;
constexpr uint64_t kQuantMax = kQuantMask; // 1023
} // namespace

void StarDatabase::load_pattern_catalog(const std::string &pattern_file) {
  // mmap-based load (Phase 3f.4). 21 MB file; previous ifstream read took
  // tens of ms on cold cache. Header is 5 × int32 (20 bytes), then a
  // densely-packed array of StarPattern (24 bytes each).
  auto [bytes, size] = mmap_file(pattern_file);
  const size_t hdr_size = 5 * sizeof(int32_t);
  if (size < hdr_size)
    throw std::runtime_error(
        "Pattern catalog truncated (no header): " + pattern_file);

  int32_t magic = 0, fov_bin_deg = 0, k_nearest = 0, quant_bits = 0,
          num_patterns = 0;
  std::memcpy(&magic, bytes + 0 * sizeof(int32_t), sizeof(int32_t));
  std::memcpy(&fov_bin_deg, bytes + 1 * sizeof(int32_t), sizeof(int32_t));
  std::memcpy(&k_nearest, bytes + 2 * sizeof(int32_t), sizeof(int32_t));
  std::memcpy(&quant_bits, bytes + 3 * sizeof(int32_t), sizeof(int32_t));
  std::memcpy(&num_patterns, bytes + 4 * sizeof(int32_t), sizeof(int32_t));

  if (magic != kPatternMagic)
    throw std::runtime_error("Pattern catalog bad magic: " + pattern_file);
  if (quant_bits != kQuantBits)
    throw std::runtime_error("Pattern catalog quant_bits mismatch (expected " +
                             std::to_string(kQuantBits) + ", got " +
                             std::to_string(quant_bits) + ")");
  if (num_patterns < 0)
    throw std::runtime_error("Pattern catalog negative num_patterns");

  const size_t body_bytes =
      hdr_size + static_cast<size_t>(num_patterns) * sizeof(StarPattern);
  if (size < body_bytes)
    throw std::runtime_error("Pattern catalog truncated (header says " +
                             std::to_string(num_patterns) +
                             " patterns, file too short): " + pattern_file);

  pattern_fov_bin_deg_ = fov_bin_deg;
  pattern_k_nearest_ = k_nearest;
  pattern_quant_bits_ = quant_bits;

  // 20-byte header, 24-byte StarPattern → array starts at byte 20 which is
  // 4-byte aligned but not 8-byte aligned. uint64_t loads on the targets we
  // support (x86-64, arm64) tolerate this; if a strict-alignment target is
  // added, replace with memcpy-per-record.
  patterns_ptr_ = reinterpret_cast<const StarPattern *>(bytes + hdr_size);
  patterns_count_ = static_cast<size_t>(num_patterns);

  // Generator already sorts ascending by key, but verify cheaply to fail
  // loudly on a corrupted/mis-ordered file rather than silently returning
  // wrong results from std::lower_bound.
  if (!std::is_sorted(patterns_ptr_, patterns_ptr_ + patterns_count_,
                      [](const StarPattern &a, const StarPattern &b) {
                        return a.key < b.key;
                      })) {
    throw std::runtime_error(
        "Pattern catalog not sorted ascending by key: " + pattern_file);
  }

  std::cout << "Loaded pattern catalog: " << pattern_file << " ("
            << num_patterns << " patterns, fov_bin=" << fov_bin_deg
            << "°, k=" << k_nearest << ")\n";
}

std::vector<StarPattern> StarDatabase::find_pattern(uint64_t key) const {
  std::vector<StarPattern> out;
  if (patterns_count_ == 0)
    return out;

  // The vector is sorted ascending by key. equal_range gives the half-open
  // [lo, hi) range whose keys equal `key`.
  auto cmp_key = [](const StarPattern &p, uint64_t k) { return p.key < k; };
  auto cmp_key_rev = [](uint64_t k, const StarPattern &p) { return k < p.key; };
  const StarPattern *begin = patterns_ptr_;
  const StarPattern *end = patterns_ptr_ + patterns_count_;
  auto lo = std::lower_bound(begin, end, key, cmp_key);
  auto hi = std::upper_bound(lo, end, key, cmp_key_rev);
  out.reserve(static_cast<size_t>(hi - lo));
  for (auto it = lo; it != hi; ++it)
    out.push_back(*it);
  return out;
}

std::vector<StarPattern>
StarDatabase::find_pattern_tolerant(uint64_t key) const {
  std::vector<StarPattern> out;
  if (patterns_count_ == 0)
    return out;

  // Unpack key into 5 × 10-bit components.
  std::array<int, 5> q = {
      static_cast<int>((key >> 0) & kQuantMask),
      static_cast<int>((key >> 10) & kQuantMask),
      static_cast<int>((key >> 20) & kQuantMask),
      static_cast<int>((key >> 30) & kQuantMask),
      static_cast<int>((key >> 40) & kQuantMask),
  };

  const StarPattern *begin = patterns_ptr_;
  const StarPattern *end = patterns_ptr_ + patterns_count_;

  // 3^5 = 243 probes. Skip probes whose perturbed value would underflow
  // below 0 or overflow above 1023. The exact-key probe (offsets all 0)
  // is included, so this strictly returns a superset of find_pattern(key).
  std::array<int, 5> off{};
  for (off[0] = -1; off[0] <= 1; ++off[0]) {
    int v0 = q[0] + off[0];
    if (v0 < 0 || v0 > static_cast<int>(kQuantMax))
      continue;
    for (off[1] = -1; off[1] <= 1; ++off[1]) {
      int v1 = q[1] + off[1];
      if (v1 < 0 || v1 > static_cast<int>(kQuantMax))
        continue;
      for (off[2] = -1; off[2] <= 1; ++off[2]) {
        int v2 = q[2] + off[2];
        if (v2 < 0 || v2 > static_cast<int>(kQuantMax))
          continue;
        for (off[3] = -1; off[3] <= 1; ++off[3]) {
          int v3 = q[3] + off[3];
          if (v3 < 0 || v3 > static_cast<int>(kQuantMax))
            continue;
          for (off[4] = -1; off[4] <= 1; ++off[4]) {
            int v4 = q[4] + off[4];
            if (v4 < 0 || v4 > static_cast<int>(kQuantMax))
              continue;
            uint64_t probe = (static_cast<uint64_t>(v4) << 40) |
                             (static_cast<uint64_t>(v3) << 30) |
                             (static_cast<uint64_t>(v2) << 20) |
                             (static_cast<uint64_t>(v1) << 10) |
                             (static_cast<uint64_t>(v0));
            auto cmp_key = [](const StarPattern &p, uint64_t k) {
              return p.key < k;
            };
            auto cmp_key_rev = [](uint64_t k, const StarPattern &p) {
              return k < p.key;
            };
            auto lo = std::lower_bound(begin, end, probe, cmp_key);
            auto hi = std::upper_bound(lo, end, probe, cmp_key_rev);
            for (auto it = lo; it != hi; ++it)
              out.push_back(*it);
          }
        }
      }
    }
  }
  return out;
}
