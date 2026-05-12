#include "catalog.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

StarDatabase::StarDatabase(const std::string &star_file,
                           const std::string &pair_file) {
  // Read stars
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

  // Read pairs
  std::ifstream fs_pairs(pair_file, std::ios::binary);
  if (!fs_pairs)
    throw std::runtime_error("Could not open pair file");

  int num_pairs;
  fs_pairs.read(reinterpret_cast<char *>(&num_pairs), sizeof(int));

  pairs.resize(num_pairs);
  for (int i = 0; i < num_pairs; ++i) {
    double cos_val;
    int id1, id2;
    fs_pairs.read(reinterpret_cast<char *>(&cos_val), sizeof(double));
    fs_pairs.read(reinterpret_cast<char *>(&id1), sizeof(int));
    fs_pairs.read(reinterpret_cast<char *>(&id2), sizeof(int));
    pairs[i] = {cos_val, id1, id2};
  }

  // Build per-star partner index for pyramid-style identification expansion.
  per_star_partners.reserve(num_stars);
  for (const auto &p : pairs) {
    per_star_partners[p.id1].emplace_back(p.cos_val, p.id2);
    per_star_partners[p.id2].emplace_back(p.cos_val, p.id1);
  }
  for (auto &kv : per_star_partners) {
    std::sort(kv.second.begin(), kv.second.end());
  }

  // --- Load Mortari k-vector index (optional) ---
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
  std::ifstream fs_kvec(kvec_file, std::ios::binary);
  if (fs_kvec) {
    int M = 0;
    fs_kvec.read(reinterpret_cast<char *>(&M), sizeof(int));
    double y_min = 0.0, y_max = 0.0, dq = 0.0;
    fs_kvec.read(reinterpret_cast<char *>(&y_min), sizeof(double));
    fs_kvec.read(reinterpret_cast<char *>(&y_max), sizeof(double));
    fs_kvec.read(reinterpret_cast<char *>(&dq), sizeof(double));
    if (M > 0 && fs_kvec) {
      std::vector<int> K(static_cast<size_t>(M) + 1);
      fs_kvec.read(reinterpret_cast<char *>(K.data()),
                   static_cast<std::streamsize>((static_cast<size_t>(M) + 1) *
                                                sizeof(int)));
      if (fs_kvec) {
        kvec_K = std::move(K);
        kvec_y_min = y_min;
        kvec_y_max = y_max;
        kvec_dq = dq;
        kvec_M = M;
        std::cout << "Loaded k-vector index: " << kvec_M
                  << " bins, dq=" << kvec_dq << "\n";
      } else {
        std::cerr << "Warning: k-vector file " << kvec_file
                  << " truncated; falling back to binary search.\n";
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
  auto it = per_star_partners.find(hip);
  if (it == per_star_partners.end())
    return result;
  const auto &lst = it->second;
  // Sorted ascending by cos. Find range [cos_target - tol, cos_target + tol].
  auto lo = std::lower_bound(
      lst.begin(), lst.end(),
      std::make_pair(cos_target - cos_tolerance, std::numeric_limits<int>::min()));
  auto hi = std::upper_bound(
      lst.begin(), lst.end(),
      std::make_pair(cos_target + cos_tolerance, std::numeric_limits<int>::max()));
  result.reserve(static_cast<size_t>(hi - lo));
  for (auto p = lo; p != hi; ++p)
    result.push_back(p->second);
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

  auto it_begin = std::lower_bound(pairs.begin(), pairs.end(),
                                   cos_target + cos_tolerance, comp);
  auto it_end = std::lower_bound(pairs.begin(), pairs.end(),
                                 cos_target - cos_tolerance, comp);

  return std::vector<CatalogPair>(it_begin, it_end);
}

std::vector<CatalogPair>
StarDatabase::find_pairs_kvec(double cos_target,
                              double cos_tolerance) const {
  // Fall back to binary search when the k-vector index is unavailable.
  if (kvec_K.empty() || kvec_M <= 0 || kvec_dq <= 0.0 || pairs.empty())
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
  int start = kvec_K[i_high];
  int end = kvec_K[i_low];
  if (start < 0)
    start = 0;
  const int P = static_cast<int>(pairs.size());
  if (end >= P)
    end = P - 1;
  if (end < 0 || start > end || start >= P)
    return {};

  // Trim the high-cos slack from the front: advance `start` while
  // pairs[start].cos_val > cos_high. Bounded by ceil(dq * (M+1) / dq) = O(1)
  // entries in practice because dq << tolerance.
  while (start <= end && pairs[start].cos_val > cos_high)
    ++start;
  // Trim the low-cos slack from the back: retreat `end` while
  // pairs[end].cos_val < cos_low.
  while (end >= start && pairs[end].cos_val < cos_low)
    --end;
  if (start > end)
    return {};

  return std::vector<CatalogPair>(pairs.begin() + start,
                                  pairs.begin() + end + 1);
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
  std::ifstream fs(pattern_file, std::ios::binary);
  if (!fs)
    throw std::runtime_error("Could not open pattern catalog: " + pattern_file);

  int32_t magic = 0;
  int32_t fov_bin_deg = 0;
  int32_t k_nearest = 0;
  int32_t quant_bits = 0;
  int32_t num_patterns = 0;
  fs.read(reinterpret_cast<char *>(&magic), sizeof(int32_t));
  fs.read(reinterpret_cast<char *>(&fov_bin_deg), sizeof(int32_t));
  fs.read(reinterpret_cast<char *>(&k_nearest), sizeof(int32_t));
  fs.read(reinterpret_cast<char *>(&quant_bits), sizeof(int32_t));
  fs.read(reinterpret_cast<char *>(&num_patterns), sizeof(int32_t));
  if (!fs)
    throw std::runtime_error("Pattern catalog header read failed: " +
                             pattern_file);
  if (magic != kPatternMagic)
    throw std::runtime_error("Pattern catalog bad magic: " + pattern_file);
  if (quant_bits != kQuantBits)
    throw std::runtime_error("Pattern catalog quant_bits mismatch (expected " +
                             std::to_string(kQuantBits) + ", got " +
                             std::to_string(quant_bits) + ")");
  if (num_patterns < 0)
    throw std::runtime_error("Pattern catalog negative num_patterns");

  pattern_fov_bin_deg_ = fov_bin_deg;
  pattern_k_nearest_ = k_nearest;
  pattern_quant_bits_ = quant_bits;

  patterns_.clear();
  patterns_.resize(static_cast<size_t>(num_patterns));
  for (int i = 0; i < num_patterns; ++i) {
    StarPattern p;
    fs.read(reinterpret_cast<char *>(&p.key), sizeof(uint64_t));
    fs.read(reinterpret_cast<char *>(&p.hips[0]), sizeof(int32_t) * 4);
    if (!fs)
      throw std::runtime_error(
          "Pattern catalog body short at index " + std::to_string(i) + " of " +
          std::to_string(num_patterns));
    patterns_[i] = p;
  }

  // Generator already sorts ascending by key, but verify cheaply to fail
  // loudly on a corrupted/mis-ordered file rather than silently returning
  // wrong results from std::lower_bound.
  if (!std::is_sorted(patterns_.begin(), patterns_.end(),
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
  if (patterns_.empty())
    return out;

  // The vector is sorted ascending by key. equal_range gives the half-open
  // [lo, hi) range whose keys equal `key`.
  auto cmp_key = [](const StarPattern &p, uint64_t k) { return p.key < k; };
  auto cmp_key_rev = [](uint64_t k, const StarPattern &p) { return k < p.key; };
  auto lo = std::lower_bound(patterns_.begin(), patterns_.end(), key, cmp_key);
  auto hi = std::upper_bound(lo, patterns_.end(), key, cmp_key_rev);
  out.reserve(static_cast<size_t>(hi - lo));
  for (auto it = lo; it != hi; ++it)
    out.push_back(*it);
  return out;
}

std::vector<StarPattern>
StarDatabase::find_pattern_tolerant(uint64_t key) const {
  std::vector<StarPattern> out;
  if (patterns_.empty())
    return out;

  // Unpack key into 5 × 10-bit components.
  std::array<int, 5> q = {
      static_cast<int>((key >> 0) & kQuantMask),
      static_cast<int>((key >> 10) & kQuantMask),
      static_cast<int>((key >> 20) & kQuantMask),
      static_cast<int>((key >> 30) & kQuantMask),
      static_cast<int>((key >> 40) & kQuantMask),
  };

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
            auto lo = std::lower_bound(patterns_.begin(), patterns_.end(),
                                       probe, cmp_key);
            auto hi =
                std::upper_bound(lo, patterns_.end(), probe, cmp_key_rev);
            for (auto it = lo; it != hi; ++it)
              out.push_back(*it);
          }
        }
      }
    }
  }
  return out;
}
