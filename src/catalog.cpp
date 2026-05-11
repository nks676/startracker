#include "catalog.h"
#include <algorithm>
#include <cmath>
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
