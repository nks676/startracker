#include "catalog.h"
#include <algorithm>
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

CatalogStar StarDatabase::get_star(int hip_id) const {
  auto it = star_map.find(hip_id);
  if (it != star_map.end()) {
    return it->second;
  }
  throw std::runtime_error("Star ID not found");
}
