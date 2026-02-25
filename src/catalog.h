#pragma once
#include <string>
#include <unordered_map>
#include <vector>

struct CatalogStar {
  int hip_id;
  double x, y, z;
};

struct CatalogPair {
  double cos_val;
  int id1, id2;
};

class StarDatabase {
public:
  StarDatabase(const std::string &star_file, const std::string &pair_file);

  // Find all catalog pairs within a tolerance around a cosine distance
  std::vector<CatalogPair> find_pairs(double cos_target,
                                      double cos_tolerance) const;

  // Retrieve star by HIP ID
  CatalogStar get_star(int hip_id) const;

private:
  std::unordered_map<int, CatalogStar> star_map;
  std::vector<CatalogPair> pairs; // Sorted descending by cos_val
};
