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
  // (binary-search implementation, kept for parity tests).
  std::vector<CatalogPair> find_pairs(double cos_target,
                                      double cos_tolerance) const;

  // Same semantics as find_pairs, but uses Mortari's k-vector for O(1)
  // boundary lookup. Falls back to find_pairs if the k-vector index is not
  // loaded.
  std::vector<CatalogPair> find_pairs_kvec(double cos_target,
                                           double cos_tolerance) const;

  // Find all HIPs C such that cos(hip, C) is within tolerance of cos_target.
  // Used by pyramid-style identification to extend a seed (A, B) to a third
  // star without scanning the full pair list.
  std::vector<int> find_partners(int hip, double cos_target,
                                 double cos_tolerance) const;

  // Retrieve star by HIP ID
  CatalogStar get_star(int hip_id) const;

private:
  std::unordered_map<int, CatalogStar> star_map;
  std::vector<CatalogPair> pairs; // Sorted descending by cos_val

  // Per-star index: hip -> sorted (cos_distance, partner_hip) entries.
  // Built once in the constructor from `pairs`. Memory: 2 * |pairs| entries.
  std::unordered_map<int, std::vector<std::pair<double, int>>>
      per_star_partners;

  // Mortari k-vector index over `pairs`. `kvec_K[i]` is the largest index j in
  // the descending pairs array such that pairs[j].cos_val >= y_min + i * dq.
  // Empty if the index file was not present (find_pairs_kvec then falls back).
  std::vector<int> kvec_K;
  double kvec_y_min = 0.0;
  double kvec_y_max = 0.0;
  double kvec_dq = 0.0;
  int kvec_M = 0;
};
