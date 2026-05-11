#include "identification.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance) {
  int N = image_stars.size();
  if (N < 3)
    return {};

  // === Camera-frame unit vectors ===
  std::vector<std::array<double, 3>> v_cam(N);
  for (int i = 0; i < N; ++i) {
    double vx = (image_stars[i].x - camera.center_x) / camera.focal_x;
    double vy = (image_stars[i].y - camera.center_y) / camera.focal_y;
    double vz = 1.0;
    double n = std::sqrt(vx * vx + vy * vy + vz * vz);
    v_cam[i] = {vx / n, vy / n, vz / n};
  }

  auto dot3 = [](const std::array<double, 3> &a,
                 const std::array<double, 3> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };

  // === Pyramid (geometric voting) identification ===
  //
  // Vote-based prefiltering fails on dense catalogs: random false matches in
  // dense sky regions accumulate enough votes to bury the correct HIPs.
  // Instead, search the catalog pair database directly:
  //   1. For each centroid pair (i, j), find catalog pairs (A, B) whose
  //      cosine distance matches the observed pair angle (db.find_pairs).
  //   2. For each candidate seed (A, B), try i=A,j=B AND i=B,j=A.
  //   3. Expand each seed by finding, for every other centroid k, a HIP C
  //      such that cos(A,C) matches obs(i,k) AND cos(B,C) matches obs(j,k).
  //      Use db.find_partners for fast O(log P) lookup of candidate C.
  //   4. Keep the seed that produces the most inliers — that's the globally
  //      self-consistent assignment.

  // Cache catalog unit vectors (one per HIP) to avoid repeated map lookups.
  std::unordered_map<int, std::array<double, 3>> cat_vec_cache;
  auto cat_vec = [&](int hip) -> const std::array<double, 3> & {
    auto it = cat_vec_cache.find(hip);
    if (it != cat_vec_cache.end())
      return it->second;
    CatalogStar s = db.get_star(hip);
    return cat_vec_cache[hip] = {s.x, s.y, s.z};
  };

  std::vector<int> best_assignment(N, -1);
  int best_inliers = 0;

  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      double obs_ij = dot3(v_cam[i], v_cam[j]);
      auto seeds = db.find_pairs(obs_ij, cos_tolerance);
      for (const auto &seed : seeds) {
        // Try both orientations: i=id1/j=id2 AND i=id2/j=id1.
        for (int swap = 0; swap < 2; ++swap) {
          int A = swap ? seed.id2 : seed.id1;
          int B = swap ? seed.id1 : seed.id2;
          if (A == B)
            continue;

          const auto &va = cat_vec(A);
          const auto &vb = cat_vec(B);

          std::vector<int> assign(N, -1);
          std::unordered_set<int> used;
          assign[i] = A;
          assign[j] = B;
          used.insert(A);
          used.insert(B);
          int inliers = 2;

          for (int k = 0; k < N; ++k) {
            if (k == i || k == j)
              continue;
            double obs_ik = dot3(v_cam[i], v_cam[k]);
            double obs_jk = dot3(v_cam[j], v_cam[k]);
            // Candidates C such that cos(A, C) ≈ obs_ik.
            auto candidates_C = db.find_partners(A, obs_ik, cos_tolerance);
            int best_C = -1;
            double best_err = cos_tolerance;
            for (int C : candidates_C) {
              if (used.count(C))
                continue;
              const auto &vc = cat_vec(C);
              double cat_bc = dot3(vb, vc);
              double err = std::abs(cat_bc - obs_jk);
              if (err < best_err) {
                best_err = err;
                best_C = C;
              }
            }
            if (best_C != -1) {
              assign[k] = best_C;
              used.insert(best_C);
              inliers++;
            }
          }

          if (inliers > best_inliers) {
            best_inliers = inliers;
            best_assignment = std::move(assign);
          }
        }
      }
    }
  }

  // Require enough inliers to outrank random coincidence.
  const int min_inliers = std::min(N, std::max(4, N / 4));
  if (best_inliers < min_inliers)
    return {};

  std::vector<IdentifiedStar> identified;
  identified.reserve(best_inliers);
  for (int i = 0; i < N; ++i) {
    if (best_assignment[i] == -1)
      continue;
    IdentifiedStar is;
    is.image_idx = i;
    is.catalog_hip_id = best_assignment[i];
    is.v_cam[0] = v_cam[i][0];
    is.v_cam[1] = v_cam[i][1];
    is.v_cam[2] = v_cam[i][2];
    identified.push_back(is);
  }
  return identified;
}
