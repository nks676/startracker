#include "identification.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

void project(const CameraModel &cam, const double v[3], double &px,
             double &py) {
  // Pinhole step: normalise to z=1 plane.
  const double x = v[0] / v[2];
  const double y = v[1] / v[2];

  const double r2 = x * x + y * y;
  const double radial =
      1.0 + cam.k1 * r2 + cam.k2 * r2 * r2 + cam.k3 * r2 * r2 * r2;
  const double x_d =
      x * radial + 2.0 * cam.p1 * x * y + cam.p2 * (r2 + 2.0 * x * x);
  const double y_d =
      y * radial + cam.p1 * (r2 + 2.0 * y * y) + 2.0 * cam.p2 * x * y;

  px = cam.focal_x * x_d + cam.center_x;
  py = cam.focal_y * y_d + cam.center_y;
}

void undistort_to_unit_vector(const CameraModel &cam, double px, double py,
                              double v_out[3]) {
  // Normalised distorted pixel coordinates.
  const double px_norm = (px - cam.center_x) / cam.focal_x;
  const double py_norm = (py - cam.center_y) / cam.focal_y;

  // Initial guess: assume zero distortion. With k1..p2 = 0 this is the exact
  // answer and the loop body is a no-op (radial == 1, dx = dy = 0), so the
  // result is bitwise identical to the legacy pinhole code path.
  double x = px_norm;
  double y = py_norm;

  for (int iter = 0; iter < 10; ++iter) {
    const double r2 = x * x + y * y;
    const double radial =
        1.0 + cam.k1 * r2 + cam.k2 * r2 * r2 + cam.k3 * r2 * r2 * r2;
    const double dx =
        2.0 * cam.p1 * x * y + cam.p2 * (r2 + 2.0 * x * x);
    const double dy =
        cam.p1 * (r2 + 2.0 * y * y) + 2.0 * cam.p2 * x * y;
    x = (px_norm - dx) / radial;
    y = (py_norm - dy) / radial;
  }

  const double n = std::sqrt(x * x + y * y + 1.0);
  v_out[0] = x / n;
  v_out[1] = y / n;
  v_out[2] = 1.0 / n;
}

std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance) {
  int N = image_stars.size();
  if (N < 3)
    return {};

  // === Camera-frame unit vectors ===
  // Routes through the camera model so Brown-Conrady distortion (if any) is
  // inverted before identification. With default zero coefficients this
  // reduces to the original pinhole math.
  std::vector<std::array<double, 3>> v_cam(N);
  for (int i = 0; i < N; ++i) {
    double v[3];
    undistort_to_unit_vector(camera, image_stars[i].x, image_stars[i].y, v);
    v_cam[i] = {v[0], v[1], v[2]};
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
  // Returns by value: holding a reference into an unordered_map across a
  // subsequent insert would be UB on rehash.
  std::unordered_map<int, std::array<double, 3>> cat_vec_cache;
  auto cat_vec = [&](int hip) -> std::array<double, 3> {
    auto it = cat_vec_cache.find(hip);
    if (it != cat_vec_cache.end())
      return it->second;
    CatalogStar s = db.get_star(hip);
    std::array<double, 3> v = {s.x, s.y, s.z};
    cat_vec_cache[hip] = v;
    return v;
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

          const auto va = cat_vec(A);
          const auto vb = cat_vec(B);

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
              const auto vc = cat_vec(C);
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
