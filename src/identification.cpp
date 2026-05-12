#include "identification.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

namespace {

inline double dot3(const std::array<double, 3> &a,
                   const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute camera-frame unit vectors from pixel centroids by inverting the
// camera distortion model.
std::vector<std::array<double, 3>>
compute_v_cam(const std::vector<StarCentroid> &image_stars,
              const CameraModel &camera) {
  std::vector<std::array<double, 3>> v_cam(image_stars.size());
  for (size_t i = 0; i < image_stars.size(); ++i) {
    double v[3];
    undistort_to_unit_vector(camera, image_stars[i].x, image_stars[i].y, v);
    v_cam[i] = {v[0], v[1], v[2]};
  }
  return v_cam;
}

// Pyramid (geometric voting) identification core.
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
//
// Returns the best assignment (image index -> HIP, or -1) and inlier count.
std::pair<std::vector<int>, int>
run_pyramid(const std::vector<std::array<double, 3>> &v_cam,
            const StarDatabase &db, double cos_tolerance,
            const std::function<std::array<double, 3>(int)> &cat_vec) {
  const int N = static_cast<int>(v_cam.size());
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

  return {std::move(best_assignment), best_inliers};
}

// Estimate the FOV scale-factor s = (true focal) / (assumed focal) from a
// coarse identification result. For every pair of assigned image stars we
// compare the observed inter-star angle (from v_cam, which was computed with
// the assumed focal) against the catalog angle, and take the median of the
// ratios. If the assumed focal is too small, observed unit vectors splay
// outward, so observed_angle > catalog_angle and s > 1.
double estimate_scale_factor(
    const std::vector<int> &assignment,
    const std::vector<std::array<double, 3>> &v_cam,
    const StarDatabase &db,
    const std::function<std::array<double, 3>(int)> &cat_vec) {
  (void)db; // catalog vectors come through `cat_vec`; reserved for future use
  std::vector<double> ratios;
  const int N = static_cast<int>(assignment.size());
  ratios.reserve(static_cast<size_t>(N) * static_cast<size_t>(N) / 2);
  for (int i = 0; i < N; ++i) {
    if (assignment[i] < 0)
      continue;
    const auto vi = cat_vec(assignment[i]);
    for (int j = i + 1; j < N; ++j) {
      if (assignment[j] < 0)
        continue;
      double obs = std::clamp(dot3(v_cam[i], v_cam[j]), -1.0, 1.0);
      const auto vj = cat_vec(assignment[j]);
      double cat = std::clamp(dot3(vi, vj), -1.0, 1.0);
      double angle_obs = std::acos(obs);
      double angle_cat = std::acos(cat);
      // Skip degenerate (near-zero) catalog angles — division blows up.
      if (angle_cat < 1e-9)
        continue;
      ratios.push_back(angle_obs / angle_cat);
    }
  }
  if (ratios.empty())
    return 1.0;
  size_t mid = ratios.size() / 2;
  std::nth_element(ratios.begin(), ratios.begin() + mid, ratios.end());
  return ratios[mid];
}

// Cross-verification: drop any assigned star whose median pairwise residual
// against the other assigned stars exceeds 3 * cos_tolerance. Repeats until
// stable (max 3 iterations). Mutates `assignment` and returns the new inlier
// count.
int cross_verify(std::vector<int> &assignment,
                 const std::vector<std::array<double, 3>> &v_cam,
                 double cos_tolerance,
                 const std::function<std::array<double, 3>(int)> &cat_vec) {
  const int N = static_cast<int>(assignment.size());
  const double cutoff = 3.0 * cos_tolerance;

  int inliers = 0;
  for (int i = 0; i < N; ++i)
    if (assignment[i] >= 0)
      ++inliers;

  for (int pass = 0; pass < 3; ++pass) {
    bool changed = false;

    // Snapshot of currently-assigned image indices for this pass.
    std::vector<int> assigned_idx;
    assigned_idx.reserve(static_cast<size_t>(inliers));
    for (int i = 0; i < N; ++i)
      if (assignment[i] >= 0)
        assigned_idx.push_back(i);

    if (assigned_idx.size() < 3) // need at least 2 partners to form a residual
      break;

    // For each star, compute its median |obs - cat| against the others.
    std::vector<bool> drop(N, false);
    std::vector<double> resid;
    resid.reserve(assigned_idx.size());
    for (int i : assigned_idx) {
      resid.clear();
      const auto vi_cat = cat_vec(assignment[i]);
      for (int j : assigned_idx) {
        if (j == i)
          continue;
        double obs = dot3(v_cam[i], v_cam[j]);
        const auto vj_cat = cat_vec(assignment[j]);
        double cat = dot3(vi_cat, vj_cat);
        resid.push_back(std::abs(obs - cat));
      }
      size_t mid = resid.size() / 2;
      std::nth_element(resid.begin(), resid.begin() + mid, resid.end());
      double med = resid[mid];
      if (med > cutoff) {
        drop[i] = true;
        changed = true;
      }
    }

    if (!changed)
      break;
    for (int i = 0; i < N; ++i) {
      if (drop[i]) {
        assignment[i] = -1;
        --inliers;
      }
    }
  }

  return inliers;
}

} // namespace

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
  std::vector<std::array<double, 3>> v_cam = compute_v_cam(image_stars, camera);

  // Cache catalog unit vectors (one per HIP) to avoid repeated map lookups.
  // Returns by value: holding a reference into an unordered_map across a
  // subsequent insert would be UB on rehash.
  std::unordered_map<int, std::array<double, 3>> cat_vec_cache;
  std::function<std::array<double, 3>(int)> cat_vec =
      [&](int hip) -> std::array<double, 3> {
    auto it = cat_vec_cache.find(hip);
    if (it != cat_vec_cache.end())
      return it->second;
    CatalogStar s = db.get_star(hip);
    std::array<double, 3> v = {s.x, s.y, s.z};
    cat_vec_cache[hip] = v;
    return v;
  };

  // Threshold for declaring a result acceptable. Both passes share it.
  const int min_inliers = std::min(N, std::max(4, N / 4));

  // === 3b.0b: Coarse-refine-reidentify ===
  //
  // First try the standard tight pyramid pass. If it returns enough inliers
  // the camera is already well-calibrated and no refinement is needed —
  // this keeps the well-calibrated path zero-cost.
  //
  // If the tight pass fails to gather enough inliers, the most common
  // cause for real cameras is a small FOV calibration drift (e.g. ~0.4%
  // on the tetra3 IMX265 alt60 fixture) that pushes correct pair cosines
  // just outside the default cos_tolerance of 1e-5. Run a coarse pass with
  // a looser tolerance to recover those pairs, estimate a single focal
  // scale factor from the matched-pair geometry, rescale the camera, and
  // rerun the tight pyramid against the corrected v_cam.
  //
  // To keep the coarse pass tractable on dense scenes (10x looser tolerance
  // expands the candidate set quadratically), we cap it at the first
  // COARSE_N centroids — that's enough to lock onto the constellation and
  // estimate a robust scale factor, but small enough that the per-pair
  // catalog search doesn't blow up.
  std::vector<int> best_assignment;
  int best_inliers;
  std::tie(best_assignment, best_inliers) =
      run_pyramid(v_cam, db, cos_tolerance, cat_vec);

  if (best_inliers < min_inliers) {
    constexpr int COARSE_N = 10;
    const double cos_tolerance_coarse =
        std::max(cos_tolerance * 10.0, 1e-4);
    const int M = std::min(N, COARSE_N);

    // Run coarse pyramid on the first M centroids only.
    std::vector<std::array<double, 3>> v_cam_coarse(v_cam.begin(),
                                                    v_cam.begin() + M);
    std::vector<int> coarse_assignment;
    int coarse_inliers;
    std::tie(coarse_assignment, coarse_inliers) =
        run_pyramid(v_cam_coarse, db, cos_tolerance_coarse, cat_vec);

    if (coarse_inliers >= std::min(M, std::max(4, M / 4))) {
      double s = estimate_scale_factor(coarse_assignment, v_cam_coarse, db,
                                        cat_vec);
      if (std::getenv("STARTRACKER_DEBUG_REFINE")) {
        std::fprintf(stderr, "[identify_stars] refine: tight inliers=%d, "
                            "coarse inliers=%d/%d, scale s=%.6f\n",
                    best_inliers, coarse_inliers, M, s);
      }
      if (std::abs(s - 1.0) > 1e-4) {
        CameraModel refined = camera;
        refined.focal_x *= s;
        refined.focal_y *= s;
        auto v_cam_refined = compute_v_cam(image_stars, refined);

        std::vector<int> refined_assignment;
        int refined_inliers;
        std::tie(refined_assignment, refined_inliers) =
            run_pyramid(v_cam_refined, db, cos_tolerance, cat_vec);

        if (refined_inliers > best_inliers) {
          v_cam = std::move(v_cam_refined);
          best_assignment = std::move(refined_assignment);
          best_inliers = refined_inliers;
        }
      }
    }
  }

  // Require enough inliers to outrank random coincidence.
  if (best_inliers < min_inliers)
    return {};

  // === 3b.3: Cross-verification ===
  //
  // For each currently-assigned star, compute its median |obs - cat|
  // residual against the other assigned stars. Drop stars whose median
  // residual exceeds 3 * cos_tolerance — those are almost certainly
  // mis-labelled. Iterate until stable (max 3 passes).
  best_inliers = cross_verify(best_assignment, v_cam, cos_tolerance, cat_vec);
  if (best_inliers < min_inliers)
    return {};

  std::vector<IdentifiedStar> identified;
  identified.reserve(static_cast<size_t>(best_inliers));
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
