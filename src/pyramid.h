#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "camera_model.h"
#include "catalog.h"
#include "identification.h" // IdentifiedStar struct
#include "image_processing.h"

// Phase 3g.1 (completion): pyramid (geometric voting) identification path,
// extracted from identification.cpp. Header-only because every entry point
// is templated on `CatVec` (the lambda passed by identify_stars). Pyramid is
// now the fallback when the pattern-hash path misses; identify_stars_pyramid
// runs the original 3a-era logic plus the 3b.0b coarse-refine-reidentify
// FOV-scale absorption pass plus the 3b.3 cross-verification.

namespace pyramid_detail {

inline double dot3(const std::array<double, 3> &a,
                   const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute camera-frame unit vectors from pixel centroids by inverting the
// camera distortion model.
inline std::vector<std::array<double, 3>>
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

} // namespace pyramid_detail

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
//   4. Keep the seed that produces the most inliers.
//
// Returns the best assignment (image index -> HIP, or -1) and inlier count.
//
// Phase 3g.5: outer loops exit early once best_inliers ≥ max(12, 0.7·N) — the
// marginal value of continuing is near zero, and the dense-frame worst case
// gets dramatically faster.
template <typename CatVec>
std::pair<std::vector<int>, int>
run_pyramid(const std::vector<std::array<double, 3>> &v_cam,
            const StarDatabase &db, double cos_tolerance,
            const CatVec &cat_vec) {
  using pyramid_detail::dot3;
  const int N = static_cast<int>(v_cam.size());
  std::vector<int> best_assignment(N, -1);
  int best_inliers = 0;

  const int kEarlyExitInliers = std::max(12, (7 * N) / 10);

  for (int i = 0; i < N; ++i) {
    if (best_inliers >= kEarlyExitInliers)
      break;
    for (int j = i + 1; j < N; ++j) {
      double obs_ij = dot3(v_cam[i], v_cam[j]);
      auto seeds = db.find_pairs_kvec(obs_ij, cos_tolerance);
      for (const auto &seed : seeds) {
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
      if (best_inliers >= kEarlyExitInliers)
        break;
    }
  }

  return {std::move(best_assignment), best_inliers};
}

// Estimate the FOV scale-factor s = (true focal) / (assumed focal) from a
// coarse identification result. For every pair of assigned image stars we
// compare the observed inter-star angle (from v_cam, computed with the
// assumed focal) against the catalog angle, and take the median of the
// ratios. If the assumed focal is too small, observed unit vectors splay
// outward, so observed_angle > catalog_angle and s > 1.
template <typename CatVec>
double estimate_scale_factor(
    const std::vector<int> &assignment,
    const std::vector<std::array<double, 3>> &v_cam,
    const StarDatabase &db, const CatVec &cat_vec) {
  using pyramid_detail::dot3;
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
template <typename CatVec>
int cross_verify(std::vector<int> &assignment,
                 const std::vector<std::array<double, 3>> &v_cam,
                 double cos_tolerance, const CatVec &cat_vec) {
  using pyramid_detail::dot3;
  const int N = static_cast<int>(assignment.size());
  const double cutoff = 3.0 * cos_tolerance;

  int inliers = 0;
  for (int i = 0; i < N; ++i)
    if (assignment[i] >= 0)
      ++inliers;

  for (int pass = 0; pass < 3; ++pass) {
    bool changed = false;

    std::vector<int> assigned_idx;
    assigned_idx.reserve(static_cast<size_t>(inliers));
    for (int i = 0; i < N; ++i)
      if (assignment[i] >= 0)
        assigned_idx.push_back(i);

    if (assigned_idx.size() < 3)
      break;

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

// Pyramid-path identify_stars: pre-3e fallback used when the pattern-hash
// path misses (or when no pattern catalog is loaded). Runs run_pyramid,
// optional coarse-refine-reidentify under a relaxed cos_tolerance to absorb
// FOV-scale miscal (3b.0b), and cross_verify (3b.3).
template <typename CatVec>
std::vector<IdentifiedStar>
identify_stars_pyramid(const std::vector<StarCentroid> &image_stars,
                       const PinholeCamera &camera, const StarDatabase &db,
                       double cos_tolerance,
                       std::vector<std::array<double, 3>> v_cam,
                       const CatVec &cat_vec) {
  int N = image_stars.size();
  const int min_inliers = std::min(N, std::max(4, N / 4));

  // === 3b.0b: Coarse-refine-reidentify ===
  std::vector<int> best_assignment;
  int best_inliers;
  std::tie(best_assignment, best_inliers) =
      run_pyramid(v_cam, db, cos_tolerance, cat_vec);

  if (best_inliers < min_inliers) {
    constexpr int COARSE_N = 10;
    const double cos_tolerance_coarse =
        std::max(cos_tolerance * 10.0, 1e-4);
    const int M = std::min(N, COARSE_N);

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
        auto v_cam_refined = pyramid_detail::compute_v_cam(image_stars, refined);

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

  if (best_inliers < min_inliers)
    return {};

  // === 3b.3: Cross-verification ===
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
