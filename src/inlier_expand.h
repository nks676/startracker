#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

#include "catalog.h"
#include "image_processing.h"
#include "wahba.h"

// Phase 3g.1 (completion): tight-attitude inlier-expansion primitives and the
// 4+1 pattern-verification helper, extracted from identification.cpp.
// Header-only because every entry point is templated on `CatVec` (3g.2:
// a lightweight closure passed in by identify_stars) and we want the
// compiler to inline the cat_vec call across module boundaries.
//
// Public surface: top_n_indices, k_nearest_within_radius, try_verify_candidate,
// expand_inliers, expand_inliers_tight, refine_and_reexpand, plus the
// 5th-star verify cosine constant. All identifiers stay at file scope to
// match the in-place call sites that the original identification.cpp used —
// detail helpers (dot3, apply_rotation_t) hide inside an internal namespace.

namespace inlier_expand_detail {

inline double dot3(const std::array<double, 3> &a,
                   const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// v_inertial = R^T * v_camera (R sends inertial→camera).
inline std::array<double, 3>
apply_rotation_t(const double R[3][3], const std::array<double, 3> &v) {
  return {R[0][0] * v[0] + R[1][0] * v[1] + R[2][0] * v[2],
          R[0][1] * v[0] + R[1][1] * v[1] + R[2][1] * v[2],
          R[0][2] * v[0] + R[1][2] * v[1] + R[2][2] * v[2]};
}

} // namespace inlier_expand_detail

// 5th-star verify cosine tolerance (≈0.05° / cos ≈ 1 - 3.8e-7). Shared by
// try_verify_candidate, refine_and_reexpand's inlier-expansion gate, and the
// post-pattern FOV-scale refine loop inside identify_stars.
inline constexpr double kFifthStarVerifyCosTol = 4e-7;

// Build the brightest-N subset of centroid indices.
//
// Phase 3g.9: comparator is (peak desc, intensity desc, y asc, x asc) — same
// tie-break as main.cpp's CENTROID_CAP pre-trim. Without the tie-break, the
// pattern-path seed pool could shift between runs (or between 8-bit-stretched
// and 16-bit-native inputs of the same scene), routing through different
// seeds and producing the alt60 35 ms identify-stage variance.
inline std::vector<int>
top_n_indices(const std::vector<StarCentroid> &image_stars, int N) {
  std::vector<int> idx(image_stars.size());
  std::iota(idx.begin(), idx.end(), 0);
  auto better = [&](int a, int b) {
    const StarCentroid &sa = image_stars[a];
    const StarCentroid &sb = image_stars[b];
    if (sa.peak != sb.peak) return sa.peak > sb.peak;
    if (sa.intensity != sb.intensity) return sa.intensity > sb.intensity;
    if (sa.y != sb.y) return sa.y < sb.y;
    return sa.x < sb.x;
  };
  if (static_cast<int>(idx.size()) <= N) {
    std::sort(idx.begin(), idx.end(), better);
    return idx;
  }
  std::partial_sort(idx.begin(), idx.begin() + N, idx.end(), better);
  idx.resize(N);
  return idx;
}

// Given a seed centroid (global index `seed_idx`) and a list of candidate
// neighbor centroid indices, return up to `k` nearest neighbors whose
// camera-frame angular distance from the seed is ≤ `max_radius_cos`
// (i.e., dot >= cos_radius). Returns indices in ascending angular-distance
// order (largest dot first).
inline std::vector<int>
k_nearest_within_radius(int seed_idx, const std::vector<int> &neighbor_pool,
                        const std::vector<std::array<double, 3>> &v_cam,
                        double cos_radius, int k) {
  std::vector<std::pair<double, int>> dots;
  dots.reserve(neighbor_pool.size());
  for (int idx : neighbor_pool) {
    if (idx == seed_idx) continue;
    double d = v_cam[seed_idx][0] * v_cam[idx][0] +
               v_cam[seed_idx][1] * v_cam[idx][1] +
               v_cam[seed_idx][2] * v_cam[idx][2];
    if (d < cos_radius) continue;
    dots.emplace_back(d, idx);
  }
  int take = std::min(k, static_cast<int>(dots.size()));
  if (take == 0) return {};
  std::partial_sort(dots.begin(), dots.begin() + take, dots.end(),
                    [](const auto &a, const auto &b) { return a.first > b.first; });
  std::vector<int> out(take);
  for (int i = 0; i < take; ++i) out[i] = dots[i].second;
  return out;
}

// Refine an assignment by adding more inliers beyond the 4+1 verified core.
// For each currently-assigned star pair (A, image i), and for each unassigned
// centroid k, project the catalog's neighbor of A at the observed angle and
// see if any candidate HIP satisfies the geometry against B (the brightest
// matched star). Mirrors the pyramid's expansion step but starts from a
// known-good seed pair so it's much cheaper.
template <typename CatVec>
int expand_inliers(std::vector<int> &assignment,
                   const std::vector<std::array<double, 3>> &v_cam,
                   const StarDatabase &db, double cos_tolerance,
                   const CatVec &cat_vec) {
  using inlier_expand_detail::dot3;
  const int N = static_cast<int>(v_cam.size());
  int seed_i = -1, seed_j = -1;
  for (int i = 0; i < N; ++i) {
    if (assignment[i] < 0) continue;
    if (seed_i < 0) seed_i = i;
    else if (seed_j < 0) {
      seed_j = i;
      break;
    }
  }
  if (seed_j < 0) return 0;
  int A = assignment[seed_i];
  int B = assignment[seed_j];
  auto vb = cat_vec(B);
  std::unordered_set<int> used;
  int inliers = 0;
  for (int i = 0; i < N; ++i) {
    if (assignment[i] >= 0) {
      used.insert(assignment[i]);
      ++inliers;
    }
  }
  for (int k = 0; k < N; ++k) {
    if (assignment[k] >= 0) continue;
    double obs_ik = dot3(v_cam[seed_i], v_cam[k]);
    double obs_jk = dot3(v_cam[seed_j], v_cam[k]);
    auto candidates = db.find_partners(A, obs_ik, cos_tolerance);
    int best_C = -1;
    double best_err = cos_tolerance;
    for (int C : candidates) {
      if (used.count(C)) continue;
      auto vc = cat_vec(C);
      double cat_bc = dot3(vb, vc);
      double err = std::abs(cat_bc - obs_jk);
      if (err < best_err) {
        best_err = err;
        best_C = C;
      }
    }
    if (best_C != -1) {
      assignment[k] = best_C;
      used.insert(best_C);
      ++inliers;
    }
  }
  return inliers;
}

// Phase 3e.5 (Change 2): tight inlier expansion driven by a TRIAD attitude.
//
// For each currently-unassigned centroid, project it from camera frame into
// inertial frame using the TRIAD rotation `R` (camera<-inertial), then search
// the catalog for the nearest star to that predicted direction. Acceptance
// gate is the same 0.05° (cos ≥ 1 - 4e-7) threshold used by the 5th-star
// verify, so newly added stars are TIGHT — much tighter than the cos_tolerance
// gate used by the pyramid-style expand_inliers, which absorbs FOV
// miscalibration at the cost of letting wrong matches slip through.
//
// We can't iterate the full catalog directly (catalog API exposes only
// get_star / find_partners), so the "brute scan" is realised by enumerating
// candidate HIPs via find_partners using EACH currently-matched HIP as an
// anchor at the predicted cos(anchor, predicted_inertial). The union is then
// scored by direct dot-product against the predicted inertial direction.
//
// Returns the new total inlier count. `assignment` is mutated in place.
template <typename CatVec>
int expand_inliers_tight(std::vector<int> &assignment,
                         const std::vector<std::array<double, 3>> &v_cam,
                         const StarDatabase &db, const CatVec &cat_vec,
                         const double R[3][3], double accept_cos_tol) {
  using inlier_expand_detail::apply_rotation_t;
  const int N = static_cast<int>(v_cam.size());

  std::vector<std::pair<int, int>> matched;
  matched.reserve(8);
  std::unordered_set<int> used;
  used.reserve(16);
  for (int i = 0; i < N; ++i) {
    if (assignment[i] >= 0) {
      matched.emplace_back(i, assignment[i]);
      used.insert(assignment[i]);
    }
  }
  if (matched.empty())
    return 0;

  std::vector<std::array<double, 3>> anchor_vecs;
  anchor_vecs.reserve(matched.size());
  for (const auto &m : matched) {
    try {
      anchor_vecs.push_back(cat_vec(m.second));
    } catch (...) {
      anchor_vecs.push_back({0.0, 0.0, 0.0});
    }
  }

  int inliers = static_cast<int>(matched.size());

  // Wider find_partners tolerance so that ±0.05° centroid noise plus mild FOV
  // miscal doesn't drop the true catalog star out of the candidate ring. The
  // final acceptance is by direct dot-product against the predicted inertial
  // direction, not by this lookup tolerance.
  constexpr double kPartnerLookupTol = 1e-3;

  for (int k = 0; k < N; ++k) {
    if (assignment[k] >= 0)
      continue;

    auto v_iner = apply_rotation_t(R, v_cam[k]);

    int best_C = -1;
    double best_dot = -2.0;
    for (size_t ai = 0; ai < matched.size(); ++ai) {
      const auto &va = anchor_vecs[ai];
      double cos_pred =
          va[0] * v_iner[0] + va[1] * v_iner[1] + va[2] * v_iner[2];
      if (cos_pred < -1.0) cos_pred = -1.0;
      if (cos_pred > 1.0) cos_pred = 1.0;
      auto candidates =
          db.find_partners(matched[ai].second, cos_pred, kPartnerLookupTol);
      for (int C : candidates) {
        if (used.count(C))
          continue;
        std::array<double, 3> vc;
        try {
          vc = cat_vec(C);
        } catch (...) {
          continue;
        }
        double d = vc[0] * v_iner[0] + vc[1] * v_iner[1] + vc[2] * v_iner[2];
        if (d > best_dot) {
          best_dot = d;
          best_C = C;
        }
      }
    }

    if (best_C != -1 && best_dot >= 1.0 - accept_cos_tol) {
      assignment[k] = best_C;
      used.insert(best_C);
      ++inliers;
    }
  }

  return inliers;
}

// Phase 3e.5: replace the TRIAD-on-the-original-4 attitude with a QUEST
// attitude over the full current inlier set, then re-run tight inlier
// expansion using that refined attitude. One pass of refine→re-expand is
// usually enough; the second pass typically adds 0 new inliers.
//
// Mutates `assignment` and writes the refined rotation into `R_out` (which
// also serves as the seed attitude on entry). Returns the final inlier count.
//
// Uses wahba::quest_R with the tight identification-side settings (30 Newton
// iters / 1e-14 tol) — the public attitude path in estimation.cpp uses 20 /
// 1e-12.
template <typename CatVec>
int refine_and_reexpand(std::vector<int> &assignment,
                        const std::vector<std::array<double, 3>> &v_cam,
                        const StarDatabase &db, const CatVec &cat_vec,
                        double R_out[3][3], double accept_cos_tol) {
  const int N = static_cast<int>(v_cam.size());

  int inliers = expand_inliers_tight(assignment, v_cam, db, cat_vec, R_out,
                                      accept_cos_tol);
  if (inliers < 3) return inliers;

  std::vector<std::array<double, 3>> v_cam_set, v_iner_set;
  v_cam_set.reserve(static_cast<size_t>(inliers));
  v_iner_set.reserve(static_cast<size_t>(inliers));
  for (int i = 0; i < N; ++i) {
    if (assignment[i] < 0) continue;
    v_cam_set.push_back(v_cam[i]);
    try {
      v_iner_set.push_back(cat_vec(assignment[i]));
    } catch (...) {
      v_cam_set.pop_back();
    }
  }
  double R_refined[3][3];
  if (!wahba::quest_R(v_cam_set, v_iner_set, R_refined,
                      /*max_iters=*/30, /*tol_lambda=*/1e-14)) {
    return inliers;
  }
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = R_refined[r][c];

  int new_inliers = expand_inliers_tight(assignment, v_cam, db, cat_vec,
                                          R_out, accept_cos_tol);
  return new_inliers;
}

// Verify a candidate 4-tuple pattern match by:
//   1. Pair observed[centroid_idx_canonical[i]] (i in 0..3) with hips[i].
//   2. Solve TRIAD on the pair with widest catalog separation (best
//      numerical conditioning).
//   3. Check all 6 candidate pair cosines are self-consistent against the
//      catalog. Reject without going to the 5th-star step if any pair is
//      grossly off.
//   4. Project each unmatched candidate centroid into inertial and verify
//      at least one lands within `verify_cos_tol` of a catalog star (looked
//      up via db.find_partners off matched HIP 0). On first match, accept.
//
// Returns the (partial) assignment image_idx -> HIP, with -1 for unassigned,
// or an empty vector on rejection. On success, `R_out` receives the TRIAD-
// derived inertial->camera rotation so the caller can drive inlier expansion
// without re-solving.
template <typename CatVec>
std::vector<int>
try_verify_candidate(const std::array<int, 4> &centroid_idx_canonical,
                     const std::array<int, 4> &hips,
                     const std::vector<int> &verify_pool,
                     const std::vector<std::array<double, 3>> &v_cam,
                     const StarDatabase &db, const CatVec &cat_vec,
                     double pair_cos_tol, double verify_cos_tol,
                     double R_out[3][3]) {
  using inlier_expand_detail::apply_rotation_t;
  using inlier_expand_detail::dot3;
  const int N = static_cast<int>(v_cam.size());
  std::vector<int> assignment(N, -1);

  std::array<int, 4> centroid_idx = centroid_idx_canonical;
  std::array<std::array<double, 3>, 4> v_obs{};
  std::array<std::array<double, 3>, 4> v_cat{};
  for (int i = 0; i < 4; ++i) {
    v_obs[i] = v_cam[centroid_idx[i]];
    try {
      v_cat[i] = cat_vec(hips[i]);
    } catch (...) {
      return {};
    }
  }

  // Pair-consistency gate (6 inter-pair cosines vs catalog within pair_cos_tol).
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      double obs = dot3(v_obs[i], v_obs[j]);
      double cat = dot3(v_cat[i], v_cat[j]);
      if (std::abs(obs - cat) > pair_cos_tol) {
        return {};
      }
    }
  }

  // Pick the two pairs (i, j) with the widest catalog separation (smallest
  // cosine) for TRIAD numerical stability.
  int best_a = 0, best_b = 1;
  double smallest_cos = 2.0;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      double c = dot3(v_cat[i], v_cat[j]);
      if (c < smallest_cos) {
        smallest_cos = c;
        best_a = i;
        best_b = j;
      }
    }
  }
  double R[3][3];
  wahba::triad(v_obs[best_a], v_obs[best_b], v_cat[best_a], v_cat[best_b], R);
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = R[r][c];

  for (int i = 0; i < 4; ++i) {
    assignment[centroid_idx[i]] = hips[i];
  }

  // 5th-star verify: brightest unmatched centroid in verify_pool, projected
  // into inertial via TRIAD, looked up against matched HIP 0's partner ring.
  bool fifth_verified = false;
  int fifth_centroid = -1;
  int fifth_hip = -1;
  std::unordered_set<int> matched_centroids(centroid_idx.begin(),
                                            centroid_idx.end());
  for (int pool_local = 0;
       pool_local < static_cast<int>(verify_pool.size()) && !fifth_verified;
       ++pool_local) {
    int ci = verify_pool[pool_local];
    if (matched_centroids.count(ci)) continue;

    auto v_cam_5 = v_cam[ci];
    auto v_iner_5 = apply_rotation_t(R, v_cam_5);

    int anchor_hip = hips[0];
    auto v_anchor = v_cat[0];
    double cos_pred = dot3(v_anchor, v_iner_5);
    if (cos_pred < -1.0) cos_pred = -1.0;
    if (cos_pred > 1.0) cos_pred = 1.0;
    auto candidates = db.find_partners(anchor_hip, cos_pred, 1e-3);
    int best_C = -1;
    double best_dot = -2.0;
    for (int C : candidates) {
      if (C == hips[0] || C == hips[1] || C == hips[2] || C == hips[3])
        continue;
      std::array<double, 3> vc;
      try {
        vc = cat_vec(C);
      } catch (...) {
        continue;
      }
      double d = dot3(vc, v_iner_5);
      if (d > best_dot) {
        best_dot = d;
        best_C = C;
      }
    }
    if (best_C != -1 && best_dot >= 1.0 - verify_cos_tol) {
      fifth_verified = true;
      fifth_centroid = ci;
      fifth_hip = best_C;
    }
  }

  if (!fifth_verified) {
    return {};
  }

  assignment[fifth_centroid] = fifth_hip;
  return assignment;
}
