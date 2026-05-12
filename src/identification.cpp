#include "identification.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <functional>
#include <numeric>
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
      auto seeds = db.find_pairs_kvec(obs_ij, cos_tolerance);
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

// === Phase 3e.3: pattern-hash identification helpers ===

// Edge enumeration in (i, j) local-index order. Must match
// tools/generate_catalog.py:EDGE_PAIRS exactly so canonical-order bit-flips
// between generator and runtime are impossible.
constexpr std::array<std::pair<int, int>, 6> kEdgePairs = {{
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
}};
constexpr int kQuantBits = 10;
constexpr int kQuantScale = (1 << kQuantBits) - 1; // 1023

// Recover the horizontal half-FOV (radians) from the camera model. Used as
// the upper bound on the 4-star pattern radius — the pattern catalog at FOV
// bin = 2·half-FOV stores only patterns with all 4 stars within that radius
// of the seed.
inline double camera_half_fov_rad(const CameraModel &camera) {
  return std::atan(camera.frame_width / (2.0 * camera.focal_x));
}


// Solve TRIAD on two camera-frame / inertial-frame vector pairs. Returns the
// 3x3 inertial->camera rotation matrix R such that v_cam ≈ R * v_inertial.
// (Mirrors the convention used elsewhere: R = M_W * M_V^T.)
void triad_rotation(const std::array<double, 3> &W1,
                    const std::array<double, 3> &W2,
                    const std::array<double, 3> &V1,
                    const std::array<double, 3> &V2,
                    double R_out[3][3]) {
  auto cross = [](const std::array<double, 3> &a,
                  const std::array<double, 3> &b) {
    return std::array<double, 3>{a[1] * b[2] - a[2] * b[1],
                                 a[2] * b[0] - a[0] * b[2],
                                 a[0] * b[1] - a[1] * b[0]};
  };
  auto norm = [](const std::array<double, 3> &v) {
    double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return std::array<double, 3>{v[0] / n, v[1] / n, v[2] / n};
  };

  auto t1W = W1;
  auto t2W = norm(cross(W1, W2));
  auto t3W = cross(t1W, t2W);
  auto t1V = V1;
  auto t2V = norm(cross(V1, V2));
  auto t3V = cross(t1V, t2V);

  double M_W[3][3] = {{t1W[0], t2W[0], t3W[0]},
                      {t1W[1], t2W[1], t3W[1]},
                      {t1W[2], t2W[2], t3W[2]}};
  double M_V_T[3][3] = {{t1V[0], t1V[1], t1V[2]},
                        {t2V[0], t2V[1], t2V[2]},
                        {t3V[0], t3V[1], t3V[2]}};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = M_W[r][0] * M_V_T[0][c] + M_W[r][1] * M_V_T[1][c] +
                    M_W[r][2] * M_V_T[2][c];
}

inline std::array<double, 3>
apply_rotation_t(const double R[3][3], const std::array<double, 3> &v) {
  // v_inertial = R^T * v_camera (R sends inertial→camera; the test uses R*v
  // direction in the pyramid code path, so we expose its transpose here).
  return {R[0][0] * v[0] + R[1][0] * v[1] + R[2][0] * v[2],
          R[0][1] * v[0] + R[1][1] * v[1] + R[2][1] * v[2],
          R[0][2] * v[0] + R[1][2] * v[1] + R[2][2] * v[2]};
}

} // namespace

// Phase 3e.5 (Change 1): generate the full set of plausible canonical keys
// for a 4-star pattern under centroid noise.
//
// The canonical-order key generation sorts the 6 inter-star angular
// distances. When two adjacent sorted distances are within a noise tolerance,
// the sort can flip them under realistic centroid noise — and that flip
// shifts the per-star edge signature, which shifts the canonical permutation,
// which shifts the 5-quantized-ratio packing into a non-adjacent hash bucket.
// `find_pattern_tolerant`'s ±1-quantum probing can't recover from such
// shifts (the slot ordering changes by up to ±N quanta where N is the
// difference between adjacent ratios).
//
// Fix: enumerate the small set of plausible edge sort orderings produced by
// flipping pairs of adjacent ranks whose underlying distances are within
// `noise_tol` (in radians). For each ordering, run the rest of the canonical
// procedure (signatures, canonical permutation, quantize, pack) and emit the
// resulting key plus the centroid→hip slot mapping. Duplicates are filtered
// at the call site (cheap hash-set dedupe).
//
// For typical 4-star patterns 0–2 adjacent pairs are close enough to flip,
// so the output is 1–4 keys. Each is a 243-probe `find_pattern_tolerant`
// lookup — still well within the µs budget.
//
// On a noise-free input, this returns exactly one key, equal to
// `pattern_key_canonical(vecs4, ids4, _)` — so it is a drop-in superset.
// Declaration is in identification.h; defined below.

uint64_t
pattern_key_canonical(const std::array<std::array<double, 3>, 4> &vecs4,
                      const std::array<int, 4> &ids4,
                      std::array<int, 4> &out_canonical) {
  // 1. Six pairwise angular distances (radians, acos of clamped dot).
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    double d = vecs4[i][0] * vecs4[j][0] + vecs4[i][1] * vecs4[j][1] +
               vecs4[i][2] * vecs4[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }

  // 2. Sort the 6 edges by (distance, EDGE_PAIRS) ascending. Tie-break by
  // the EDGE_PAIRS (i, j) tuple — equivalent to Python's
  // sorted(..., key=lambda e: (dists[e], EDGE_PAIRS[e])).
  std::array<int, 6> edge_order = {0, 1, 2, 3, 4, 5};
  std::sort(edge_order.begin(), edge_order.end(), [&](int a, int b) {
    if (dists[a] != dists[b]) return dists[a] < dists[b];
    // Lex compare on EDGE_PAIRS tuple.
    if (kEdgePairs[a].first != kEdgePairs[b].first)
      return kEdgePairs[a].first < kEdgePairs[b].first;
    return kEdgePairs[a].second < kEdgePairs[b].second;
  });

  // 3. Per-star edge-rank signature: sorted 3-tuple of ranks (in 0..5) of the
  // 3 edges incident to each star.
  std::array<std::array<int, 3>, 4> sig{};
  std::array<int, 4> sig_fill = {0, 0, 0, 0};
  for (int rank = 0; rank < 6; ++rank) {
    int e = edge_order[rank];
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    sig[i][sig_fill[i]++] = rank;
    sig[j][sig_fill[j]++] = rank;
  }
  for (int s = 0; s < 4; ++s)
    std::sort(sig[s].begin(), sig[s].end());

  // 4. Sort the 4 local-star indices ascending by (signature, input_id).
  std::array<int, 4> canonical = {0, 1, 2, 3};
  std::sort(canonical.begin(), canonical.end(), [&](int a, int b) {
    for (int k = 0; k < 3; ++k) {
      if (sig[a][k] != sig[b][k]) return sig[a][k] < sig[b][k];
    }
    return ids4[a] < ids4[b];
  });
  out_canonical = canonical;

  // 5. Quantize the 5 smallest sorted distances normalized by the largest.
  // sorted_dists[5] = dists[edge_order[5]] is the largest. In a non-degenerate
  // pattern it is strictly positive; we still guard for safety.
  double largest = dists[edge_order[5]];
  if (!(largest > 0.0)) {
    // Degenerate (coincident stars): return 0 key. Pattern lookups won't
    // match anything useful, but at least no NaN / div-by-zero.
    return 0ULL;
  }
  std::array<uint64_t, 5> q{};
  for (int k = 0; k < 5; ++k) {
    double r = dists[edge_order[k]] / largest;
    double v = std::round(r * static_cast<double>(kQuantScale));
    if (v < 0.0) v = 0.0;
    if (v > static_cast<double>(kQuantScale)) v = static_cast<double>(kQuantScale);
    q[k] = static_cast<uint64_t>(v);
  }

  return (q[4] << 40) | (q[3] << 30) | (q[2] << 20) | (q[1] << 10) | q[0];
}

// Compute the canonical (key, centroid_canonical_order) for a 4-star pattern
// GIVEN a specific edge ordering. Returns {0, ...} on degenerate input.
// Mirrors steps 3–6 of pattern_key_canonical; the only difference is that
// the caller supplies edge_order rather than the sort-determined ordering.
static std::pair<uint64_t, std::array<int, 4>>
pattern_key_with_edge_order(const std::array<double, 6> &dists,
                            const std::array<int, 6> &edge_order,
                            const std::array<int, 4> &ids4) {
  std::array<std::array<int, 3>, 4> sig{};
  std::array<int, 4> sig_fill = {0, 0, 0, 0};
  for (int rank = 0; rank < 6; ++rank) {
    int e = edge_order[rank];
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    sig[i][sig_fill[i]++] = rank;
    sig[j][sig_fill[j]++] = rank;
  }
  for (int s = 0; s < 4; ++s)
    std::sort(sig[s].begin(), sig[s].end());

  std::array<int, 4> canonical = {0, 1, 2, 3};
  std::sort(canonical.begin(), canonical.end(), [&](int a, int b) {
    for (int k = 0; k < 3; ++k) {
      if (sig[a][k] != sig[b][k]) return sig[a][k] < sig[b][k];
    }
    return ids4[a] < ids4[b];
  });

  double largest = dists[edge_order[5]];
  if (!(largest > 0.0)) return {0ULL, canonical};

  std::array<uint64_t, 5> q{};
  for (int k = 0; k < 5; ++k) {
    double r = dists[edge_order[k]] / largest;
    double v = std::round(r * static_cast<double>(kQuantScale));
    if (v < 0.0) v = 0.0;
    if (v > static_cast<double>(kQuantScale)) v = static_cast<double>(kQuantScale);
    q[k] = static_cast<uint64_t>(v);
  }
  uint64_t key =
      (q[4] << 40) | (q[3] << 30) | (q[2] << 20) | (q[1] << 10) | q[0];
  return {key, canonical};
}

std::vector<std::pair<uint64_t, std::array<int, 4>>>
pattern_keys_noise_robust(const std::array<std::array<double, 3>, 4> &vecs4,
                          const std::array<int, 4> &ids4, double noise_tol) {
  // Compute the 6 pairwise distances (matches step 1 of pattern_key_canonical).
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    double d = vecs4[i][0] * vecs4[j][0] + vecs4[i][1] * vecs4[j][1] +
               vecs4[i][2] * vecs4[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }

  // Build the base sorted edge order (same tie-break as canonical).
  std::array<int, 6> base_order = {0, 1, 2, 3, 4, 5};
  std::sort(base_order.begin(), base_order.end(), [&](int a, int b) {
    if (dists[a] != dists[b]) return dists[a] < dists[b];
    if (kEdgePairs[a].first != kEdgePairs[b].first)
      return kEdgePairs[a].first < kEdgePairs[b].first;
    return kEdgePairs[a].second < kEdgePairs[b].second;
  });

  // Identify uncertain adjacent pairs (ranks i, i+1) whose distance gap is
  // within noise_tol. At most ~3 in practice for a noisy 4-star scene.
  std::vector<int> uncertain_pos;
  uncertain_pos.reserve(5);
  for (int i = 0; i < 5; ++i) {
    double gap = dists[base_order[i + 1]] - dists[base_order[i]];
    if (gap < noise_tol) uncertain_pos.push_back(i);
  }

  // Enumerate 2^k orderings (k = uncertain count). For each subset of
  // uncertain positions, swap the corresponding adjacent ranks in the order.
  std::vector<std::pair<uint64_t, std::array<int, 4>>> out;
  std::unordered_set<uint64_t> seen;
  const size_t num_subsets = static_cast<size_t>(1) << uncertain_pos.size();
  out.reserve(num_subsets);
  for (size_t mask = 0; mask < num_subsets; ++mask) {
    // Apply the swap pattern atop base_order. Swaps may overlap (adjacent
    // positions); we apply them left-to-right so the result is a valid
    // ordering. Overlap is handled by the dedupe below.
    std::array<int, 6> order = base_order;
    for (size_t bit = 0; bit < uncertain_pos.size(); ++bit) {
      if (mask & (1u << bit)) {
        int pos = uncertain_pos[bit];
        std::swap(order[pos], order[pos + 1]);
      }
    }
    auto kp = pattern_key_with_edge_order(dists, order, ids4);
    if (!seen.insert(kp.first).second) continue;
    out.emplace_back(kp);
  }
  return out;
}

namespace {

// Build the brightest-N subset of centroid indices, preserving the input
// order among the selected indices. Inputs are assumed to already be ranked
// by peak (the main binary does this before calling identify_stars). For
// safety we re-sort by peak descending and keep the top N.
std::vector<int> top_n_indices(const std::vector<StarCentroid> &image_stars,
                               int N) {
  std::vector<int> idx(image_stars.size());
  std::iota(idx.begin(), idx.end(), 0);
  if (static_cast<int>(idx.size()) <= N) {
    // Sort all by peak desc to stabilize "brightest first" iteration order
    // even if caller didn't pre-rank.
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
      return image_stars[a].peak > image_stars[b].peak;
    });
    return idx;
  }
  std::partial_sort(idx.begin(), idx.begin() + N, idx.end(), [&](int a, int b) {
    return image_stars[a].peak > image_stars[b].peak;
  });
  idx.resize(N);
  return idx;
}

// Given a seed centroid (global index `seed_idx`) and a list of candidate
// neighbor centroid indices, return up to `k` nearest neighbors whose
// camera-frame angular distance from the seed is ≤ `max_radius_cos`
// (i.e., dot >= cos_radius). Returns indices in ascending angular-distance
// order (largest dot first). May return fewer than k if there aren't enough
// in radius.
std::vector<int> k_nearest_within_radius(
    int seed_idx, const std::vector<int> &neighbor_pool,
    const std::vector<std::array<double, 3>> &v_cam, double cos_radius,
    int k) {
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
int expand_inliers(std::vector<int> &assignment,
                   const std::vector<std::array<double, 3>> &v_cam,
                   const StarDatabase &db, double cos_tolerance,
                   const std::function<std::array<double, 3>(int)> &cat_vec) {
  const int N = static_cast<int>(v_cam.size());
  // Pick the first two assigned stars as the seed pair.
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
// We can't iterate the full catalog directly from identification.cpp (catalog
// API exposes only get_star / find_partners), so the "brute scan" is realised
// by enumerating candidate HIPs via find_partners using EACH currently-matched
// HIP as an anchor at the predicted cos(anchor, predicted_inertial). The union
// of those candidate sets is then scored by direct dot-product against the
// predicted inertial direction. With 4+ anchors and pattern-path-tight
// cosines, the true catalog star is reliably in the union; the tight 0.05°
// gate filters out the rest.
//
// Returns the new total inlier count. `assignment` is mutated in place.
int expand_inliers_tight(std::vector<int> &assignment,
                         const std::vector<std::array<double, 3>> &v_cam,
                         const StarDatabase &db,
                         const std::function<std::array<double, 3>(int)> &cat_vec,
                         const double R[3][3], double accept_cos_tol) {
  const int N = static_cast<int>(v_cam.size());

  // Snapshot currently-matched (image_idx, hip) pairs.
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

  // Cache anchor inertial vectors.
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

    // Predict inertial direction: v_iner = R^T * v_cam.
    auto v_iner = apply_rotation_t(R, v_cam[k]);

    // Score candidates pulled from every matched HIP's partner ring at the
    // predicted angle. A real catalog star will show up across multiple
    // anchors; spurious candidates show up against only one.
    int best_C = -1;
    double best_dot = -2.0;
    for (size_t ai = 0; ai < matched.size(); ++ai) {
      const auto &va = anchor_vecs[ai];
      double cos_pred = va[0] * v_iner[0] + va[1] * v_iner[1] + va[2] * v_iner[2];
      if (cos_pred < -1.0) cos_pred = -1.0;
      if (cos_pred > 1.0) cos_pred = 1.0;
      auto candidates = db.find_partners(matched[ai].second, cos_pred,
                                          kPartnerLookupTol);
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

// Phase 3e.5 (Change 2): solve a small Wahba problem and produce the
// inertial->camera rotation matrix. Mirrors estimation.cpp's QUEST but kept
// local so identification.cpp doesn't need to pull in estimation.h. Used only
// to refine the TRIAD attitude after inlier expansion — the public API
// (identify_stars → IdentifiedStar list → estimation.cpp) still does the
// official QUEST downstream.
//
// Returns true on convergence; R_out is unmodified on failure.
bool quest_attitude_local(
    const std::vector<std::array<double, 3>> &v_cam_set,
    const std::vector<std::array<double, 3>> &v_iner_set, double R_out[3][3]) {
  const int M = static_cast<int>(v_cam_set.size());
  if (M < 2) return false;
  const double w = 1.0 / static_cast<double>(M);

  double B[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int k = 0; k < M; ++k) {
    const auto &b = v_cam_set[k];
    const auto &r = v_iner_set[k];
    for (int a = 0; a < 3; ++a)
      for (int bb = 0; bb < 3; ++bb)
        B[a][bb] += w * b[a] * r[bb];
  }

  double S[3][3];
  for (int a = 0; a < 3; ++a)
    for (int bb = 0; bb < 3; ++bb)
      S[a][bb] = B[a][bb] + B[bb][a];

  double sigma = B[0][0] + B[1][1] + B[2][2];
  double Z[3] = {B[2][1] - B[1][2], B[0][2] - B[2][0], B[1][0] - B[0][1]};

  auto det3 = [](const double M3[3][3]) {
    return M3[0][0] * (M3[1][1] * M3[2][2] - M3[1][2] * M3[2][1]) -
           M3[0][1] * (M3[1][0] * M3[2][2] - M3[1][2] * M3[2][0]) +
           M3[0][2] * (M3[1][0] * M3[2][1] - M3[1][1] * M3[2][0]);
  };
  double trace_adj_S = (S[1][1] * S[2][2] - S[1][2] * S[2][1]) +
                       (S[0][0] * S[2][2] - S[0][2] * S[2][0]) +
                       (S[0][0] * S[1][1] - S[0][1] * S[1][0]);
  double detS = det3(S);
  double ZtZ = Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2];
  double SZ[3] = {S[0][0] * Z[0] + S[0][1] * Z[1] + S[0][2] * Z[2],
                  S[1][0] * Z[0] + S[1][1] * Z[1] + S[1][2] * Z[2],
                  S[2][0] * Z[0] + S[2][1] * Z[1] + S[2][2] * Z[2]};
  double ZtSZ = Z[0] * SZ[0] + Z[1] * SZ[1] + Z[2] * SZ[2];
  double S2[3][3];
  for (int a = 0; a < 3; ++a)
    for (int bb = 0; bb < 3; ++bb)
      S2[a][bb] = S[a][0] * S[0][bb] + S[a][1] * S[1][bb] + S[a][2] * S[2][bb];
  double S2Z[3] = {S2[0][0] * Z[0] + S2[0][1] * Z[1] + S2[0][2] * Z[2],
                   S2[1][0] * Z[0] + S2[1][1] * Z[1] + S2[1][2] * Z[2],
                   S2[2][0] * Z[0] + S2[2][1] * Z[1] + S2[2][2] * Z[2]};
  double ZtS2Z = Z[0] * S2Z[0] + Z[1] * S2Z[1] + Z[2] * S2Z[2];

  double a = sigma * sigma - trace_adj_S;
  double bcoef = sigma * sigma + ZtZ;
  double c = detS + ZtSZ;
  double dcoef = ZtS2Z;
  double coef_lam2 = -(a + bcoef);
  double coef_lam1 = -c;
  double coef_lam0 = a * bcoef + c * sigma - dcoef;

  double lambda = 1.0;
  bool converged = false;
  for (int it = 0; it < 30; ++it) {
    double f = lambda * lambda * lambda * lambda +
               coef_lam2 * lambda * lambda + coef_lam1 * lambda + coef_lam0;
    double fp = 4.0 * lambda * lambda * lambda + 2.0 * coef_lam2 * lambda +
                coef_lam1;
    if (std::abs(fp) < 1e-30) break;
    double delta = f / fp;
    lambda -= delta;
    if (std::abs(delta) < 1e-14) {
      converged = true;
      break;
    }
  }
  if (!converged || !std::isfinite(lambda)) return false;

  double alpha = lambda * lambda - sigma * sigma + trace_adj_S;
  double beta = lambda - sigma;
  double gamma = (lambda + sigma) * alpha - detS;
  double Mq[3][3];
  for (int aa = 0; aa < 3; ++aa)
    for (int bb = 0; bb < 3; ++bb) {
      Mq[aa][bb] = beta * S[aa][bb] + S2[aa][bb];
      if (aa == bb) Mq[aa][bb] += alpha;
    }
  double X[3] = {Mq[0][0] * Z[0] + Mq[0][1] * Z[1] + Mq[0][2] * Z[2],
                 Mq[1][0] * Z[0] + Mq[1][1] * Z[1] + Mq[1][2] * Z[2],
                 Mq[2][0] * Z[0] + Mq[2][1] * Z[1] + Mq[2][2] * Z[2]};
  double norm_sq = gamma * gamma + X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
  if (!std::isfinite(norm_sq) || norm_sq <= 0.0) return false;
  double inv_n = 1.0 / std::sqrt(norm_sq);
  double qx = X[0] * inv_n, qy = X[1] * inv_n, qz = X[2] * inv_n;
  double qw = gamma * inv_n;

  // Convert (qx, qy, qz, qw) — the inertial->camera quaternion — to a 3x3
  // rotation matrix mapping inertial vectors into the camera frame.
  R_out[0][0] = 1.0 - 2.0 * (qy * qy + qz * qz);
  R_out[0][1] = 2.0 * (qx * qy - qz * qw);
  R_out[0][2] = 2.0 * (qx * qz + qy * qw);
  R_out[1][0] = 2.0 * (qx * qy + qz * qw);
  R_out[1][1] = 1.0 - 2.0 * (qx * qx + qz * qz);
  R_out[1][2] = 2.0 * (qy * qz - qx * qw);
  R_out[2][0] = 2.0 * (qx * qz - qy * qw);
  R_out[2][1] = 2.0 * (qy * qz + qx * qw);
  R_out[2][2] = 1.0 - 2.0 * (qx * qx + qy * qy);
  return true;
}

// Phase 3e.5: replace the TRIAD-on-the-original-4 attitude with a QUEST
// attitude over the full current inlier set, then re-run tight inlier
// expansion using that refined attitude. One pass of refine→re-expand is
// usually enough: the second pass typically adds 0 new inliers.
//
// Mutates `assignment` and writes the refined rotation into `R_out` (which
// also serves as the seed attitude on entry). Returns the final inlier count.
int refine_and_reexpand(std::vector<int> &assignment,
                        const std::vector<std::array<double, 3>> &v_cam,
                        const StarDatabase &db,
                        const std::function<std::array<double, 3>(int)> &cat_vec,
                        double R_out[3][3], double accept_cos_tol) {
  const int N = static_cast<int>(v_cam.size());

  // First pass: expand using the TRIAD seed attitude already in R_out.
  int inliers = expand_inliers_tight(assignment, v_cam, db, cat_vec, R_out,
                                      accept_cos_tol);
  if (inliers < 3) return inliers;

  // QUEST refine on the current set.
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
  if (!quest_attitude_local(v_cam_set, v_iner_set, R_refined)) {
    return inliers;
  }
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = R_refined[r][c];

  // Second pass with refined attitude; this can rescue inliers that were
  // just outside the gate under the loose TRIAD-on-2-vectors seed attitude.
  int new_inliers = expand_inliers_tight(assignment, v_cam, db, cat_vec,
                                          R_out, accept_cos_tol);
  return new_inliers;
}

// Try to verify a candidate 4-tuple pattern match by:
//   1. Pair observed[centroid_idx_canonical[i]] (i in 0..3) with hips[i].
//   2. Solve TRIAD on the pair with widest catalog separation (best numerical
//      conditioning).
//   3. Check all 6 candidate pair cosines are self-consistent against the
//      catalog. If even one pair is grossly off, reject without going to the
//      5th-star step.
//   4. Project each of the unmatched candidate centroids into inertial and
//      verify at least one lands within `verify_cos_tol` of a catalog star
//      (looked up via db.find_partners of one of the matched HIPs). On first
//      match, accept.
//
// Returns the (partial) assignment image_idx -> HIP, with -1 for unassigned,
// or nullopt-equivalent (empty vector) if rejection. On success, `R_out`
// receives the TRIAD-derived inertial->camera rotation so the caller can
// drive the inlier-expansion step without re-solving.
std::vector<int>
try_verify_candidate(const std::array<int, 4> &centroid_idx_canonical,
                     const std::array<int, 4> &hips,
                     const std::vector<int> &verify_pool,
                     const std::vector<std::array<double, 3>> &v_cam,
                     const StarDatabase &db,
                     const std::function<std::array<double, 3>(int)> &cat_vec,
                     double pair_cos_tol, double verify_cos_tol,
                     double R_out[3][3]) {
  const int N = static_cast<int>(v_cam.size());
  std::vector<int> assignment(N, -1);

  std::array<int, 4> centroid_idx = centroid_idx_canonical;
  std::array<std::array<double, 3>, 4> v_obs{};
  std::array<std::array<double, 3>, 4> v_cat{};
  for (int i = 0; i < 4; ++i) {
    v_obs[i] = v_cam[centroid_idx[i]];
    // db.get_star will throw if HIP not in catalog — shouldn't happen for a
    // valid pattern record, but we don't want a throw to leak out.
    try {
      v_cat[i] = cat_vec(hips[i]);
    } catch (...) {
      return {};
    }
  }

  // Pair-consistency gate. The 6 inter-pair cosines should match the catalog
  // within `pair_cos_tol`. This is the same self-consistency check the
  // pyramid does, but with the canonical pairing fixed by the hash so we
  // don't need to permute.
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      double obs = dot3(v_obs[i], v_obs[j]);
      double cat = dot3(v_cat[i], v_cat[j]);
      if (std::abs(obs - cat) > pair_cos_tol) {
        return {}; // reject
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
  triad_rotation(v_obs[best_a], v_obs[best_b], v_cat[best_a], v_cat[best_b], R);
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = R[r][c];

  // Provisional 4-star assignment.
  for (int i = 0; i < 4; ++i) {
    assignment[centroid_idx[i]] = hips[i];
  }

  // Try to verify a 5th star: pick the brightest unmatched centroid in the
  // verification pool, project it into the inertial frame, and look up its
  // catalog neighbor against any one of the 4 matched HIPs.
  //
  // The pool is expected to be ordered brightest-first, so this gives the
  // strongest unmatched centroid the first chance.
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

    // Look up candidate HIPs at the observed angle from matched HIP 0.
    // (Matched HIP 0 is the brightest star of the 4 in canonical order in
    // generic cases; either way the geometry is the same.)
    int anchor_hip = hips[0];
    auto v_anchor = v_cat[0];
    double cos_pred = dot3(v_anchor, v_iner_5);
    if (cos_pred < -1.0) cos_pred = -1.0;
    if (cos_pred > 1.0) cos_pred = 1.0;
    // The partner search is exact-cosine-distance; tolerance covers both
    // FOV miscal and centroid noise. Use a moderately loose tolerance so
    // we don't miss the correct partner under realistic noise — final
    // gating is by angular residual below.
    auto candidates = db.find_partners(anchor_hip, cos_pred, 1e-3);
    int best_C = -1;
    double best_dot = -2.0; // largest dot = closest match
    for (int C : candidates) {
      // Skip the 4 already-matched HIPs.
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
    return {}; // reject — 4 stars alone are too weak to commit to a candidate
  }

  assignment[fifth_centroid] = fifth_hip;
  return assignment;
}

// === Phase 3e.3: pattern-path identify_stars implementation. ===
//
// Returns the (filled) assignment image_idx -> HIP, or an empty vector if no
// pattern attempt produced a verified attitude. On success, the assignment
// covers the 4 pattern stars + the verified 5th star (other stars may still
// be -1; caller / expansion step extends from there).
//
// Logs "[3e] pattern path failed → pyramid fallback" if no candidate
// verifies after the first kPatternAttempts seed centroids — that lets us
// count fallback frequency in production.

// Number of seed centroids we try (brightest N by peak). Originally 4 per
// plan 3e.3. Diagnostics in the 3e.5 follow-up (DEBUG_PATTERN_STATS) showed
// the failure mode on noisy synthetic scenes: when one of the top-4 brightest
// centroids isn't a catalog star (or its nearest detected neighbors aren't a
// subset of its catalog top-8 neighbors), all 4 seeds can fail in lock-step
// even though dimmer-seed alternatives would succeed. Bumping to 10 raises
// the per-frame pattern-path hit rate substantially on Monte Carlo without
// measurably increasing the success-case latency (the path returns on first
// verify, so extra seeds are only exercised on the fallback boundary).
constexpr int kPatternAttempts = 10;
// Per-seed, enumerate C(kPatternKNearest, 3) triples drawn from the seed's
// kPatternKNearest nearest in-radius neighbors. The catalog generator stores
// patterns built from each catalog star's top-8 nearest within FOV/2, so
// kPatternKNearest=8 mirrors the catalog's coverage. We bump to 10 to
// absorb cases where centroid noise / vmag-cutoff shifts which neighbor
// ranks highest.
constexpr int kPatternKNearest = 10;
// Verification pool: where the 5th star is hunted. The seed's 3 nearest
// neighbors are drawn from ALL detected centroids (not just bright ones) —
// the catalog only contains 4-tuples within FOV/2 of the seed, and at wide
// FOVs the brightest stars are sparse enough that a bright-only pool can
// be empty inside that radius. The verify pool is intentionally narrower
// because verification should anchor on real stars (not noise).
constexpr int kVerifyPoolSize = 12;
// 5th-star angular acceptance threshold ≈ 0.05° per the plan. cos(0.05°) ≈
// 1 - 3.8e-7, so we require dot >= 1 - 4e-7.
constexpr double kFifthStarVerifyCosTol = 4e-7;
// Inlier-expansion acceptance threshold. Same 0.05° gate the brief mandates
// for the post-verify inlier expansion step (cos(0.05°) ≈ 1 - 3.8e-7). Wide
// enough to absorb sub-degree FOV miscal that the TRIAD seed attitude can't
// correct, while still rejecting wrong-star matches (which would land 10×
// further). The follow-up FOV-scale absorption in identify_stars() tightens
// the implied tolerance once the camera focal is rescaled.
constexpr double kInlierExpansionCosTol = 4e-7;
// Pair-consistency cosine tolerance for the 4 pattern stars. Cos noise
// across a 5-20° pattern at our centroid noise level is ~1e-5 — we use 5e-4
// which still rejects gross mismatches but is loose enough to absorb FOV
// calibration drift (1e-3 was the alt60 worst case).
constexpr double kPatternPairCosTol = 5e-4;
// Safety margin on the catalog FOV radius: the catalog stores patterns
// within FOV_bin/2 of each seed. We don't know FOV_bin directly (catalog.h
// doesn't expose it), but the binary picks FOV_bin to bracket the camera
// FOV. For real-image fixtures the camera FOV is ~11° and bin=10° (radius
// 5°); for FOV=20° MC bin=20° (radius 10°). cos(half_FOV) is a safe upper
// bound on the catalog radius — and any quad of stars within that radius is
// either in the catalog or geometrically equivalent to one that is. We pull
// in a small margin (95%) to absorb FOV-bin/camera-FOV mismatch.
constexpr double kPatternRadiusSafetyMargin = 0.95;

std::vector<int>
identify_stars_pattern(const std::vector<StarCentroid> &image_stars,
                       const PinholeCamera &camera,
                       const std::vector<std::array<double, 3>> &v_cam,
                       const StarDatabase &db,
                       const std::function<std::array<double, 3>(int)> &cat_vec) {
  const int N = static_cast<int>(image_stars.size());
  if (N < 5)
    return {}; // pattern path requires 4 + 1 verification star

  // Pattern catalogs only contain 4-star tuples with all 4 stars within
  // FOV_bin/2 of the brightest. Use the camera half-FOV as a proxy for
  // FOV_bin/2 (main.cpp picks FOV_bin to bracket the camera FOV).
  const double half_fov_rad = camera_half_fov_rad(camera) *
                              kPatternRadiusSafetyMargin;
  const double cos_radius = std::cos(half_fov_rad);

  // Rank all centroids by peak descending. The first kSeedPoolSize are seed
  // candidates; the full ranked list is the neighbor pool (the seed's 3
  // nearest within FOV/2 may be a mag-5 star next to the bright seed, not
  // necessarily another mag-2 star). Verification pool is the top-N
  // brightest so we don't anchor 5th-star verification on faint noise.
  auto ranked = top_n_indices(image_stars, static_cast<int>(image_stars.size()));
  if (static_cast<int>(ranked.size()) < 4) return {};

  // Neighbor pool = ALL detected centroids (peak-ordered, brightest first).
  const std::vector<int> &neighbor_pool = ranked;
  std::vector<int> verify_pool(
      ranked.begin(),
      ranked.begin() + std::min<int>(ranked.size(), kVerifyPoolSize));

  const int seeds = std::min(kPatternAttempts, static_cast<int>(ranked.size()));
  for (int s = 0; s < seeds; ++s) {
    int seed_idx = ranked[s];

    // Find the seed's nearest in-radius neighbors. The catalog stores
    // patterns built from each catalog star's top-8 nearest within FOV/2;
    // mirror that on the observation side by enumerating triples of the
    // seed's kPatternKNearest nearest in-radius detected centroids.
    auto knn = k_nearest_within_radius(seed_idx, neighbor_pool, v_cam,
                                       cos_radius, kPatternKNearest);
    if (knn.size() < 3) continue; // not enough neighbors → next seed

    // Enumerate C(knn.size(), 3) triples. With kPatternKNearest=10 this is
    // up to 120 hashes per seed, but most catalog probes return zero or
    // few candidates — the hot path is the hash + tolerant probe, both
    // cheap.
    const int K = static_cast<int>(knn.size());
    for (int a = 0; a < K - 2; ++a) {
      for (int b = a + 1; b < K - 1; ++b) {
        for (int c = b + 1; c < K; ++c) {
          std::array<int, 4> centroid_indices = {seed_idx, knn[a], knn[b],
                                                  knn[c]};

          // Phase 3e.5 (Change 1): noise-robust multi-key probing.
          //
          // The canonical-order key sorts the 6 inter-star angles. Under
          // centroid noise, two near-equal angles can flip rank, shifting the
          // canonical permutation and bumping the key to a distant bucket
          // that ±1 tolerant probing can't reach. The catalog only stores the
          // canonical key for the noise-free geometry, so we enumerate the
          // 2^k orderings produced by flipping each uncertain adjacent-rank
          // pair (where k is the number of uncertain pairs ≤ 5; in practice
          // 1–3). One of those orderings will match the catalog's canonical
          // key. Probe each with find_pattern_tolerant for residual ±1
          // quantization slack.
          //
          // The 24-input-permutation variant tetra3 uses doesn't apply to our
          // signature-based canonical_order: input permutation is invariant
          // by construction here (per-star signatures are geometry-only, and
          // both catalog and query tie-break on the input identifier — HIP
          // and centroid index respectively — so all 24 input permutations
          // produce the same canonical key). The rank-flip mechanism in
          // pattern_keys_noise_robust is what actually addresses the noise
          // failure mode.
          //
          // noise_tol = 2 mrad ≈ 7 arcmin. At synthetic 5-px-noise input on a
          // 5000 px focal, per-pair angular noise is ~1.4e-3 rad; 1.5σ ≈ 2e-3
          // catches the common rank-flip cases without generating an
          // unreasonably large alternate-key set. Real images (lower per-star
          // noise) end up with the same key as the canonical, which is the
          // single-key path.
          std::array<std::array<double, 3>, 4> vecs4{};
          std::array<int, 4> ids4{};
          for (int i = 0; i < 4; ++i) {
            vecs4[i] = v_cam[centroid_indices[i]];
            ids4[i] = centroid_indices[i];
          }
          constexpr double kPatternRankFlipNoiseRad = 2e-3;
          auto raw_keys =
              pattern_keys_noise_robust(vecs4, ids4, kPatternRankFlipNoiseRad);
          std::vector<std::pair<uint64_t, std::array<int, 4>>> probe_keys;
          probe_keys.reserve(raw_keys.size());
          for (const auto &kp : raw_keys) {
            std::array<int, 4> centroid_canonical{};
            for (int i = 0; i < 4; ++i)
              centroid_canonical[i] = centroid_indices[kp.second[i]];
            probe_keys.emplace_back(kp.first, centroid_canonical);
          }
          for (const auto &kp : probe_keys) {
            uint64_t key = kp.first;
            const auto &centroid_canonical = kp.second;
            auto candidates = db.find_pattern_tolerant(key);
            if (candidates.empty()) continue;

            for (const auto &cand : candidates) {
              std::array<int, 4> hips = {cand.hips[0], cand.hips[1],
                                          cand.hips[2], cand.hips[3]};
              double R_triad[3][3] = {{0}};
              auto assignment = try_verify_candidate(
                  centroid_canonical, hips, verify_pool, v_cam, db, cat_vec,
                  kPatternPairCosTol, kFifthStarVerifyCosTol, R_triad);
              if (assignment.empty()) continue;

              // Change 2: tight inlier expansion + QUEST refine. Projects all
              // remaining centroids into inertial using the TRIAD attitude,
              // matches them against the catalog with a 0.05° gate, then
              // refines the attitude via QUEST and re-expands. Bumps the
              // inlier list from 5 (TRIAD on 4 + 5th-star verify) to N
              // (typically 8–25 on a real scene), which lets the downstream
              // QUEST in estimation.cpp converge to ~0.001° accuracy.
              int n_before = 0;
              for (int ii = 0; ii < (int)assignment.size(); ++ii)
                if (assignment[ii] >= 0) ++n_before;
              int n_after = refine_and_reexpand(assignment, v_cam, db, cat_vec,
                                                R_triad, kInlierExpansionCosTol);
              if (std::getenv("STARTRACKER_DEBUG_3E5")) {
                std::fprintf(stderr,
                             "[3e.5] pattern verify: %d -> %d inliers\n",
                             n_before, n_after);
                for (int ii = 0; ii < (int)assignment.size(); ++ii) {
                  if (assignment[ii] >= 0)
                    std::fprintf(stderr, "  centroid %d -> HIP %d\n", ii,
                                 assignment[ii]);
                }
              }
              return assignment;
            }
          }
        }
      }
    }
  }

  return {};
}

// Original pyramid-path implementation (kept verbatim under a new name so we
// can fall back to it on pattern-path failure or when the catalog isn't
// loaded). The body is the pre-3e identify_stars logic.
std::vector<IdentifiedStar>
identify_stars_pyramid(const std::vector<StarCentroid> &image_stars,
                       const PinholeCamera &camera, const StarDatabase &db,
                       double cos_tolerance,
                       std::vector<std::array<double, 3>> v_cam,
                       const std::function<std::array<double, 3>(int)> &cat_vec) {
  int N = image_stars.size();

  // Threshold for declaring a result acceptable. Both passes share it.
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

} // namespace

std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance) {
  int N = image_stars.size();
  if (N < 3)
    return {};

  // === Camera-frame unit vectors ===
  std::vector<std::array<double, 3>> v_cam = compute_v_cam(image_stars, camera);

  // Cache catalog unit vectors (one per HIP) to avoid repeated map lookups.
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

  // === Phase 3e.3: pattern-hash path ===
  //
  // Try the pattern path first if the catalog is loaded. Pattern path is
  // O(K * candidates_per_key) where K ≤ 4 seed attempts × 243 probes; each
  // hit is one TRIAD + a partner search. Compared to the pyramid's O(N²)
  // seed enumeration over N≤50 centroids, this is ~3 orders of magnitude
  // faster when it succeeds. Falls back to pyramid on miss.
  if (db.has_pattern_catalog() && N >= 5) {
    auto pattern_assignment =
        identify_stars_pattern(image_stars, camera, v_cam, db, cat_vec);
    if (!pattern_assignment.empty()) {
      // The pattern path already performed tight inlier expansion + QUEST
      // refine + re-expand inside identify_stars_pattern (Phase 3e.5 change
      // 2). All inliers in the returned assignment are within 0.05° of their
      // catalog star against the refined attitude, which is tighter than
      // cross_verify's 3*cos_tolerance gate; running cross_verify here would
      // over-reject on miscalibrated cameras.

      // Phase 3e.5: FOV-scale absorption. Mirrors the pyramid path's
      // coarse-refine-reidentify step (3b.0b). The inlier set gives a stable
      // estimate of the focal scale factor s = true_focal / assumed_focal.
      // Real cameras commonly drift 0.1–0.5% from their nominal FOV; the
      // pattern path's TRIAD/QUEST attitude inherits that drift, capping
      // accuracy at ~(s-1) × FOV/2. Rescaling and re-expanding lets the
      // attitude break through that floor and reach centroid-noise limits.
      //
      // Apply iteratively: each pass collapses ~90% of the residual scale
      // error, so 2–3 passes converge to within the noise floor. We stop
      // once the per-pass adjustment falls below the centroid-noise floor
      // (~1e-4 ratio) or we hit the iteration cap.
      CameraModel refined = camera;
      double cumulative_scale = 1.0;
      for (int iter = 0; iter < 3; ++iter) {
        double scale = estimate_scale_factor(pattern_assignment, v_cam, db,
                                              cat_vec);
        if (!std::isfinite(scale)) break;
        // Trigger only when the per-pass residual is well above per-star
        // centroid noise (~1e-4 for our focal). Below that we'd just inject
        // bias rather than remove it.
        if (std::abs(scale - 1.0) <= 1e-4) break;

        refined.focal_x *= scale;
        refined.focal_y *= scale;
        cumulative_scale *= scale;
        auto v_cam_refined = compute_v_cam(image_stars, refined);

        // Re-expand inliers under the refined camera. Solve TRIAD on the two
        // matched stars with the widest catalog separation, then run the
        // tight-expansion + QUEST-refine loop. The same routine that built
        // the initial inlier set is reused; under a properly scaled camera,
        // it typically converges to all visible cataloged stars.
        std::vector<int> refined_assignment(N, -1);
        int n_matched = 0;
        for (int i = 0; i < N; ++i) {
          if (pattern_assignment[i] < 0) continue;
          refined_assignment[i] = pattern_assignment[i];
          ++n_matched;
        }
        if (n_matched < 2) break;
        // Build a seed TRIAD attitude from the two matched stars with the
        // widest catalog separation (best conditioning).
        int best_i = -1, best_j = -1;
        double smallest_cos = 2.0;
        std::vector<int> matched_idx;
        matched_idx.reserve(n_matched);
        for (int i = 0; i < N; ++i)
          if (refined_assignment[i] >= 0) matched_idx.push_back(i);
        for (size_t ai = 0; ai < matched_idx.size(); ++ai) {
          auto vi = cat_vec(refined_assignment[matched_idx[ai]]);
          for (size_t aj = ai + 1; aj < matched_idx.size(); ++aj) {
            auto vj = cat_vec(refined_assignment[matched_idx[aj]]);
            double cc = vi[0] * vj[0] + vi[1] * vj[1] + vi[2] * vj[2];
            if (cc < smallest_cos) {
              smallest_cos = cc;
              best_i = matched_idx[ai];
              best_j = matched_idx[aj];
            }
          }
        }
        if (best_i < 0 || best_j < 0) break;
        std::array<double, 3> W1 = v_cam_refined[best_i];
        std::array<double, 3> W2 = v_cam_refined[best_j];
        std::array<double, 3> V1 = cat_vec(refined_assignment[best_i]);
        std::array<double, 3> V2 = cat_vec(refined_assignment[best_j]);
        double R_seed[3][3];
        triad_rotation(W1, W2, V1, V2, R_seed);
        refine_and_reexpand(refined_assignment, v_cam_refined, db, cat_vec,
                             R_seed, kFifthStarVerifyCosTol);

        // Count inliers under refined camera; commit to refined if it didn't
        // strictly lose ground.
        int refined_n = 0;
        for (int i = 0; i < N; ++i)
          if (refined_assignment[i] >= 0) ++refined_n;
        if (refined_n < n_matched) break; // refinement removed inliers — back out
        if (std::getenv("STARTRACKER_DEBUG_3E5")) {
          std::fprintf(stderr,
                       "[3e.5] FOV iter %d: s=%.6f (cumulative %.6f), "
                       "%d -> %d inliers\n",
                       iter, scale, cumulative_scale, n_matched, refined_n);
        }
        pattern_assignment = std::move(refined_assignment);
        v_cam = std::move(v_cam_refined);
      }

      int inliers = 0;
      for (int i = 0; i < N; ++i)
        if (pattern_assignment[i] >= 0) ++inliers;
      if (inliers >= 4) {
        std::vector<IdentifiedStar> identified;
        identified.reserve(static_cast<size_t>(inliers));
        for (int i = 0; i < N; ++i) {
          if (pattern_assignment[i] < 0) continue;
          IdentifiedStar is;
          is.image_idx = i;
          is.catalog_hip_id = pattern_assignment[i];
          is.v_cam[0] = v_cam[i][0];
          is.v_cam[1] = v_cam[i][1];
          is.v_cam[2] = v_cam[i][2];
          identified.push_back(is);
        }
        return identified;
      }
      // Fall through to pyramid — pattern path matched but produced too
      // few inliers; treat as a soft fallback.
    }
    std::fprintf(stderr, "[3e] pattern path failed → pyramid fallback\n");
  }

  // === Fallback: pyramid + coarse-refine-reidentify + cross-verify ===
  return identify_stars_pyramid(image_stars, camera, db, cos_tolerance,
                                std::move(v_cam), cat_vec);
}
