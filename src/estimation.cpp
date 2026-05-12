#include "estimation.h"
#include "wahba.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

struct Vec3 {
  double x, y, z;
};

inline double dot(const Vec3 &a, const Vec3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

// TRIAD fallback. Solves Wahba on the two stars with the largest angular
// separation among the supplied (already pair-validated) candidates. The
// `best_i` / `best_j` indices are produced by the pair-validation loop in
// `estimate_attitude` so we don't re-walk the O(N^2) check here.
//
// Phase 3g.1: TRIAD algebra moved to wahba::triad; here we just shuffle types
// and convert the resulting rotation matrix to a Quaternion in the project's
// public convention (x, y, z, w; w >= 0).
static Quaternion estimate_attitude_triad(const std::vector<IdentifiedStar> &stars,
                                          const StarDatabase &db, int best_i,
                                          int best_j) {
  const std::array<double, 3> W1 = {stars[best_i].v_cam[0],
                                    stars[best_i].v_cam[1],
                                    stars[best_i].v_cam[2]};
  const std::array<double, 3> W2 = {stars[best_j].v_cam[0],
                                    stars[best_j].v_cam[1],
                                    stars[best_j].v_cam[2]};
  const CatalogStar cat1 = db.get_star(stars[best_i].catalog_hip_id);
  const CatalogStar cat2 = db.get_star(stars[best_j].catalog_hip_id);
  const std::array<double, 3> V1 = {cat1.x, cat1.y, cat1.z};
  const std::array<double, 3> V2 = {cat2.x, cat2.y, cat2.z};

  double R[3][3];
  wahba::triad(W1, W2, V1, V2, R);

  // Rotation matrix -> quaternion (x, y, z, w). Standard 4-branch
  // numerically-stable formula; preserves the (qw >= 0) convention used by
  // the IdentityRotation test.
  double tr = R[0][0] + R[1][1] + R[2][2];
  Quaternion q;
  if (tr > 0) {
    double S = std::sqrt(tr + 1.0) * 2; // S=4*qw
    q.w = 0.25 * S;
    q.x = (R[2][1] - R[1][2]) / S;
    q.y = (R[0][2] - R[2][0]) / S;
    q.z = (R[1][0] - R[0][1]) / S;
  } else if ((R[0][0] > R[1][1]) && (R[0][0] > R[2][2])) {
    double S = std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2; // S=4*qx
    q.w = (R[2][1] - R[1][2]) / S;
    q.x = 0.25 * S;
    q.y = (R[0][1] + R[1][0]) / S;
    q.z = (R[0][2] + R[2][0]) / S;
  } else if (R[1][1] > R[2][2]) {
    double S = std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2; // S=4*qy
    q.w = (R[0][2] - R[2][0]) / S;
    q.x = (R[0][1] + R[1][0]) / S;
    q.y = 0.25 * S;
    q.z = (R[1][2] + R[2][1]) / S;
  } else {
    double S = std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2; // S=4*qz
    q.w = (R[1][0] - R[0][1]) / S;
    q.x = (R[0][2] + R[2][0]) / S;
    q.y = (R[1][2] + R[2][1]) / S;
    q.z = 0.25 * S;
  }
  double norm = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
  q.x /= norm;
  q.y /= norm;
  q.z /= norm;
  q.w /= norm;
  return q;
}

Quaternion estimate_attitude(const std::vector<IdentifiedStar> &stars,
                             const StarDatabase &db) {
  if (stars.size() < 2) {
    throw std::runtime_error("Need at least 2 stars for attitude estimation");
  }

  // Pair-validation loop: rejects obviously misidentified stars by checking
  // |dot(W_i, W_j) - dot(V_i, V_j)| against ~3e-3. Collect the set of "valid"
  // stars (any star appearing in at least one consistent pair) and also track
  // the widest-separation valid pair as a TRIAD fallback anchor.
  const int N = (int)stars.size();
  std::vector<bool> in_valid_pair(N, false);
  int best_i = -1, best_j = -1;
  double max_sep = -1.0;
  double best_diff = 1e9;
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      Vec3 w_i = {stars[i].v_cam[0], stars[i].v_cam[1], stars[i].v_cam[2]};
      Vec3 w_j = {stars[j].v_cam[0], stars[j].v_cam[1], stars[j].v_cam[2]};
      double dW = dot(w_i, w_j);

      CatalogStar c_i = db.get_star(stars[i].catalog_hip_id);
      CatalogStar c_j = db.get_star(stars[j].catalog_hip_id);
      Vec3 v_i = {c_i.x, c_i.y, c_i.z};
      Vec3 v_j = {c_j.x, c_j.y, c_j.z};
      double dV = dot(v_i, v_j);

      double diff = std::abs(dW - dV);
      if (diff < best_diff)
        best_diff = diff;

      // 3e-3 cos tolerance ≈ 0.17° – wide enough for real centroid noise while
      // still rejecting badly misidentified stars (whose diffs would be >> 0.01)
      if (diff < 3e-3) {
        in_valid_pair[i] = true;
        in_valid_pair[j] = true;
        double sep = 1.0 - dW;
        if (sep > max_sep) {
          max_sep = sep;
          best_i = i;
          best_j = j;
        }
      }
    }
  }

  if (best_i == -1 || best_j == -1) {
    std::cerr << "QUEST/TRIAD: no valid pair found. Best |dW-dV| across "
              << stars.size() << " stars = " << best_diff << "\n";
    throw std::runtime_error("No valid pairs found for attitude estimation");
  }

  // Build the QUEST input set from stars that appear in a valid pair.
  std::vector<int> valid_idx;
  valid_idx.reserve(N);
  for (int i = 0; i < N; ++i) {
    if (in_valid_pair[i])
      valid_idx.push_back(i);
  }

  // Need at least 3 valid stars for QUEST to be worth running; otherwise TRIAD
  // is the optimal 2-star solution.
  if ((int)valid_idx.size() < 3) {
    return estimate_attitude_triad(stars, db, best_i, best_j);
  }

  // --- QUEST (Shuster & Oh 1981) ---
  // Phase 3g.1: numerical core moved to wahba::quest. The public path keeps
  // 20 Newton iters + 1e-12 tolerance (was the case before the extraction);
  // the identification-side refine pass uses 30 / 1e-14 because its post-
  // expansion gates are tighter.
  std::vector<std::array<double, 3>> v_cam_set;
  std::vector<std::array<double, 3>> v_iner_set;
  v_cam_set.reserve(valid_idx.size());
  v_iner_set.reserve(valid_idx.size());
  for (int i : valid_idx) {
    v_cam_set.push_back(
        {stars[i].v_cam[0], stars[i].v_cam[1], stars[i].v_cam[2]});
    CatalogStar c = db.get_star(stars[i].catalog_hip_id);
    v_iner_set.push_back({c.x, c.y, c.z});
  }

  double q_xyzw[4];
  if (!wahba::quest(v_cam_set, v_iner_set, q_xyzw)) {
    std::cerr << "QUEST: fallback to TRIAD\n";
    return estimate_attitude_triad(stars, db, best_i, best_j);
  }

  Quaternion q;
  q.x = q_xyzw[0];
  q.y = q_xyzw[1];
  q.z = q_xyzw[2];
  q.w = q_xyzw[3]; // wahba::quest already canonicalises w >= 0
  return q;
}
