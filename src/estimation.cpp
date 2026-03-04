#include "estimation.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

struct Vec3 {
  double x, y, z;
};

inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double dot(const Vec3 &a, const Vec3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 normalize(const Vec3 &v) {
  double mag = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  return {v.x / mag, v.y / mag, v.z / mag};
}

Quaternion estimate_attitude(const std::vector<IdentifiedStar> &stars,
                             const StarDatabase &db) {
  if (stars.size() < 2) {
    throw std::runtime_error("Need at least 2 stars for TRIAD");
  }

  // Pick the two stars with the largest angular separation that ALSO form a
  // valid pair (camera dot product consistent with catalog dot product).
  // Tolerance of 1e-3 accounts for ~0.05 deg centroid noise.
  int best_i = -1, best_j = -1;
  double max_sep = -1.0;
  double best_diff = 1e9;
  for (int i = 0; i < (int)stars.size(); ++i) {
    for (int j = i + 1; j < (int)stars.size(); ++j) {
      // Camera vectors
      Vec3 w_i = {stars[i].v_cam[0], stars[i].v_cam[1], stars[i].v_cam[2]};
      Vec3 w_j = {stars[j].v_cam[0], stars[j].v_cam[1], stars[j].v_cam[2]};
      double dW = dot(w_i, w_j);

      // Inertial vectors
      CatalogStar c_i = db.get_star(stars[i].catalog_hip_id);
      CatalogStar c_j = db.get_star(stars[j].catalog_hip_id);
      Vec3 v_i = {c_i.x, c_i.y, c_i.z};
      Vec3 v_j = {c_j.x, c_j.y, c_j.z};
      double dV = dot(v_i, v_j);

      double diff = std::abs(dW - dV);
      if (diff < best_diff)
        best_diff = diff; // track for diagnostics

      // 3e-3 cos tolerance ≈ 0.17° – wide enough for real centroid noise while
      // still rejecting badly misidentified stars (whose diffs would be >> 0.01)
      if (diff < 3e-3) {
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
    std::cerr << "TRIAD: no valid pair found. Best |dW-dV| across "
              << stars.size() << " stars = " << best_diff << "\n";
    throw std::runtime_error("No valid pairs found for TRIAD");
  }

  // Camera vectors (W)
  Vec3 W1 = {stars[best_i].v_cam[0], stars[best_i].v_cam[1],
             stars[best_i].v_cam[2]};
  Vec3 W2 = {stars[best_j].v_cam[0], stars[best_j].v_cam[1],
             stars[best_j].v_cam[2]};

  // Inertial vectors (V) from catalog
  CatalogStar cat1 = db.get_star(stars[best_i].catalog_hip_id);
  CatalogStar cat2 = db.get_star(stars[best_j].catalog_hip_id);

  Vec3 V1 = {cat1.x, cat1.y, cat1.z};
  Vec3 V2 = {cat2.x, cat2.y, cat2.z};

  // TRIAD for Camera frame (W)
  Vec3 t1W = W1;
  Vec3 t2W = normalize(cross(W1, W2));
  Vec3 t3W = cross(t1W, t2W);

  // TRIAD for Inertial frame (V)
  Vec3 t1V = V1;
  Vec3 t2V = normalize(cross(V1, V2));
  Vec3 t3V = cross(t1V, t2V);

  // R = M_W * M_V^T
  double R[3][3];
  double M_W[3][3] = {
      {t1W.x, t2W.x, t3W.x}, {t1W.y, t2W.y, t3W.y}, {t1W.z, t2W.z, t3W.z}};

  double M_V_T[3][3] = {
      {t1V.x, t1V.y, t1V.z}, {t2V.x, t2V.y, t2V.z}, {t3V.x, t3V.y, t3V.z}};

  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      R[r][c] = M_W[r][0] * M_V_T[0][c] + M_W[r][1] * M_V_T[1][c] +
                M_W[r][2] * M_V_T[2][c];
    }
  }

  // Convert R (Rotation Matrix) to Quaternion [x,y,z,w]
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

  // Normalize to be safe
  double norm = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
  q.x /= norm;
  q.y /= norm;
  q.z /= norm;
  q.w /= norm;

  return q;
}
