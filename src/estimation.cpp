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

// TRIAD fallback. Solves Wahba on the two stars with the largest angular
// separation among the supplied (already pair-validated) candidates. The
// `best_i` / `best_j` indices are produced by the pair-validation loop in
// `estimate_attitude` so we don't re-walk the O(N^2) check here.
static Quaternion estimate_attitude_triad(const std::vector<IdentifiedStar> &stars,
                                          const StarDatabase &db, int best_i,
                                          int best_j) {
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
  // Build attitude profile matrix B = Σ w_i · b_i · r_i^T
  // b_i = camera-frame unit vector, r_i = inertial unit vector.
  // Equal weights w_i = 1/M (so Σ w_i = 1, which is the noise-free max eigenvalue).
  const int M = (int)valid_idx.size();
  const double w = 1.0 / (double)M;

  double B[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int k = 0; k < M; ++k) {
    int i = valid_idx[k];
    double b[3] = {stars[i].v_cam[0], stars[i].v_cam[1], stars[i].v_cam[2]};
    CatalogStar c = db.get_star(stars[i].catalog_hip_id);
    double r[3] = {c.x, c.y, c.z};
    for (int a = 0; a < 3; ++a) {
      for (int bb = 0; bb < 3; ++bb) {
        B[a][bb] += w * b[a] * r[bb];
      }
    }
  }

  // S = B + B^T
  double S[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int bb = 0; bb < 3; ++bb) {
      S[a][bb] = B[a][bb] + B[bb][a];
    }
  }

  // σ = trace(B)
  double sigma = B[0][0] + B[1][1] + B[2][2];

  // Z = vex(B - B^T) = [B21-B12, B02-B20, B10-B01]^T (Shuster's convention so
  // that the resulting optimal quaternion's vector part aligns with the axis
  // of rotation; verified on the IdentityRotation / Known90DegZ tests).
  double Z[3] = {B[2][1] - B[1][2], B[0][2] - B[2][0], B[1][0] - B[0][1]};

  // Helpers for 3x3 matrix algebra
  auto det3 = [](const double M3[3][3]) {
    return M3[0][0] * (M3[1][1] * M3[2][2] - M3[1][2] * M3[2][1]) -
           M3[0][1] * (M3[1][0] * M3[2][2] - M3[1][2] * M3[2][0]) +
           M3[0][2] * (M3[1][0] * M3[2][1] - M3[1][1] * M3[2][0]);
  };

  // trace(adj(S)) = sum of 2x2 principal minors of S
  double trace_adj_S = (S[1][1] * S[2][2] - S[1][2] * S[2][1]) +
                       (S[0][0] * S[2][2] - S[0][2] * S[2][0]) +
                       (S[0][0] * S[1][1] - S[0][1] * S[1][0]);

  double detS = det3(S);

  // Z^T Z
  double ZtZ = Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2];

  // S*Z
  double SZ[3] = {
      S[0][0] * Z[0] + S[0][1] * Z[1] + S[0][2] * Z[2],
      S[1][0] * Z[0] + S[1][1] * Z[1] + S[1][2] * Z[2],
      S[2][0] * Z[0] + S[2][1] * Z[1] + S[2][2] * Z[2],
  };

  // Z^T S Z
  double ZtSZ = Z[0] * SZ[0] + Z[1] * SZ[1] + Z[2] * SZ[2];

  // S^2
  double S2[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int bb = 0; bb < 3; ++bb) {
      S2[a][bb] = S[a][0] * S[0][bb] + S[a][1] * S[1][bb] + S[a][2] * S[2][bb];
    }
  }

  // S^2 * Z
  double S2Z[3] = {
      S2[0][0] * Z[0] + S2[0][1] * Z[1] + S2[0][2] * Z[2],
      S2[1][0] * Z[0] + S2[1][1] * Z[1] + S2[1][2] * Z[2],
      S2[2][0] * Z[0] + S2[2][1] * Z[1] + S2[2][2] * Z[2],
  };

  // Z^T S^2 Z
  double ZtS2Z = Z[0] * S2Z[0] + Z[1] * S2Z[1] + Z[2] * S2Z[2];

  // Characteristic-equation coefficients (Shuster):
  //   λ^4 - (a + b) λ^2 - c λ + (a·b + c·σ - d) = 0
  double a = sigma * sigma - trace_adj_S;
  double b = sigma * sigma + ZtZ;
  double c = detS + ZtSZ;
  double d = ZtS2Z;

  double coef_lam2 = -(a + b);
  double coef_lam1 = -c;
  double coef_lam0 = a * b + c * sigma - d;

  // Newton-Raphson for λ_max. Seed at Σ w_i = 1 (the noise-free eigenvalue).
  double lambda = 1.0;
  bool converged = false;
  for (int it = 0; it < 20; ++it) {
    double f = lambda * lambda * lambda * lambda +
               coef_lam2 * lambda * lambda + coef_lam1 * lambda + coef_lam0;
    double fp = 4.0 * lambda * lambda * lambda + 2.0 * coef_lam2 * lambda +
                coef_lam1;
    if (std::abs(fp) < 1e-30) {
      break; // singular derivative
    }
    double delta = f / fp;
    lambda -= delta;
    if (std::abs(delta) < 1e-12) {
      converged = true;
      break;
    }
  }

  if (!converged || !std::isfinite(lambda)) {
    std::cerr << "QUEST: fallback to TRIAD\n";
    return estimate_attitude_triad(stars, db, best_i, best_j);
  }

  // Shuster's closed-form optimal quaternion:
  //   α = λ^2 - σ^2 + trace(adj(S))
  //   β = λ - σ
  //   γ = (λ + σ)·α - det(S)
  //   X = (α·I + β·S + S^2) · Z
  //   q_xyz = X,  q_w = γ,  then normalize.
  double alpha = lambda * lambda - sigma * sigma + trace_adj_S;
  double beta = lambda - sigma;
  double gamma = (lambda + sigma) * alpha - detS;

  // M = α·I + β·S + S^2
  double Mq[3][3];
  for (int aa = 0; aa < 3; ++aa) {
    for (int bb = 0; bb < 3; ++bb) {
      Mq[aa][bb] = beta * S[aa][bb] + S2[aa][bb];
      if (aa == bb)
        Mq[aa][bb] += alpha;
    }
  }

  double X[3] = {
      Mq[0][0] * Z[0] + Mq[0][1] * Z[1] + Mq[0][2] * Z[2],
      Mq[1][0] * Z[0] + Mq[1][1] * Z[1] + Mq[1][2] * Z[2],
      Mq[2][0] * Z[0] + Mq[2][1] * Z[1] + Mq[2][2] * Z[2],
  };

  double norm_sq = gamma * gamma + X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
  if (!std::isfinite(norm_sq) || norm_sq <= 0.0) {
    std::cerr << "QUEST: fallback to TRIAD\n";
    return estimate_attitude_triad(stars, db, best_i, best_j);
  }
  double inv_n = 1.0 / std::sqrt(norm_sq);

  Quaternion q;
  q.x = X[0] * inv_n;
  q.y = X[1] * inv_n;
  q.z = X[2] * inv_n;
  q.w = gamma * inv_n;

  // Canonicalize sign so q.w >= 0 (matches TRIAD output convention on the
  // IdentityRotation test where the expected quaternion is [0,0,0,+1]).
  if (q.w < 0.0) {
    q.x = -q.x;
    q.y = -q.y;
    q.z = -q.z;
    q.w = -q.w;
  }

  if (!std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z) ||
      !std::isfinite(q.w)) {
    std::cerr << "QUEST: fallback to TRIAD\n";
    return estimate_attitude_triad(stars, db, best_i, best_j);
  }

  return q;
}
