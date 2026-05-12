#include "wahba.h"
#include <cmath>

namespace wahba {

void triad(const std::array<double, 3> &W1,
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

  // R = M_W * M_V^T sends V-frame vectors to W-frame: v_W = R · v_V.
  const double M_W[3][3] = {{t1W[0], t2W[0], t3W[0]},
                            {t1W[1], t2W[1], t3W[1]},
                            {t1W[2], t2W[2], t3W[2]}};
  const double M_V_T[3][3] = {{t1V[0], t1V[1], t1V[2]},
                              {t2V[0], t2V[1], t2V[2]},
                              {t3V[0], t3V[1], t3V[2]}};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      R_out[r][c] = M_W[r][0] * M_V_T[0][c] + M_W[r][1] * M_V_T[1][c] +
                    M_W[r][2] * M_V_T[2][c];
}

bool quest(const std::vector<std::array<double, 3>> &v_cam,
           const std::vector<std::array<double, 3>> &v_iner,
           double q_out[4],
           int max_iters,
           double tol_lambda) {
  const int M = static_cast<int>(v_cam.size());
  if (M < 2 || static_cast<int>(v_iner.size()) != M) return false;
  const double w = 1.0 / static_cast<double>(M);

  // Attitude profile matrix B = Σ w · b · r^T.
  double B[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int k = 0; k < M; ++k) {
    const auto &b = v_cam[k];
    const auto &r = v_iner[k];
    for (int a = 0; a < 3; ++a)
      for (int bb = 0; bb < 3; ++bb)
        B[a][bb] += w * b[a] * r[bb];
  }

  // S = B + B^T, σ = trace(B), Z = vex(B - B^T).
  double S[3][3];
  for (int a = 0; a < 3; ++a)
    for (int bb = 0; bb < 3; ++bb)
      S[a][bb] = B[a][bb] + B[bb][a];

  const double sigma = B[0][0] + B[1][1] + B[2][2];
  const double Z[3] = {B[2][1] - B[1][2], B[0][2] - B[2][0],
                       B[1][0] - B[0][1]};

  auto det3 = [](const double M3[3][3]) {
    return M3[0][0] * (M3[1][1] * M3[2][2] - M3[1][2] * M3[2][1]) -
           M3[0][1] * (M3[1][0] * M3[2][2] - M3[1][2] * M3[2][0]) +
           M3[0][2] * (M3[1][0] * M3[2][1] - M3[1][1] * M3[2][0]);
  };

  const double trace_adj_S =
      (S[1][1] * S[2][2] - S[1][2] * S[2][1]) +
      (S[0][0] * S[2][2] - S[0][2] * S[2][0]) +
      (S[0][0] * S[1][1] - S[0][1] * S[1][0]);
  const double detS = det3(S);
  const double ZtZ = Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2];
  const double SZ[3] = {S[0][0] * Z[0] + S[0][1] * Z[1] + S[0][2] * Z[2],
                        S[1][0] * Z[0] + S[1][1] * Z[1] + S[1][2] * Z[2],
                        S[2][0] * Z[0] + S[2][1] * Z[1] + S[2][2] * Z[2]};
  const double ZtSZ = Z[0] * SZ[0] + Z[1] * SZ[1] + Z[2] * SZ[2];

  double S2[3][3];
  for (int a = 0; a < 3; ++a)
    for (int bb = 0; bb < 3; ++bb)
      S2[a][bb] = S[a][0] * S[0][bb] + S[a][1] * S[1][bb] + S[a][2] * S[2][bb];
  const double S2Z[3] = {
      S2[0][0] * Z[0] + S2[0][1] * Z[1] + S2[0][2] * Z[2],
      S2[1][0] * Z[0] + S2[1][1] * Z[1] + S2[1][2] * Z[2],
      S2[2][0] * Z[0] + S2[2][1] * Z[1] + S2[2][2] * Z[2]};
  const double ZtS2Z = Z[0] * S2Z[0] + Z[1] * S2Z[1] + Z[2] * S2Z[2];

  // Characteristic-equation coefficients (Shuster):
  //   λ^4 - (a + b) λ^2 - c λ + (a·b + c·σ - d) = 0
  const double a = sigma * sigma - trace_adj_S;
  const double bcoef = sigma * sigma + ZtZ;
  const double c = detS + ZtSZ;
  const double dcoef = ZtS2Z;
  const double coef_lam2 = -(a + bcoef);
  const double coef_lam1 = -c;
  const double coef_lam0 = a * bcoef + c * sigma - dcoef;

  // Newton-Raphson for λ_max. Seed at Σ w_i = 1 (noise-free eigenvalue).
  double lambda = 1.0;
  bool converged = false;
  for (int it = 0; it < max_iters; ++it) {
    const double f = lambda * lambda * lambda * lambda +
                     coef_lam2 * lambda * lambda + coef_lam1 * lambda +
                     coef_lam0;
    const double fp = 4.0 * lambda * lambda * lambda +
                      2.0 * coef_lam2 * lambda + coef_lam1;
    if (std::abs(fp) < 1e-30) break; // singular derivative
    const double delta = f / fp;
    lambda -= delta;
    if (std::abs(delta) < tol_lambda) {
      converged = true;
      break;
    }
  }
  if (!converged || !std::isfinite(lambda)) return false;

  // Shuster's closed-form optimal quaternion:
  //   α = λ^2 - σ^2 + trace(adj(S))
  //   β = λ - σ
  //   γ = (λ + σ)·α - det(S)
  //   X = (α·I + β·S + S^2) · Z
  //   q_xyz = X,  q_w = γ, normalize.
  const double alpha = lambda * lambda - sigma * sigma + trace_adj_S;
  const double beta = lambda - sigma;
  const double gamma = (lambda + sigma) * alpha - detS;
  double Mq[3][3];
  for (int aa = 0; aa < 3; ++aa)
    for (int bb = 0; bb < 3; ++bb) {
      Mq[aa][bb] = beta * S[aa][bb] + S2[aa][bb];
      if (aa == bb) Mq[aa][bb] += alpha;
    }
  const double X[3] = {Mq[0][0] * Z[0] + Mq[0][1] * Z[1] + Mq[0][2] * Z[2],
                       Mq[1][0] * Z[0] + Mq[1][1] * Z[1] + Mq[1][2] * Z[2],
                       Mq[2][0] * Z[0] + Mq[2][1] * Z[1] + Mq[2][2] * Z[2]};
  const double norm_sq =
      gamma * gamma + X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
  if (!std::isfinite(norm_sq) || norm_sq <= 0.0) return false;
  const double inv_n = 1.0 / std::sqrt(norm_sq);

  double qx = X[0] * inv_n;
  double qy = X[1] * inv_n;
  double qz = X[2] * inv_n;
  double qw = gamma * inv_n;

  // Canonicalise q.w >= 0 (matches the public TRIAD output convention on the
  // IdentityRotation test where the expected quaternion is [0, 0, 0, +1]).
  if (qw < 0.0) {
    qx = -qx;
    qy = -qy;
    qz = -qz;
    qw = -qw;
  }
  if (!std::isfinite(qx) || !std::isfinite(qy) || !std::isfinite(qz) ||
      !std::isfinite(qw))
    return false;

  q_out[0] = qx;
  q_out[1] = qy;
  q_out[2] = qz;
  q_out[3] = qw;
  return true;
}

void quat_to_rotation(const double q[4], double R_out[3][3]) {
  const double qx = q[0], qy = q[1], qz = q[2], qw = q[3];
  R_out[0][0] = 1.0 - 2.0 * (qy * qy + qz * qz);
  R_out[0][1] = 2.0 * (qx * qy - qz * qw);
  R_out[0][2] = 2.0 * (qx * qz + qy * qw);
  R_out[1][0] = 2.0 * (qx * qy + qz * qw);
  R_out[1][1] = 1.0 - 2.0 * (qx * qx + qz * qz);
  R_out[1][2] = 2.0 * (qy * qz - qx * qw);
  R_out[2][0] = 2.0 * (qx * qz - qy * qw);
  R_out[2][1] = 2.0 * (qy * qz + qx * qw);
  R_out[2][2] = 1.0 - 2.0 * (qx * qx + qy * qy);
}

bool quest_R(const std::vector<std::array<double, 3>> &v_cam,
             const std::vector<std::array<double, 3>> &v_iner,
             double R_out[3][3],
             int max_iters,
             double tol_lambda) {
  double q[4];
  if (!quest(v_cam, v_iner, q, max_iters, tol_lambda)) return false;
  quat_to_rotation(q, R_out);
  return true;
}

} // namespace wahba
