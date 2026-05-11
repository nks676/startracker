#include "estimation.h"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <sstream>

// --- 2d. Estimation Tests ---
// Uses a small synthetic catalog (written to temp files) so these tests
// are fully self-contained — no real Hipparcos data needed.

// Helper: write a minimal binary catalog with a few stars at known positions
static std::string write_temp_stars(
    const std::vector<std::tuple<int, double, double, double>> &stars) {
  std::string path = "/tmp/test_stars.bin";
  std::ofstream f(path, std::ios::binary);
  int n = stars.size();
  f.write(reinterpret_cast<char *>(&n), sizeof(int));
  for (auto &[id, x, y, z] : stars) {
    int hip = id;
    double vx = x, vy = y, vz = z;
    f.write(reinterpret_cast<char *>(&hip), sizeof(int));
    f.write(reinterpret_cast<char *>(&vx), sizeof(double));
    f.write(reinterpret_cast<char *>(&vy), sizeof(double));
    f.write(reinterpret_cast<char *>(&vz), sizeof(double));
  }
  return path;
}

static std::string write_empty_pairs() {
  std::string path = "/tmp/test_pairs.bin";
  std::ofstream f(path, std::ios::binary);
  int n = 0;
  f.write(reinterpret_cast<char *>(&n), sizeof(int));
  return path;
}

// Helper: build IdentifiedStar from camera-frame vector and catalog HIP
static IdentifiedStar make_star(int idx, int hip, double cx, double cy,
                                double cz) {
  IdentifiedStar s;
  s.image_idx = idx;
  s.catalog_hip_id = hip;
  s.v_cam[0] = cx;
  s.v_cam[1] = cy;
  s.v_cam[2] = cz;
  return s;
}

// Helper: compute angular error between two quaternions (in degrees)
static double quat_angle_error(const Quaternion &a, const Quaternion &b) {
  // dot product of quaternions
  double d = std::abs(a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w);
  if (d > 1.0)
    d = 1.0;
  return 2.0 * std::acos(d) * 180.0 / M_PI;
}

TEST(EstimationTest, IdentityRotation) {
  // Two stars with camera vectors = inertial vectors (no rotation → R = I)
  // Star 1 at [1, 0, 0], Star 2 at [0, 1, 0]
  auto star_path = write_temp_stars({
      {1, 1.0, 0.0, 0.0},
      {2, 0.0, 1.0, 0.0},
      {3, 0.0, 0.0, 1.0},
  });
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  std::vector<IdentifiedStar> stars = {
      make_star(0, 1, 1.0, 0.0, 0.0), // cam = inertial
      make_star(1, 2, 0.0, 1.0, 0.0),
      make_star(2, 3, 0.0, 0.0, 1.0),
  };

  Quaternion q = estimate_attitude(stars, db);

  // Identity quaternion: [0, 0, 0, 1]
  EXPECT_NEAR(std::abs(q.w), 1.0, 0.001)
      << "q = [" << q.x << ", " << q.y << ", " << q.z << ", " << q.w << "]";
  double err = quat_angle_error(q, {0, 0, 0, 1});
  EXPECT_LT(err, 0.1) << "Angular error: " << err << " degrees";
}

TEST(EstimationTest, Known90DegZRotation) {
  // Rotation of 90° around Z axis:
  //   R = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
  //   quaternion = [0, 0, sin(π/4), cos(π/4)] = [0, 0, 0.7071, 0.7071]
  //
  // If v_inertial = [1,0,0] then v_cam = R * v_inertial = [0,1,0]
  // If v_inertial = [0,1,0] then v_cam = R * v_inertial = [-1,0,0]
  // If v_inertial = [0,0,1] then v_cam = R * v_inertial = [0,0,1]

  auto star_path = write_temp_stars({
      {1, 1.0, 0.0, 0.0},
      {2, 0.0, 1.0, 0.0},
      {3, 0.0, 0.0, 1.0},
  });
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  std::vector<IdentifiedStar> stars = {
      make_star(0, 1, 0.0, 1.0, 0.0),  // R * [1,0,0] = [0,1,0]
      make_star(1, 2, -1.0, 0.0, 0.0), // R * [0,1,0] = [-1,0,0]
      make_star(2, 3, 0.0, 0.0, 1.0),  // R * [0,0,1] = [0,0,1]
  };

  Quaternion q = estimate_attitude(stars, db);

  // Expected: [0, 0, sin(45°), cos(45°)]
  double s45 = std::sin(M_PI / 4.0);
  double c45 = std::cos(M_PI / 4.0);
  Quaternion expected = {0, 0, s45, c45};

  double err = quat_angle_error(q, expected);
  EXPECT_LT(err, 0.1) << "Angular error: " << err << " degrees";
}

TEST(EstimationTest, RandomRotationConsistency) {
  // Generate a random rotation, apply it, run TRIAD, verify low error.
  std::mt19937 rng(12345);
  std::normal_distribution<double> dist(0, 1);

  // Random unit quaternion
  double qx = dist(rng), qy = dist(rng), qz = dist(rng), qw = dist(rng);
  double norm = std::sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
  qx /= norm;
  qy /= norm;
  qz /= norm;
  qw /= norm;
  if (qw < 0) {
    qx = -qx;
    qy = -qy;
    qz = -qz;
    qw = -qw;
  } // canonical

  // Quaternion to rotation matrix
  double R[3][3] = {
      {1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qz * qw),
       2 * (qx * qz + qy * qw)},
      {2 * (qx * qy + qz * qw), 1 - 2 * (qx * qx + qz * qz),
       2 * (qy * qz - qx * qw)},
      {2 * (qx * qz - qy * qw), 2 * (qy * qz + qx * qw),
       1 - 2 * (qx * qx + qy * qy)},
  };

  // Three orthogonal inertial vectors
  double v_inertial[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };

  // Rotate to get camera vectors
  double v_cam[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      v_cam[i][j] = R[j][0] * v_inertial[i][0] + R[j][1] * v_inertial[i][1] +
                    R[j][2] * v_inertial[i][2];
    }
  }

  auto star_path = write_temp_stars({
      {1, 1.0, 0.0, 0.0},
      {2, 0.0, 1.0, 0.0},
      {3, 0.0, 0.0, 1.0},
  });
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  std::vector<IdentifiedStar> stars = {
      make_star(0, 1, v_cam[0][0], v_cam[0][1], v_cam[0][2]),
      make_star(1, 2, v_cam[1][0], v_cam[1][1], v_cam[1][2]),
      make_star(2, 3, v_cam[2][0], v_cam[2][1], v_cam[2][2]),
  };

  Quaternion q = estimate_attitude(stars, db);
  Quaternion expected = {qx, qy, qz, qw};
  double err = quat_angle_error(q, expected);

  EXPECT_LT(err, 0.01) << "Angular error: " << err << "° (expected q=[" << qx
                       << "," << qy << "," << qz << "," << qw << "])";
}

// --- QUEST-specific tests (Phase 3b.2) ---

// Apply quaternion q (as a rotation matrix) to a 3-vector.
static void rotate(const Quaternion &q, const double v[3], double out[3]) {
  double xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
  double xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
  double wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
  double R[3][3] = {
      {1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy)},
      {2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx)},
      {2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)},
  };
  for (int i = 0; i < 3; ++i) {
    out[i] = R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2];
  }
}

static void normalize3(double v[3]) {
  double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= n;
  v[1] /= n;
  v[2] /= n;
}

// QuestBeatsTriadWithNoise: 5 stars with ~30 arcsec jitter; QUEST should
// outperform a 2-star TRIAD solution on average across noise realisations.
// We reuse the public entry point for QUEST; to get a comparable TRIAD-only
// result we feed just 2 stars (which triggers the TRIAD fallback path).
TEST(EstimationTest, QuestBeatsTriadWithNoise) {
  // Five inertial directions spread across the sky (well-conditioned set)
  std::vector<std::tuple<int, double, double, double>> star_defs = {
      {1, 1.0, 0.0, 0.0},
      {2, 0.0, 1.0, 0.0},
      {3, 0.0, 0.0, 1.0},
      {4, 0.7071, 0.7071, 0.0},
      {5, 0.5774, 0.5774, 0.5774},
  };
  for (auto &t : star_defs) {
    double vx = std::get<1>(t), vy = std::get<2>(t), vz = std::get<3>(t);
    double n = std::sqrt(vx * vx + vy * vy + vz * vz);
    std::get<1>(t) = vx / n;
    std::get<2>(t) = vy / n;
    std::get<3>(t) = vz / n;
  }

  auto star_path = write_temp_stars(star_defs);
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  // Use a small, fixed rotation
  Quaternion q_truth = {0.05, 0.08, 0.06, 0.99};
  double qn = std::sqrt(q_truth.x * q_truth.x + q_truth.y * q_truth.y +
                        q_truth.z * q_truth.z + q_truth.w * q_truth.w);
  q_truth.x /= qn;
  q_truth.y /= qn;
  q_truth.z /= qn;
  q_truth.w /= qn;

  std::mt19937 rng(98765);
  // ~30 arcsec ≈ 1.45e-4 rad jitter on each component
  std::normal_distribution<double> noise(0.0, 1.0e-4);

  // Average error across multiple noise realisations (single-trial comparison
  // is noisy; with M=20 trials QUEST's variance reduction over TRIAD is
  // statistically detectable).
  const int kTrials = 20;
  double sum_err_quest = 0.0, sum_err_triad = 0.0;
  for (int trial = 0; trial < kTrials; ++trial) {
    std::vector<IdentifiedStar> stars5;
    for (size_t i = 0; i < star_defs.size(); ++i) {
      auto &t = star_defs[i];
      double v_in[3] = {std::get<1>(t), std::get<2>(t), std::get<3>(t)};
      double v_cam[3];
      rotate(q_truth, v_in, v_cam);
      v_cam[0] += noise(rng);
      v_cam[1] += noise(rng);
      v_cam[2] += noise(rng);
      normalize3(v_cam);
      stars5.push_back(
          make_star((int)i, std::get<0>(t), v_cam[0], v_cam[1], v_cam[2]));
    }
    std::vector<IdentifiedStar> stars2 = {stars5[0], stars5[1]};

    Quaternion q_quest = estimate_attitude(stars5, db);
    Quaternion q_triad = estimate_attitude(stars2, db);

    sum_err_quest += quat_angle_error(q_quest, q_truth);
    sum_err_triad += quat_angle_error(q_triad, q_truth);
  }

  double avg_err_quest = sum_err_quest / kTrials;
  double avg_err_triad = sum_err_triad / kTrials;

  EXPECT_LT(avg_err_quest, avg_err_triad)
      << "avg QUEST err=" << avg_err_quest << "°, avg TRIAD err=" << avg_err_triad
      << "° over " << kTrials << " trials";
  EXPECT_LT(avg_err_quest, 0.05)
      << "avg QUEST err=" << avg_err_quest
      << "° should be small with 1e-4 rad noise";
}

// QuestFallbackOnTwoStars: 2-star input must still return a finite, sensible
// quaternion (TRIAD fallback path).
TEST(EstimationTest, QuestFallbackOnTwoStars) {
  auto star_path = write_temp_stars({
      {1, 1.0, 0.0, 0.0},
      {2, 0.0, 1.0, 0.0},
  });
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  // Identity rotation: cam == inertial.
  std::vector<IdentifiedStar> stars = {
      make_star(0, 1, 1.0, 0.0, 0.0),
      make_star(1, 2, 0.0, 1.0, 0.0),
  };

  Quaternion q = estimate_attitude(stars, db);

  // Finite components
  EXPECT_TRUE(std::isfinite(q.x));
  EXPECT_TRUE(std::isfinite(q.y));
  EXPECT_TRUE(std::isfinite(q.z));
  EXPECT_TRUE(std::isfinite(q.w));

  // Within a reasonable tolerance of identity
  double err = quat_angle_error(q, {0, 0, 0, 1});
  EXPECT_LT(err, 0.1) << "TRIAD fallback err = " << err << "°";
}

// QuestStableNearCollinear: 3 stars clustered within ~5° of a great circle.
// QUEST should not produce NaN; it may converge OR may trip the TRIAD
// fallback, but the returned quaternion must be finite and reasonable.
TEST(EstimationTest, QuestStableNearCollinear) {
  // Three stars whose inertial vectors all lie within ~5° of the x-axis great
  // circle (i.e., y ≈ 0, x and z varying).
  std::vector<std::tuple<int, double, double, double>> defs;
  // Star A: pure +x
  defs.push_back({1, 1.0, 0.0, 0.0});
  // Star B: rotated 10° toward +z (still near x-axis)
  double ang_b = 10.0 * M_PI / 180.0;
  defs.push_back({2, std::cos(ang_b), 0.0, std::sin(ang_b)});
  // Star C: rotated 20° toward +z, with a tiny y offset (~3°) so it isn't
  // exactly on the line A–B but still within 5° of the great circle in y.
  double ang_c = 20.0 * M_PI / 180.0;
  double tilt = 3.0 * M_PI / 180.0;
  double cx = std::cos(ang_c) * std::cos(tilt);
  double cy = std::sin(tilt);
  double cz = std::sin(ang_c) * std::cos(tilt);
  defs.push_back({3, cx, cy, cz});

  auto star_path = write_temp_stars(defs);
  auto pair_path = write_empty_pairs();
  StarDatabase db(star_path, pair_path);

  // Small rotation truth.
  Quaternion q_truth = {0.02, 0.0, 0.03, 1.0};
  double n = std::sqrt(q_truth.x * q_truth.x + q_truth.y * q_truth.y +
                       q_truth.z * q_truth.z + q_truth.w * q_truth.w);
  q_truth.x /= n;
  q_truth.y /= n;
  q_truth.z /= n;
  q_truth.w /= n;

  std::vector<IdentifiedStar> stars;
  for (size_t i = 0; i < defs.size(); ++i) {
    auto &t = defs[i];
    double v_in[3] = {std::get<1>(t), std::get<2>(t), std::get<3>(t)};
    double v_cam[3];
    rotate(q_truth, v_in, v_cam);
    normalize3(v_cam);
    stars.push_back(
        make_star((int)i, std::get<0>(t), v_cam[0], v_cam[1], v_cam[2]));
  }

  Quaternion q = estimate_attitude(stars, db);

  // No NaN/Inf
  EXPECT_TRUE(std::isfinite(q.x)) << "q.x = " << q.x;
  EXPECT_TRUE(std::isfinite(q.y)) << "q.y = " << q.y;
  EXPECT_TRUE(std::isfinite(q.z)) << "q.z = " << q.z;
  EXPECT_TRUE(std::isfinite(q.w)) << "q.w = " << q.w;

  // Whether QUEST or TRIAD fallback ran, the answer should be reasonable
  // (well below 1°)
  double err = quat_angle_error(q, q_truth);
  EXPECT_LT(err, 1.0) << "Angular error: " << err << "°";
}
