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
