#include "identification.h"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <set>

// --- 3a.1: Brown-Conrady camera distortion model tests ---

namespace {

CameraModel make_distorted_cam(double fov_deg, double k1, double k2, double p1,
                               double p2, double k3) {
  CameraModel cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = cam.frame_width / 2.0;
  cam.center_y = cam.frame_height / 2.0;
  cam.focal_x = cam.frame_width / (2.0 * std::tan(fov_deg * M_PI / 360.0));
  cam.focal_y = cam.focal_x;
  cam.k1 = k1;
  cam.k2 = k2;
  cam.k3 = k3;
  cam.p1 = p1;
  cam.p2 = p2;
  return cam;
}

} // namespace

// 1. Round-trip identity: undistort(px, py) followed by project() should
// reproduce the original pixel within tight numerical tolerance for a
// moderately distorted camera.
TEST(CameraModelTest, UndistortProjectRoundTrip) {
  CameraModel cam =
      make_distorted_cam(/*fov_deg=*/20.0, /*k1=*/-0.3, /*k2=*/0.1,
                         /*p1=*/0.0, /*p2=*/0.0, /*k3=*/0.0);

  std::mt19937 rng(0xC0FFEE);
  std::uniform_real_distribution<double> ux(50.0, cam.frame_width - 50.0);
  std::uniform_real_distribution<double> uy(50.0, cam.frame_height - 50.0);

  for (int i = 0; i < 20; ++i) {
    double px = ux(rng);
    double py = uy(rng);

    double v[3];
    undistort_to_unit_vector(cam, px, py, v);

    double px2, py2;
    project(cam, v, px2, py2);

    EXPECT_NEAR(px, px2, 1e-6) << "iteration " << i;
    EXPECT_NEAR(py, py2, 1e-6) << "iteration " << i;
  }
}

// 2. Zero-distortion equivalence: undistort_to_unit_vector with all
// coefficients = 0 must match the textbook pinhole inverse exactly (to
// floating-point precision).
TEST(CameraModelTest, ZeroDistortionEqualsPinhole) {
  CameraModel cam = make_distorted_cam(20.0, 0, 0, 0, 0, 0);

  std::mt19937 rng(0xBADF00D);
  std::uniform_real_distribution<double> ux(0.0, cam.frame_width);
  std::uniform_real_distribution<double> uy(0.0, cam.frame_height);

  for (int i = 0; i < 20; ++i) {
    double px = ux(rng);
    double py = uy(rng);

    // Legacy inline math from identification.cpp pre-3a.1.
    double vx = (px - cam.center_x) / cam.focal_x;
    double vy = (py - cam.center_y) / cam.focal_y;
    double vz = 1.0;
    double n = std::sqrt(vx * vx + vy * vy + vz * vz);
    double ref[3] = {vx / n, vy / n, vz / n};

    double v[3];
    undistort_to_unit_vector(cam, px, py, v);

    EXPECT_NEAR(v[0], ref[0], 1e-12);
    EXPECT_NEAR(v[1], ref[1], 1e-12);
    EXPECT_NEAR(v[2], ref[2], 1e-12);
  }
}

// 3. End-to-end distorted identification: project real catalog stars through
// a known rotation + distortion, then verify identify_stars recovers their
// HIPs when given the matching CameraModel. This is the killer test — it
// proves the inverse model is correct enough to feed the voting pipeline.
class DistortedIdentificationTest : public ::testing::Test {
protected:
  static constexpr const char *STAR_FILE = "../data/catalog_stars.bin";
  static constexpr const char *PAIR_FILE = "../data/catalog_pairs.bin";

  void SetUp() override {
    std::ifstream fs(STAR_FILE);
    if (!fs.good()) {
      GTEST_SKIP() << "Catalog files not found — run generate_catalog.py first";
    }
    db = std::make_unique<StarDatabase>(STAR_FILE, PAIR_FILE);
  }

  std::unique_ptr<StarDatabase> db;
};

TEST_F(DistortedIdentificationTest, BarrelDistortedStarsIdentified) {
  // Wide-FoV camera with barrel distortion typical of cheap CS-mount lenses.
  CameraModel cam =
      make_distorted_cam(/*fov_deg=*/90.0, /*k1=*/-0.3, /*k2=*/0.1,
                         /*p1=*/0.0, /*p2=*/0.0, /*k3=*/0.0);

  // Same well-known Ursa Major / Minor stars used by IdentificationTest.
  int test_hips[] = {
      11767, 82396, 75097, 72607, 54061, 53910,
      58001, 59774, 62956, 65378, 67301,
  };

  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;

  for (int hip : test_hips) {
    try {
      CatalogStar s = db->get_star(hip);
      double v[3] = {s.x, s.y, s.z};
      if (v[2] <= 0)
        continue;
      double px, py;
      project(cam, v, px, py);
      if (px < 0 || px >= cam.frame_width || py < 0 || py >= cam.frame_height)
        continue;
      StarCentroid c;
      c.x = px;
      c.y = py;
      c.intensity = 1000.0;
      centroids.push_back(c);
      expected_hips.push_back(hip);
    } catch (...) {
      continue;
    }
  }

  ASSERT_GE(centroids.size(), 4u);

  auto identified = identify_stars(centroids, cam, *db, 1e-5);
  ASSERT_GE(identified.size(), 3u);

  std::set<int> identified_hips;
  for (const auto &s : identified)
    identified_hips.insert(s.catalog_hip_id);

  int correct = 0;
  for (int hip : expected_hips)
    if (identified_hips.count(hip))
      correct++;
  EXPECT_GE(correct, 3)
      << "Expected at least 3 correctly identified stars under k1=-0.3 barrel "
         "distortion";
}
