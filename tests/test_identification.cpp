#include "identification.h"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <set>

// --- 2c. Identification Tests ---
// Uses the real Hipparcos catalog to test the full voting pipeline.

static const char *STAR_FILE = "../data/catalog_stars.bin";
static const char *PAIR_FILE = "../data/catalog_pairs.bin";

class IdentificationTest : public ::testing::Test {
protected:
  void SetUp() override {
    std::ifstream fs(STAR_FILE);
    if (!fs.good()) {
      GTEST_SKIP() << "Catalog files not found — run generate_catalog.py first";
    }
    db = std::make_unique<StarDatabase>(STAR_FILE, PAIR_FILE);
  }

  std::unique_ptr<StarDatabase> db;
};

// Helper: project a catalog star's inertial unit vector through a rotation
// to get a camera-frame centroid in pixel coordinates.
static StarCentroid project_star(const CatalogStar &star,
                                 const PinholeCamera &cam,
                                 const double R[3][3]) {
  // v_cam = R * v_inertial
  double vx = R[0][0] * star.x + R[0][1] * star.y + R[0][2] * star.z;
  double vy = R[1][0] * star.x + R[1][1] * star.y + R[1][2] * star.z;
  double vz = R[2][0] * star.x + R[2][1] * star.y + R[2][2] * star.z;

  StarCentroid c;
  c.x = -1;
  c.y = -1;
  c.intensity = 0;

  if (vz <= 0)
    return c; // behind camera

  c.x = cam.focal_x * (vx / vz) + cam.center_x;
  c.y = cam.focal_y * (vy / vz) + cam.center_y;
  c.intensity = 1000.0;
  return c;
}

TEST_F(IdentificationTest, KnownStarsIdentityRotation) {
  // Use identity rotation (boresight = +Z = north celestial pole).
  // Use a wide 90° FoV so we capture plenty of bright stars in the region.
  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(90.0 * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;

  // Identity rotation matrix
  double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  // Bright stars in or near Ursa Major / Ursa Minor (Dec > 45°, visible with
  // boresight +Z and 90° FoV covering Dec > 45°)
  int test_hips[] = {
      11767, // Polaris (α UMi)     Dec 89°
      82396, // Kochab (β UMi)      Dec 74°
      75097, // Pherkad (γ UMi)     Dec 72°
      72607, // δ UMi               Dec 87°
      54061, // Dubhe (α UMa)       Dec 62°
      53910, // Merak (β UMa)       Dec 56°
      58001, // Phecda (γ UMa)      Dec 54°
      59774, // Megrez (δ UMa)      Dec 57°
      62956, // Alioth (ε UMa)      Dec 56°
      65378, // Mizar (ζ UMa)       Dec 55°
      67301, // Alkaid (η UMa)      Dec 49°
  };

  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;

  for (int hip : test_hips) {
    try {
      CatalogStar s = db->get_star(hip);
      if (s.z <= 0)
        continue; // behind camera

      StarCentroid c = project_star(s, cam, R);
      if (c.x >= 0 && c.x < cam.frame_width && c.y >= 0 &&
          c.y < cam.frame_height) {
        centroids.push_back(c);
        expected_hips.push_back(hip);
      }
    } catch (...) {
      continue;
    }
  }

  ASSERT_GE(centroids.size(), 3u)
      << "Need at least 3 visible stars for identification";

  double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam, *db, cos_tol);

  EXPECT_GE(identified.size(), 2u);

  // Check that identified stars map to the correct HIPs
  std::set<int> identified_hips;
  for (const auto &s : identified) {
    identified_hips.insert(s.catalog_hip_id);
  }

  int correct_count = 0;
  for (int hip : expected_hips) {
    if (identified_hips.count(hip)) {
      correct_count++;
    }
  }
  EXPECT_GE(correct_count, 2)
      << "Expected at least 2 correctly identified stars";
}

TEST_F(IdentificationTest, FalseStarExcluded) {
  // Same setup as above but with an added bogus centroid
  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(90.0 * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;

  double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  int test_hips[] = {
      11767, 82396, 75097, 72607, 54061, 53910,
      58001, 59774, 62956, 65378, 67301,
  };

  std::vector<StarCentroid> centroids;

  for (int hip : test_hips) {
    try {
      CatalogStar s = db->get_star(hip);
      if (s.z <= 0)
        continue;

      StarCentroid c = project_star(s, cam, R);
      if (c.x >= 0 && c.x < cam.frame_width && c.y >= 0 &&
          c.y < cam.frame_height) {
        centroids.push_back(c);
      }
    } catch (...) {
      continue;
    }
  }

  if (centroids.size() < 3) {
    GTEST_SKIP() << "Not enough test stars found";
  }

  // Add a bogus centroid that doesn't correspond to any real star
  StarCentroid bogus;
  bogus.x = 100.0;
  bogus.y = 100.0;
  bogus.intensity = 5000.0;
  int bogus_idx = centroids.size();
  centroids.push_back(bogus);

  double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam, *db, cos_tol);

  // The bogus star should NOT be in the identified list
  for (const auto &s : identified) {
    EXPECT_NE(s.image_idx, bogus_idx)
        << "Bogus centroid should not have been identified";
  }
}
