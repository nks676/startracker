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

// Helper: gather centroids and HIPs visible in the camera FOV given a fixed
// boresight and a starting list of candidate HIPs. Returns pairs of
// (centroid, hip) so tests can build "ground-truth" mappings.
static void gather_visible_stars(
    const std::vector<int> &candidate_hips, const StarDatabase &db,
    const PinholeCamera &cam, const double R[3][3],
    std::vector<StarCentroid> &centroids_out,
    std::vector<int> &hips_out) {
  for (int hip : candidate_hips) {
    try {
      CatalogStar s = db.get_star(hip);
      double vz = R[2][0] * s.x + R[2][1] * s.y + R[2][2] * s.z;
      if (vz <= 0)
        continue;
      StarCentroid c = project_star(s, cam, R);
      if (c.x >= 0 && c.x < cam.frame_width && c.y >= 0 &&
          c.y < cam.frame_height) {
        centroids_out.push_back(c);
        hips_out.push_back(hip);
      }
    } catch (...) {
      continue;
    }
  }
}

// 3b.0b: simulate an FOV calibration error by rendering the scene at the
// "true" focal length, then passing the identification routine a camera
// whose focal length is scaled by 1.004 (a 0.4% drift comparable to the
// real-image alt60 fixture). At default cos_tol=1e-5 the pyramid alone
// can't bridge the gap, but the coarse-refine-reidentify wrapper estimates
// the scale factor from the coarse pass and corrects the camera before the
// tight second pass.
TEST_F(IdentificationTest, RefinedFovRecoversAlt60Scenario) {
  // Mirror the alt60 ESA tetra3 fixture: 1024x768 frame, 11.5° horizontal
  // FOV, and a dense star field (matching the real image's ~25 detectable
  // stars). The "without refinement" failure mode is a 0.4% FOV calibration
  // drift between assumed and true focal length, which puts catalog-pair
  // matches just outside the default cos_tolerance of 1e-5 for pairs around
  // 5-11° apart. The coarse-refine-reidentify wrapper has to recover from
  // this without changing the public cos_tol input.
  PinholeCamera cam_truth;
  cam_truth.frame_width = 1024;
  cam_truth.frame_height = 768;
  cam_truth.center_x = 512.0;
  cam_truth.center_y = 384.0;
  const double fov_deg = 11.5;
  cam_truth.focal_x =
      1024.0 / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  cam_truth.focal_y = cam_truth.focal_x;

  // Boresight on HIP 3334 — a relatively dense catalog region (~70+ stars
  // within 5.5°), which produces a realistic 20-25 star centroid list at
  // 11.5° FOV. Density matters: the coarse pyramid pass needs many
  // small-angle pairs (whose drift stays well under 1e-4) to lock onto the
  // correct constellation.
  CatalogStar bs = db->get_star(3334);
  double c3[3] = {bs.x, bs.y, bs.z};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0;
    up[1] = 0.0;
    up[2] = 0.0;
  }
  double c1[3] = {up[1] * c3[2] - up[2] * c3[1],
                  up[2] * c3[0] - up[0] * c3[2],
                  up[0] * c3[1] - up[1] * c3[0]};
  double n1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
  for (int k = 0; k < 3; ++k)
    c1[k] /= n1;
  double c2[3] = {c3[1] * c1[2] - c3[2] * c1[1],
                  c3[2] * c1[0] - c3[0] * c1[2],
                  c3[0] * c1[1] - c3[1] * c1[0]};
  double R[3][3] = {{c1[0], c1[1], c1[2]},
                    {c2[0], c2[1], c2[2]},
                    {c3[0], c3[1], c3[2]}};

  // 25 catalog HIPs within ~6° of HIP 3334, listed in ascending order of
  // angular distance. Pre-computed offline from `data/catalog_stars.bin`
  // so the test doesn't need to scan the catalog at runtime.
  std::vector<int> candidates = {
      3334, 3649, 3058, 3821, 4292, 3030, 2876, 4422, 4383,
      2377, 3584, 3179, 4440, 4427, 4811, 4151, 4961, 2074,
      2101, 1892, 3988, 4475, 5232, 2440, 5361,
  };

  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam_truth, R, centroids, expected_hips);

  // Need a dense scene (~20+ stars) for the coarse pass to disambiguate.
  ASSERT_GE(centroids.size(), 15u)
      << "Not enough stars in synthetic FOV for refinement test";

  // Now construct the "assumed" camera: focal_x and focal_y scaled by 1.004
  // (so identification sees a 0.4% FOV calibration drift). The undistort
  // step in identify_stars will produce v_cam values whose pair angles are
  // slightly too small to match the catalog within cos_tol=1e-5.
  PinholeCamera cam_assumed = cam_truth;
  cam_assumed.focal_x *= 1.004;
  cam_assumed.focal_y *= 1.004;

  const double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam_assumed, *db, cos_tol);

  // With refinement the coarse pass (cos_tol_coarse = 1e-4) catches the
  // small-angle pairs, estimates s ~ 1/1.004 ≈ 0.996, rescales focal back
  // to truth, and the tight second pass identifies essentially all of them.
  // Without refinement: identify_stars at default cos_tol would return {}.
  ASSERT_GE(identified.size(), expected_hips.size() / 2)
      << "Coarse-refine-reidentify failed to recover at default cos_tol; "
      << "got " << identified.size() << " of " << expected_hips.size();

  std::set<int> identified_hips;
  for (const auto &s : identified)
    identified_hips.insert(s.catalog_hip_id);

  int correct = 0;
  for (int hip : expected_hips)
    if (identified_hips.count(hip))
      ++correct;
  EXPECT_GE(correct, (int)expected_hips.size() / 2)
      << "FOV refinement should produce mostly-correct HIPs; correct="
      << correct;
}

// 3b.3: cross-verification rejects a single deliberately-mislabelled star
// without removing the rest. We can't reach into identify_stars to inject a
// wrong HIP directly, so we exploit cross-verification's input format: a
// scene where one image star has been moved to a wrong sky location that
// the pyramid happens to label as a nearby (but incorrect) HIP. The cleaner
// way is to set up a synthetic "clean" identification, then re-run the
// pipeline with one centroid replaced by a centroid that corresponds to a
// different real star. We verify that the swap-victim is the only star
// missing from the result.
TEST_F(IdentificationTest, CrossVerificationRejectsMislabel) {
  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(90.0 * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;
  double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  std::vector<int> candidates = {
      11767, 82396, 75097, 72607, 54061, 53910,
      58001, 59774, 62956, 65378, 67301,
  };
  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam, R, centroids, expected_hips);
  ASSERT_GE(centroids.size(), 6u);

  // Move the last centroid to a location that does NOT correspond to any
  // catalog star projection. Coordinates well off the actual sky pattern
  // mean the pyramid won't be able to find a consistent HIP for it, so it
  // either won't appear in `identified` or, if it does (mismatched HIP),
  // cross-verification will drop it. Either way: the OTHER stars should
  // still be identified.
  int victim_idx = static_cast<int>(centroids.size()) - 1;
  centroids[victim_idx].x = 50.0;
  centroids[victim_idx].y = 950.0;

  const double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam, *db, cos_tol);

  std::set<int> identified_hips;
  bool victim_image_present = false;
  for (const auto &s : identified) {
    identified_hips.insert(s.catalog_hip_id);
    if (s.image_idx == victim_idx)
      victim_image_present = true;
  }
  // The victim's image index should not appear in the result — either the
  // pyramid didn't find a consistent HIP for it, or cross-verification
  // rejected it.
  EXPECT_FALSE(victim_image_present)
      << "Victim centroid should have been removed";
  // And the rest of the (clean) stars should still be there.
  int kept = 0;
  for (int i = 0; i < (int)expected_hips.size(); ++i) {
    if (i == victim_idx)
      continue;
    if (identified_hips.count(expected_hips[i]))
      ++kept;
  }
  EXPECT_GE(kept, (int)expected_hips.size() - 2)
      << "Cross-verification dropped too many clean stars";
}

// 3b.3: a clean scene should not trigger cross-verification rejection.
TEST_F(IdentificationTest, CrossVerificationPreservesGoodAssignments) {
  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(90.0 * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;
  double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  std::vector<int> candidates = {
      11767, 82396, 75097, 72607, 54061, 53910,
      58001, 59774, 62956, 65378, 67301,
  };
  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam, R, centroids, expected_hips);
  ASSERT_GE(centroids.size(), 8u)
      << "Need at least 8 clean stars for this scenario";

  const double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam, *db, cos_tol);
  EXPECT_GE(identified.size(), 8u)
      << "Cross-verification should not drop any clean stars";

  // Every identified HIP must be one of the ground-truth HIPs.
  std::set<int> truth(expected_hips.begin(), expected_hips.end());
  for (const auto &s : identified)
    EXPECT_TRUE(truth.count(s.catalog_hip_id))
        << "Identified an unexpected HIP " << s.catalog_hip_id;
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
