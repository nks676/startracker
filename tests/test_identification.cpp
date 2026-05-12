#include "identification.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <set>

// --- 2c. Identification Tests ---
// Uses the real Hipparcos catalog to test the full voting pipeline.

static const char *STAR_FILE = "../data/catalog_stars.bin";
static const char *PAIR_FILE = "../data/catalog_pairs.bin";
static const char *PATTERN_FILE_20 = "../data/catalog_patterns_20.bin";

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

// === Phase 3e.3: pattern-hash identification tests ===

// Read the binary catalog at `path` and pull out the pattern record at
// `target_index`. Used by PatternKeyFromObservedReproducesCatalog to anchor
// the test against a known-good (key, hips) tuple without re-running the
// generator. Header layout matches generate_pattern_catalog().
static bool read_pattern_at_index(const std::string &path, int target_index,
                                  uint64_t &out_key,
                                  std::array<int, 4> &out_hips) {
  std::ifstream fs(path, std::ios::binary);
  if (!fs) return false;
  int32_t magic = 0, fov_bin = 0, k_nearest = 0, quant_bits = 0,
          num_patterns = 0;
  fs.read(reinterpret_cast<char *>(&magic), 4);
  fs.read(reinterpret_cast<char *>(&fov_bin), 4);
  fs.read(reinterpret_cast<char *>(&k_nearest), 4);
  fs.read(reinterpret_cast<char *>(&quant_bits), 4);
  fs.read(reinterpret_cast<char *>(&num_patterns), 4);
  if (!fs || magic != 0x50415431) return false;
  if (target_index < 0 || target_index >= num_patterns) return false;
  // Record stride = 8 (key) + 4*4 (hips) = 24 bytes.
  fs.seekg(static_cast<std::streamoff>(target_index) * 24, std::ios::cur);
  fs.read(reinterpret_cast<char *>(&out_key), 8);
  int32_t hips[4] = {0, 0, 0, 0};
  fs.read(reinterpret_cast<char *>(hips), 16);
  if (!fs) return false;
  for (int i = 0; i < 4; ++i) out_hips[i] = hips[i];
  return true;
}

// 3e.3 — load-bearing correctness test. Pick a known catalog 4-tuple, project
// its 4 HIPs through an arbitrary inertial→camera rotation, then run the C++
// canonical_order key on the resulting unit vectors. The returned uint64 key
// must match the catalog entry bit-for-bit. If this fails, the pattern path
// will never match anything.
//
// Wave 1 reported index 435,036 has key=0x365a5261894e7 and
// HIPs=[57211, 57175, 56561, 57669].
TEST_F(IdentificationTest, PatternKeyFromObservedReproducesCatalog) {
  std::ifstream fs(PATTERN_FILE_20);
  if (!fs.good()) {
    GTEST_SKIP() << "Pattern catalog not present — run generate_catalog.py";
  }

  uint64_t expected_key = 0;
  std::array<int, 4> hips{};
  ASSERT_TRUE(read_pattern_at_index(PATTERN_FILE_20, 435036, expected_key, hips))
      << "Could not read pattern at index 435036";
  EXPECT_EQ(expected_key, 0x365a5261894e7ULL);
  EXPECT_EQ(hips[0], 57211);
  EXPECT_EQ(hips[1], 57175);
  EXPECT_EQ(hips[2], 56561);
  EXPECT_EQ(hips[3], 57669);

  // Inertial unit vectors for the 4 catalog stars.
  std::array<std::array<double, 3>, 4> v_inertial{};
  for (int i = 0; i < 4; ++i) {
    CatalogStar s = db->get_star(hips[i]);
    v_inertial[i] = {s.x, s.y, s.z};
  }

  // Arbitrary rotation: 17° about (1, 2, 3)/√14.
  const double angle = 17.0 * M_PI / 180.0;
  const double ax = 1.0 / std::sqrt(14.0);
  const double ay = 2.0 / std::sqrt(14.0);
  const double az = 3.0 / std::sqrt(14.0);
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  // Rodrigues' formula.
  double R[3][3] = {
      {c + ax * ax * (1 - c), ax * ay * (1 - c) - az * s,
       ax * az * (1 - c) + ay * s},
      {ay * ax * (1 - c) + az * s, c + ay * ay * (1 - c),
       ay * az * (1 - c) - ax * s},
      {az * ax * (1 - c) - ay * s, az * ay * (1 - c) + ax * s,
       c + az * az * (1 - c)},
  };

  // Apply R to each inertial vector to get a camera-frame unit vector.
  std::array<std::array<double, 3>, 4> v_cam{};
  for (int i = 0; i < 4; ++i) {
    v_cam[i] = {R[0][0] * v_inertial[i][0] + R[0][1] * v_inertial[i][1] +
                    R[0][2] * v_inertial[i][2],
                R[1][0] * v_inertial[i][0] + R[1][1] * v_inertial[i][1] +
                    R[1][2] * v_inertial[i][2],
                R[2][0] * v_inertial[i][0] + R[2][1] * v_inertial[i][1] +
                    R[2][2] * v_inertial[i][2]};
  }

  // Feed in input-local IDs in some arbitrary order — the canonical-order
  // algorithm is geometry-determined, so the resulting key must be
  // independent of input id (except for lex tie-breaks on bit-equal
  // distances, which don't occur here).
  std::array<int, 4> ids = {100, 200, 300, 400};
  std::array<int, 4> canonical_local{};
  uint64_t key = pattern_key_canonical(v_cam, ids, canonical_local);
  EXPECT_EQ(key, expected_key)
      << "C++ canonical_order_and_key disagrees with Python generator";

  // Bonus: shuffle the inputs and re-derive — same key, different canonical
  // permutation but pointing at the same set.
  std::array<std::array<double, 3>, 4> v_cam_shuf = {v_cam[2], v_cam[0],
                                                      v_cam[3], v_cam[1]};
  std::array<int, 4> ids_shuf = {300, 100, 400, 200};
  std::array<int, 4> can_shuf{};
  uint64_t key_shuf = pattern_key_canonical(v_cam_shuf, ids_shuf, can_shuf);
  EXPECT_EQ(key_shuf, expected_key)
      << "Pattern key should be invariant to input permutation";
}

// 3e.3 — pattern path end-to-end on a synthetic scene. Construct a small
// scene of catalog stars projected through a known rotation, load the
// pattern catalog, run identify_stars, and assert the HIPs come out right.
TEST_F(IdentificationTest, PatternPathIdentifiesSyntheticScene) {
  // Load the FOV-20° pattern catalog. Skip if missing.
  try {
    db->load_pattern_catalog(PATTERN_FILE_20);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "Pattern catalog not available: " << e.what();
  }
  ASSERT_TRUE(db->has_pattern_catalog());

  // 1024x1024, 20° FOV camera.
  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(20.0 * M_PI / 360.0));
  cam.focal_y = cam.focal_x;

  // Boresight on Polaris (HIP 11767, near +Z but offset enough to expose
  // pyramid-vs-pattern path differences). Use HIP 32349 (Sirius) at Dec ≈
  // -16°, RA ≈ 6h45m. To get a dense field we point at a busy region: pick
  // a star and use it as boresight.
  CatalogStar bs = db->get_star(11767); // Polaris
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
  for (int k = 0; k < 3; ++k) c1[k] /= n1;
  double c2[3] = {c3[1] * c1[2] - c3[2] * c1[1],
                  c3[2] * c1[0] - c3[0] * c1[2],
                  c3[0] * c1[1] - c3[1] * c1[0]};
  double R[3][3] = {{c1[0], c1[1], c1[2]},
                    {c2[0], c2[1], c2[2]},
                    {c3[0], c3[1], c3[2]}};

  // Stars within ~5° of Polaris (HIP 11767). Precomputed from the catalog
  // so the test doesn't have to scan it at runtime.
  std::vector<int> candidates = {
      11767, 7283,   84535, 5928,   115746, 19454, 118285,
      37391, 59767,  5372,  32948,  115550, 71725, 109694,
      10800, 25911,  85699, 109693, 89465,  21386,
  };
  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  for (int hip : candidates) {
    CatalogStar s;
    try { s = db->get_star(hip); } catch (...) { continue; }
    double vz = R[2][0] * s.x + R[2][1] * s.y + R[2][2] * s.z;
    if (vz <= 0.5) continue; // require z > 0.5 (FOV margin)
    double vx = R[0][0] * s.x + R[0][1] * s.y + R[0][2] * s.z;
    double vy = R[1][0] * s.x + R[1][1] * s.y + R[1][2] * s.z;
    StarCentroid c;
    c.x = cam.focal_x * (vx / vz) + cam.center_x;
    c.y = cam.focal_y * (vy / vz) + cam.center_y;
    c.intensity = 1000.0;
    c.peak = 250.0 - centroids.size(); // descending peak; first = brightest
    if (c.x >= 0 && c.x < cam.frame_width && c.y >= 0 &&
        c.y < cam.frame_height) {
      centroids.push_back(c);
      expected_hips.push_back(hip);
    }
  }
  ASSERT_GE(centroids.size(), 6u)
      << "Need ≥6 stars in FOV for the pattern path";

  auto identified = identify_stars(centroids, cam, *db, 1e-5);

  // The pattern path returns 4 + 1 verified + possibly more from expansion.
  ASSERT_GE(identified.size(), 4u);
  std::set<int> expected(expected_hips.begin(), expected_hips.end());
  for (const auto &s : identified) {
    EXPECT_TRUE(expected.count(s.catalog_hip_id))
        << "Pattern path produced an unexpected HIP " << s.catalog_hip_id;
  }
}

// 3e.3 — fallback path: an instance of StarDatabase that has NOT had a
// pattern catalog loaded must still produce correct results via the pyramid.
// Smoke-tests that the pattern path doesn't short-circuit on `false` from
// has_pattern_catalog().
TEST_F(IdentificationTest, PatternPathFallsBackOnMissingCatalog) {
  // Fresh DB with no pattern catalog loaded.
  StarDatabase fresh(STAR_FILE, PAIR_FILE);
  ASSERT_FALSE(fresh.has_pattern_catalog());

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
  gather_visible_stars(candidates, fresh, cam, R, centroids, expected_hips);
  ASSERT_GE(centroids.size(), 6u);

  auto identified = identify_stars(centroids, cam, fresh, 1e-5);
  // Pyramid path on this dense Ursa Major scene must still work.
  EXPECT_GE(identified.size(), 6u);
  std::set<int> truth(expected_hips.begin(), expected_hips.end());
  int correct = 0;
  for (const auto &s : identified)
    if (truth.count(s.catalog_hip_id)) ++correct;
  EXPECT_GE(correct, 6) << "Pyramid fallback should identify the bright stars";
}
