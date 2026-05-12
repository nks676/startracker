#include "estimation.h"
#include "identification.h"
#include <algorithm>
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
static const char *PATTERN_FILE_10 = "../data/catalog_patterns_10.bin";
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
// Anchored on the FIRST pattern record of the file (index 0) rather than a
// hand-picked offset, so the test stays valid through catalog regenerations
// (the K_NEAREST bump from 8 → 12 changed every index by ≈4×).
TEST_F(IdentificationTest, PatternKeyFromObservedReproducesCatalog) {
  std::ifstream fs(PATTERN_FILE_20);
  if (!fs.good()) {
    GTEST_SKIP() << "Pattern catalog not present — run generate_catalog.py";
  }

  uint64_t expected_key = 0;
  std::array<int, 4> hips{};
  ASSERT_TRUE(read_pattern_at_index(PATTERN_FILE_20, 0, expected_key, hips))
      << "Could not read pattern at index 0";
  // The first record's specific (key, HIPs) is determined by the catalog
  // generator and changes when K_NEAREST or the input star set changes; we
  // intentionally don't pin to particular values here.

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

// === Phase 3e.5: noise-robust permutation probing, sweep, post-verify accuracy ===
//
// Build a 4-star scene at 10° FOV from a known catalog pattern, project under
// a known rotation, and recover the 4 camera-frame unit vectors. Shared by
// the 3e.5 tests below.
static void build_known_4star_scene(
    const StarDatabase &db,
    const PinholeCamera &cam, double R[3][3],
    std::array<int, 4> &hips,
    std::array<std::array<double, 3>, 4> &v_cam_out) {
  // Re-use index 435036 from the FOV-20 catalog: a verified-good 4-tuple
  // whose canonical key was confirmed in PatternKeyFromObservedReproducesCatalog.
  uint64_t expected_key = 0;
  ASSERT_TRUE(read_pattern_at_index(PATTERN_FILE_20, 435036, expected_key, hips))
      << "Could not read anchor pattern from catalog";

  // Build a frame whose +Z axis is aligned with the centroid of the 4 HIPs.
  double cx = 0.0, cy = 0.0, cz = 0.0;
  for (int hip : hips) {
    CatalogStar s = db.get_star(hip);
    cx += s.x; cy += s.y; cz += s.z;
  }
  double n = std::sqrt(cx * cx + cy * cy + cz * cz);
  ASSERT_GT(n, 1e-6);
  double c3[3] = {cx / n, cy / n, cz / n};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
  }
  double c1[3] = {up[1] * c3[2] - up[2] * c3[1],
                  up[2] * c3[0] - up[0] * c3[2],
                  up[0] * c3[1] - up[1] * c3[0]};
  double n1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
  for (int k = 0; k < 3; ++k) c1[k] /= n1;
  double c2[3] = {c3[1] * c1[2] - c3[2] * c1[1],
                  c3[2] * c1[0] - c3[0] * c1[2],
                  c3[0] * c1[1] - c3[1] * c1[0]};
  R[0][0] = c1[0]; R[0][1] = c1[1]; R[0][2] = c1[2];
  R[1][0] = c2[0]; R[1][1] = c2[1]; R[1][2] = c2[2];
  R[2][0] = c3[0]; R[2][1] = c3[1]; R[2][2] = c3[2];

  (void)cam; // camera unused at this granularity — we only need v_cam
  for (int i = 0; i < 4; ++i) {
    CatalogStar s = db.get_star(hips[i]);
    v_cam_out[i] = {
        R[0][0] * s.x + R[0][1] * s.y + R[0][2] * s.z,
        R[1][0] * s.x + R[1][1] * s.y + R[1][2] * s.z,
        R[2][0] * s.x + R[2][1] * s.y + R[2][2] * s.z,
    };
  }
}

// 3e.5 — pattern_keys_noise_robust(): on a noise-free 4-tuple, the returned
// set is exactly the canonical key; when an adjacent-edge gap is below the
// noise tolerance, the set includes BOTH the canonical key and the rank-
// flipped key. Anchors the 24-perm probing behavior end-to-end at the
// key-computation layer (without going through identify_stars).
TEST_F(IdentificationTest, PermutationProbeFindsKey) {
  std::ifstream fs(PATTERN_FILE_20);
  if (!fs.good())
    GTEST_SKIP() << "Pattern catalog not present — run generate_catalog.py";

  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(20.0 * M_PI / 360.0));
  cam.focal_y = cam.focal_x;

  std::array<int, 4> hips{};
  std::array<std::array<double, 3>, 4> v_cam{};
  double R[3][3];
  build_known_4star_scene(*db, cam, R, hips, v_cam);

  // Sanity: noise-free canonical key
  std::array<int, 4> can{};
  uint64_t k_clean = pattern_key_canonical(v_cam, hips, can);
  ASSERT_NE(k_clean, 0ULL);

  // (a) On clean input with a tight noise_tol, the noise-robust set collapses
  // to {k_clean} (no adjacent gaps are within the tolerance).
  auto out_tight = pattern_keys_noise_robust(v_cam, hips, 1e-9);
  ASSERT_FALSE(out_tight.empty());
  EXPECT_EQ(out_tight.front().first, k_clean)
      << "First entry of noise-robust output must be the noise-free key";
  std::set<uint64_t> tight_keys;
  for (const auto &p : out_tight) tight_keys.insert(p.first);
  EXPECT_EQ(tight_keys.size(), 1u)
      << "On a clean input with a 1e-9 noise_tol the set should be {canonical}";

  // (b) Find the SMALLEST adjacent-distance gap across the 6 sorted edges,
  // then run with a noise_tol just above it. This forces at least one rank
  // flip and the set must contain AT LEAST 2 distinct keys; the canonical
  // key must remain in the set as the first entry.
  static constexpr int EDGE_PAIRS[6][2] = {
      {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = EDGE_PAIRS[e][0];
    int j = EDGE_PAIRS[e][1];
    double d = v_cam[i][0] * v_cam[j][0] + v_cam[i][1] * v_cam[j][1] +
               v_cam[i][2] * v_cam[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }
  std::array<double, 6> sorted_dists = dists;
  std::sort(sorted_dists.begin(), sorted_dists.end());
  double min_gap = 1e9;
  for (int i = 0; i < 5; ++i)
    min_gap = std::min(min_gap, sorted_dists[i + 1] - sorted_dists[i]);
  ASSERT_GT(min_gap, 0.0)
      << "Anchor pattern is degenerate — gap of zero between adjacent edges";

  double noise_tol = min_gap * 1.5;
  auto out_loose = pattern_keys_noise_robust(v_cam, hips, noise_tol);
  ASSERT_FALSE(out_loose.empty());
  EXPECT_EQ(out_loose.front().first, k_clean)
      << "Canonical key must remain first in the noise-robust output";
  std::set<uint64_t> loose_keys;
  for (const auto &p : out_loose) loose_keys.insert(p.first);
  EXPECT_GE(loose_keys.size(), 2u)
      << "noise_tol = 1.5x min_gap should force at least one rank flip; got "
      << loose_keys.size() << " unique keys";
  EXPECT_TRUE(loose_keys.count(k_clean))
      << "Canonical key must be in the noise-robust set";
}

// 3e.5 — NoiseRobustnessSweep: identify_stars must keep producing correct
// HIPs as centroid noise climbs from 0″ to 60″. At our 11.5° FOV (focal
// ~5083 px) and 11k px / FOV the per-pixel angle is ~8.1″, so 60″ ≈ 7 px
// jitter — enough to scramble centroids without losing identifiability.
TEST_F(IdentificationTest, NoiseRobustnessSweep) {
  try {
    db->load_pattern_catalog(PATTERN_FILE_10);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "FOV-10 pattern catalog not available: " << e.what();
  }

  PinholeCamera cam_truth;
  cam_truth.frame_width = 1024;
  cam_truth.frame_height = 1024;
  cam_truth.center_x = 512.0;
  cam_truth.center_y = 512.0;
  const double fov_deg = 11.5;
  cam_truth.focal_x =
      1024.0 / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  cam_truth.focal_y = cam_truth.focal_x;

  // Dense star field around HIP 3334 (same anchor as RefinedFovRecoversAlt60Scenario).
  CatalogStar bs = db->get_star(3334);
  double c3[3] = {bs.x, bs.y, bs.z};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
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

  std::vector<int> candidates = {
      3334, 3649, 3058, 3821, 4292, 3030, 2876, 4422, 4383,
      2377, 3584, 3179, 4440, 4427, 4811, 4151, 4961, 2074,
      2101, 1892, 3988, 4475, 5232, 2440, 5361,
  };
  std::vector<StarCentroid> base_centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam_truth, R, base_centroids,
                       expected_hips);
  ASSERT_GE(base_centroids.size(), 15u);

  // Per-pixel angle ≈ FOV_rad / frame_width; 1″ noise = (1/3600°·π/180) rad
  // → roughly 0.123 px at this configuration. We sweep up to 60″ ≈ 7.4 px.
  const std::vector<double> noise_arcsec = {0.0, 5.0, 15.0, 30.0, 60.0};
  std::mt19937 rng(20260511);
  std::set<int> truth(expected_hips.begin(), expected_hips.end());

  for (double sigma_arcsec : noise_arcsec) {
    const double sigma_px = sigma_arcsec * (M_PI / 180.0 / 3600.0) *
                            cam_truth.focal_x;
    std::normal_distribution<double> jitter(0.0, sigma_px);

    std::vector<StarCentroid> noisy = base_centroids;
    for (auto &c : noisy) {
      c.x += jitter(rng);
      c.y += jitter(rng);
    }

    auto identified = identify_stars(noisy, cam_truth, *db, 1e-5);
    int correct = 0;
    for (const auto &s : identified)
      if (truth.count(s.catalog_hip_id)) ++correct;

    // At every sigma we exercise here, the pattern path or the pyramid
    // fallback must recover at least 4 correct stars — enough for QUEST
    // to produce a valid attitude.
    EXPECT_GE(correct, 4)
        << "Noise sweep failed at sigma=" << sigma_arcsec
        << "″ — correct=" << correct << " identified=" << identified.size();
  }
}

// 3e.5 — AccuracyAfterVerify: after the pattern verify + tight-inlier-expand
// + QUEST-refine + re-expand chain (and the outer FOV-scale loop), the
// resulting IdentifiedStar list must yield an attitude within 0.05° of truth
// when fed to estimate_attitude. This anchors the end-to-end accuracy floor
// that 3e.5 was designed to deliver. We exercise it at 0″ and 5″ centroid
// noise; the former is the absolute-best case, the latter is the noise
// regime the Monte Carlo runs at.
TEST_F(IdentificationTest, AccuracyAfterVerify) {
  try {
    db->load_pattern_catalog(PATTERN_FILE_10);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "FOV-10 pattern catalog not available: " << e.what();
  }

  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  const double fov_deg = 11.5;
  cam.focal_x = 1024.0 / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;

  // Build truth rotation: +Z aligned with HIP 3334.
  CatalogStar bs = db->get_star(3334);
  double c3[3] = {bs.x, bs.y, bs.z};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
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

  std::vector<int> candidates = {
      3334, 3649, 3058, 3821, 4292, 3030, 2876, 4422, 4383,
      2377, 3584, 3179, 4440, 4427, 4811, 4151, 4961, 2074,
      2101, 1892, 3988, 4475, 5232, 2440, 5361,
  };
  std::vector<StarCentroid> base_centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam, R, base_centroids, expected_hips);
  ASSERT_GE(base_centroids.size(), 15u);

  // Truth quaternion from R: R is the inertial→camera rotation. Convert via
  // standard quaternion-from-matrix; matched against the same convention used
  // by estimate_attitude (vector part is the rotation axis).
  auto quat_from_R = [](const double M[3][3]) {
    double tr = M[0][0] + M[1][1] + M[2][2];
    Quaternion q;
    if (tr > 0) {
      double s = 0.5 / std::sqrt(tr + 1.0);
      q.w = 0.25 / s;
      q.x = (M[2][1] - M[1][2]) * s;
      q.y = (M[0][2] - M[2][0]) * s;
      q.z = (M[1][0] - M[0][1]) * s;
    } else if (M[0][0] > M[1][1] && M[0][0] > M[2][2]) {
      double s = 2.0 * std::sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]);
      q.w = (M[2][1] - M[1][2]) / s;
      q.x = 0.25 * s;
      q.y = (M[0][1] + M[1][0]) / s;
      q.z = (M[0][2] + M[2][0]) / s;
    } else if (M[1][1] > M[2][2]) {
      double s = 2.0 * std::sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]);
      q.w = (M[0][2] - M[2][0]) / s;
      q.x = (M[0][1] + M[1][0]) / s;
      q.y = 0.25 * s;
      q.z = (M[1][2] + M[2][1]) / s;
    } else {
      double s = 2.0 * std::sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]);
      q.w = (M[1][0] - M[0][1]) / s;
      q.x = (M[0][2] + M[2][0]) / s;
      q.y = (M[1][2] + M[2][1]) / s;
      q.z = 0.25 * s;
    }
    if (q.w < 0) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
    return q;
  };
  Quaternion q_truth = quat_from_R(R);

  auto attitude_err_deg = [](const Quaternion &a, const Quaternion &b) {
    double d = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    if (d < 0) d = -d;
    if (d > 1.0) d = 1.0;
    return 2.0 * std::acos(d) * 180.0 / M_PI;
  };

  const std::vector<double> noise_arcsec = {0.0, 5.0};
  std::mt19937 rng(20260511);
  for (double sigma_arcsec : noise_arcsec) {
    const double sigma_px = sigma_arcsec * (M_PI / 180.0 / 3600.0) * cam.focal_x;
    std::normal_distribution<double> jitter(0.0, sigma_px);

    std::vector<StarCentroid> noisy = base_centroids;
    for (auto &c : noisy) {
      c.x += jitter(rng);
      c.y += jitter(rng);
    }

    auto identified = identify_stars(noisy, cam, *db, 1e-5);
    ASSERT_GE(identified.size(), 4u)
        << "AccuracyAfterVerify needs ≥4 identified stars at sigma="
        << sigma_arcsec << "″";

    Quaternion q_est = estimate_attitude(identified, *db);
    double err = attitude_err_deg(q_est, q_truth);
    EXPECT_LT(err, 0.05)
        << "AccuracyAfterVerify: err=" << err << "° at sigma="
        << sigma_arcsec << "″ exceeds 0.05° gate";
  }
}

// === Phase 3g.7: adversarial unit tests ===

// 3g.7 (1) — EdgeRankFlipStress. Anchor on the same FOV-20 4-tuple as
// PermutationProbeFindsKey, project it to a 20° FOV camera, and add jitter
// in *pixel* space that's enough to flip the canonical sort order on the
// two closest sorted edges. The noise-robust path
// (pattern_keys_noise_robust, called from identify_stars with a 2 mrad
// tolerance) must produce the canonical key under one of the 2^k flipped
// orderings and recover the 4 HIPs end-to-end.
//
// identify_stars's pattern path requires N >= 5 (4 + 1 verification), so we
// pad the anchor 4-tuple with extra catalog stars projected through the
// same rotation. We pull them straight from catalog_stars.bin (whose layout
// is documented in catalog.cpp) so the test stays self-contained.
TEST_F(IdentificationTest, EdgeRankFlipStress) {
  try {
    db->load_pattern_catalog(PATTERN_FILE_20);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "FOV-20 pattern catalog not available: " << e.what();
  }

  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  cam.focal_x = 1024.0 / (2.0 * std::tan(20.0 * M_PI / 360.0));
  cam.focal_y = cam.focal_x;

  std::array<int, 4> hips{};
  std::array<std::array<double, 3>, 4> v_cam{};
  double R[3][3];
  build_known_4star_scene(*db, cam, R, hips, v_cam);

  // Find the smallest adjacent-distance gap across the 6 sorted edges. We'll
  // then nudge centroids by an amount that, in angle terms, just exceeds the
  // gap — guaranteeing a rank flip on the noise-free canonical key.
  static constexpr int EDGE_PAIRS[6][2] = {
      {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = EDGE_PAIRS[e][0];
    int j = EDGE_PAIRS[e][1];
    double d = v_cam[i][0] * v_cam[j][0] + v_cam[i][1] * v_cam[j][1] +
               v_cam[i][2] * v_cam[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }
  std::array<double, 6> sorted_dists = dists;
  std::sort(sorted_dists.begin(), sorted_dists.end());
  double min_gap = 1e9;
  for (int i = 0; i < 5; ++i)
    min_gap = std::min(min_gap, sorted_dists[i + 1] - sorted_dists[i]);
  ASSERT_GT(min_gap, 0.0)
      << "Anchor pattern is degenerate — adjacent edges have zero gap";

  // Confirm the anchor pattern's two closest edges are within the
  // noise-robust trigger threshold (2 mrad inside identify_stars). At the
  // FOV-20 index 435036 anchor used by build_known_4star_scene this is
  // ~1.85e-3 rad, so any centroid noise pushing the per-edge angle by half
  // that amount can flip the canonical sort order. If the anchor ever
  // becomes a "wider gap" pattern this assertion would catch it.
  EXPECT_LT(min_gap, 2e-3)
      << "Anchor pattern's smallest adjacent edge gap is " << min_gap
      << " rad; rank-flip stress is mild and the noise-robust path may "
         "not engage";

  // Project the 4 camera-frame unit vectors to pixel coordinates. v_cam from
  // build_known_4star_scene is already in the camera frame (rotation applied
  // to inertial unit vectors), so projection is the direct pinhole map.
  std::vector<StarCentroid> centroids(4);
  std::set<int> seen_hips;
  for (int i = 0; i < 4; ++i) {
    ASSERT_GT(v_cam[i][2], 0.0) << "Star " << i << " is behind the camera";
    centroids[i].x = cam.focal_x * (v_cam[i][0] / v_cam[i][2]) + cam.center_x;
    centroids[i].y = cam.focal_y * (v_cam[i][1] / v_cam[i][2]) + cam.center_y;
    centroids[i].intensity = 1000.0;
    centroids[i].peak = 250.0 - i; // brightness-rank for seed ordering
    seen_hips.insert(hips[i]);
  }

  // Pad with extra real catalog stars from catalog_stars.bin within the
  // frame, so identify_stars's pattern path (which needs N >= 5) can run.
  // Layout: int32 num_stars, then num_stars × {int32 hip, double x,y,z}.
  // Collect their HIPs in `extra_hips` so the truth set below covers the
  // full real-catalog content of the scene.
  std::vector<int> extra_hips;
  {
    std::ifstream fs(STAR_FILE, std::ios::binary);
    ASSERT_TRUE(fs.good()) << "catalog_stars.bin missing";
    int32_t num_stars = 0;
    fs.read(reinterpret_cast<char *>(&num_stars), 4);
    ASSERT_GT(num_stars, 0);
    int extras_added = 0;
    for (int i = 0; i < num_stars && extras_added < 6; ++i) {
      int32_t hip;
      double x, y, z;
      fs.read(reinterpret_cast<char *>(&hip), 4);
      fs.read(reinterpret_cast<char *>(&x), 8);
      fs.read(reinterpret_cast<char *>(&y), 8);
      fs.read(reinterpret_cast<char *>(&z), 8);
      if (seen_hips.count(hip)) continue;
      double vz = R[2][0] * x + R[2][1] * y + R[2][2] * z;
      if (vz <= 0.5) continue; // must be in front of camera, away from edge
      double vx = R[0][0] * x + R[0][1] * y + R[0][2] * z;
      double vy = R[1][0] * x + R[1][1] * y + R[1][2] * z;
      StarCentroid c;
      c.x = cam.focal_x * (vx / vz) + cam.center_x;
      c.y = cam.focal_y * (vy / vz) + cam.center_y;
      if (c.x < 0 || c.x >= cam.frame_width || c.y < 0 ||
          c.y >= cam.frame_height)
        continue;
      c.intensity = 500.0;
      c.peak = 100.0 - extras_added; // fainter than anchor 4 (peak 247-250)
      centroids.push_back(c);
      extra_hips.push_back(hip);
      extras_added++;
    }
    ASSERT_GE(centroids.size(), 5u)
        << "Couldn't find enough extra catalog stars near the anchor 4-tuple";
  }

  // Add pixel jitter sized to flip rank on the closest-pair edges. At this
  // FOV the per-pixel angle is ~20°/1024 ≈ 3.4e-4 rad, so 1 px of jitter
  // produces ~0.5 mrad of inter-star angle noise — enough to flip rank on
  // the closest edge pair (min_gap ≈ 1.85e-3 rad) under some RNG seeds,
  // while staying inside the post-TRIAD verify gate (4e-7 cosine ≈ 0.05°,
  // which corresponds to ~1 mrad of inertial projection error after TRIAD
  // on the noisy anchor pair). Sample a small grid of RNG seeds and require
  // at least one seed to produce a successful identification: this proves
  // the noise-robust probing reaches the canonical key under at least one
  // realistic rank-flipped ordering.
  // Truth set = anchor 4 HIPs + extra HIPs added to pad N >= 5.
  std::set<int> truth(hips.begin(), hips.end());
  for (int h : extra_hips) truth.insert(h);

  bool any_seed_ok = false;
  int best_identified = 0;
  int best_correct = 0;
  std::normal_distribution<double> jitter(0.0, 1.0);
  for (uint32_t seed : {20260512u, 20260513u, 20260514u, 20260515u, 20260516u}) {
    std::mt19937 rng(seed);
    std::vector<StarCentroid> noisy = centroids;
    // Apply noise to the first 4 centroids (the rank-flip-sensitive anchor
    // 4-tuple). The extra-padding stars stay clean so they provide stable
    // cross-verification anchors for the pipeline.
    for (int i = 0; i < 4; ++i) {
      noisy[i].x += jitter(rng);
      noisy[i].y += jitter(rng);
    }

    auto identified = identify_stars(noisy, cam, *db, 1e-5);
    int correct = 0;
    bool any_anchor = false;
    std::set<int> anchor_set(hips.begin(), hips.end());
    for (const auto &s : identified) {
      if (truth.count(s.catalog_hip_id)) ++correct;
      if (anchor_set.count(s.catalog_hip_id)) any_anchor = true;
    }
    if ((int)identified.size() > best_identified) {
      best_identified = identified.size();
      best_correct = correct;
    }
    // Require at least 4 correct identifications AND at least one of them
    // to be from the anchor 4-tuple (i.e. the rank-flipped pattern actually
    // landed on the canonical key under the noise-robust probing).
    if ((int)identified.size() >= 4 && correct >= 4 && any_anchor) {
      any_seed_ok = true;
      break;
    }
  }
  EXPECT_TRUE(any_seed_ok)
      << "EdgeRankFlipStress: noise-robust probing failed under all RNG "
         "seeds. Best run: identified=" << best_identified
      << " correct=" << best_correct
      << ". The pipeline did not recover the rank-flipped canonical key.";
}

// 3g.7 (2) — RejectsBrightNonCatalogContaminant. Build the Ursa Major scene,
// then insert one ultra-bright centroid at a pixel location that doesn't
// correspond to any real catalog star. The contaminant has the largest
// `peak`, so the seed-selection loop in identify_stars will try it first.
// Identification must reject the contaminant (it won't appear in the
// returned list) while still recovering the real stars.
TEST_F(IdentificationTest, RejectsBrightNonCatalogContaminant) {
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
  ASSERT_GE(centroids.size(), 8u) << "Need clean scene of ≥8 stars";

  // Real-star centroids get descending peaks so they have a stable
  // brightness order independent of insertion order.
  for (size_t i = 0; i < centroids.size(); ++i)
    centroids[i].peak = 200.0 - static_cast<double>(i);

  // Insert a non-catalog "bright object" at a pixel position chosen to not
  // sit within ~5° of any of the real stars. The corner (50, 950) is well
  // off-pattern for our Ursa-Major-centric +Z boresight at 90° FOV.
  StarCentroid contaminant;
  contaminant.x = 50.0;
  contaminant.y = 950.0;
  contaminant.intensity = 10000.0;
  contaminant.peak = 500.0; // largest peak → seed candidate index 0
  int contaminant_idx = static_cast<int>(centroids.size());
  centroids.push_back(contaminant);

  // Sanity: contaminant is at least 5° from each real star's unit-vector
  // direction. We back-project both to camera-frame unit vectors and check
  // their dot product against cos(5°).
  double v_cont[3];
  undistort_to_unit_vector(cam, contaminant.x, contaminant.y, v_cont);
  const double cos5 = std::cos(5.0 * M_PI / 180.0);
  for (size_t i = 0; i < centroids.size() - 1; ++i) {
    double v_real[3];
    undistort_to_unit_vector(cam, centroids[i].x, centroids[i].y, v_real);
    double d = v_cont[0] * v_real[0] + v_cont[1] * v_real[1] +
               v_cont[2] * v_real[2];
    ASSERT_LT(d, cos5)
        << "Contaminant is within 5° of real star " << i;
  }

  const double cos_tol = 1e-5;
  auto identified = identify_stars(centroids, cam, *db, cos_tol);

  // The contaminant centroid must not appear in the result.
  for (const auto &s : identified) {
    EXPECT_NE(s.image_idx, contaminant_idx)
        << "Non-catalog contaminant was matched to HIP " << s.catalog_hip_id;
  }

  // Every returned HIP must be a real catalog HIP from our truth set.
  std::set<int> truth(expected_hips.begin(), expected_hips.end());
  for (const auto &s : identified) {
    EXPECT_TRUE(truth.count(s.catalog_hip_id))
        << "Identified an unexpected (non-truth) HIP " << s.catalog_hip_id
        << " for image_idx=" << s.image_idx;
  }

  // The real stars should still be recoverable despite the contaminant.
  int correct = 0;
  for (const auto &s : identified)
    if (truth.count(s.catalog_hip_id)) ++correct;
  EXPECT_GE(correct, 6)
      << "Contaminant disrupted identification of real stars; correct="
      << correct;
}

// 3g.7 (3a) — FovBinBoundary12_5. main.cpp's FOV-to-pattern-bin rule is
// (fov_deg <= 12.5) ? 10 : (fov_deg <= 17.5) ? 15 : 20. 12.5° lands in bin
// 10; 17.5° lands in bin 15 (covered by the next test). Build a dense
// synthetic scene at exactly 12.5° FOV and confirm identification succeeds.
TEST_F(IdentificationTest, FovBinBoundary12_5) {
  // Mirror main.cpp's bin selection: 12.5° → bin 10.
  try {
    db->load_pattern_catalog(PATTERN_FILE_10);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "FOV-10 pattern catalog not available: " << e.what();
  }

  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  const double fov_deg = 12.5;
  cam.focal_x = 1024.0 / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;

  // Reuse the dense HIP-3334 anchor used by RefinedFovRecoversAlt60Scenario
  // and NoiseRobustnessSweep. ~25 candidates within 6° of boresight.
  CatalogStar bs = db->get_star(3334);
  double c3[3] = {bs.x, bs.y, bs.z};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
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

  std::vector<int> candidates = {
      3334, 3649, 3058, 3821, 4292, 3030, 2876, 4422, 4383,
      2377, 3584, 3179, 4440, 4427, 4811, 4151, 4961, 2074,
      2101, 1892, 3988, 4475, 5232, 2440, 5361,
  };
  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam, R, centroids, expected_hips);
  ASSERT_GE(centroids.size(), 6u)
      << "Need ≥6 stars in 12.5° FOV around HIP 3334 to exercise the pattern path";

  auto identified = identify_stars(centroids, cam, *db, 1e-5);
  ASSERT_GE(identified.size(), 4u)
      << "FovBinBoundary12_5: identification produced fewer than 4 matches "
      << "(got " << identified.size() << ")";

  std::set<int> truth(expected_hips.begin(), expected_hips.end());
  int correct = 0;
  for (const auto &s : identified)
    if (truth.count(s.catalog_hip_id)) ++correct;
  EXPECT_GE(correct, 4)
      << "FovBinBoundary12_5: " << correct << " correct of "
      << identified.size() << " identified";
}

// 3g.7 (3b) — FovBinBoundary17_5. 17.5° → bin 15. The FOV-15 pattern catalog
// isn't referenced elsewhere by these tests; we use the path literally.
TEST_F(IdentificationTest, FovBinBoundary17_5) {
  const char *PATTERN_FILE_15 = "../data/catalog_patterns_15.bin";
  try {
    db->load_pattern_catalog(PATTERN_FILE_15);
  } catch (const std::exception &e) {
    GTEST_SKIP() << "FOV-15 pattern catalog not available: " << e.what();
  }

  PinholeCamera cam;
  cam.frame_width = 1024;
  cam.frame_height = 1024;
  cam.center_x = 512.0;
  cam.center_y = 512.0;
  const double fov_deg = 17.5;
  cam.focal_x = 1024.0 / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  cam.focal_y = cam.focal_x;

  // Same dense HIP-3334 anchor; at 17.5° FOV we capture even more candidates.
  CatalogStar bs = db->get_star(3334);
  double c3[3] = {bs.x, bs.y, bs.z};
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(c3[0]) < 1e-6 && std::abs(c3[1]) < 1e-6) {
    up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
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

  std::vector<int> candidates = {
      3334, 3649, 3058, 3821, 4292, 3030, 2876, 4422, 4383,
      2377, 3584, 3179, 4440, 4427, 4811, 4151, 4961, 2074,
      2101, 1892, 3988, 4475, 5232, 2440, 5361,
  };
  std::vector<StarCentroid> centroids;
  std::vector<int> expected_hips;
  gather_visible_stars(candidates, *db, cam, R, centroids, expected_hips);
  ASSERT_GE(centroids.size(), 8u)
      << "Need ≥8 stars in 17.5° FOV around HIP 3334";

  auto identified = identify_stars(centroids, cam, *db, 1e-5);
  ASSERT_GE(identified.size(), 4u)
      << "FovBinBoundary17_5: identification produced fewer than 4 matches "
      << "(got " << identified.size() << ")";

  std::set<int> truth(expected_hips.begin(), expected_hips.end());
  int correct = 0;
  for (const auto &s : identified)
    if (truth.count(s.catalog_hip_id)) ++correct;
  EXPECT_GE(correct, 4)
      << "FovBinBoundary17_5: " << correct << " correct of "
      << identified.size() << " identified";
}
