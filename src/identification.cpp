#include "identification.h"
#include "inlier_expand.h"
#include "pattern_hash.h"
#include "pyramid.h"
#include "wahba.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>

// Phase 3g.1 (completion): identification.cpp now hosts just the orchestrator
// — the pattern-path `identify_stars_pattern` and the top-level
// `identify_stars` entry. The pyramid-path helpers (run_pyramid,
// estimate_scale_factor, cross_verify, identify_stars_pyramid) live in
// pyramid.h; the inlier-expansion + pattern-verify helpers
// (top_n_indices, k_nearest_within_radius, expand_inliers,
// expand_inliers_tight, refine_and_reexpand, try_verify_candidate,
// kFifthStarVerifyCosTol) live in inlier_expand.h. Both header-only because
// every entry point is templated on CatVec.

namespace {

// Recover the horizontal half-FOV (radians) from the camera model. Used as
// the upper bound on the 4-star pattern radius — the pattern catalog at FOV
// bin = 2·half-FOV stores only patterns with all 4 stars within that radius
// of the seed. Stays here because identify_stars_pattern (below) is the only
// user.
inline double camera_half_fov_rad(const CameraModel &camera) {
  return std::atan(camera.frame_width / (2.0 * camera.focal_x));
}

// === Phase 3e.3: pattern-path identify_stars implementation. ===
//
// Returns the (filled) assignment image_idx -> HIP, or an empty vector if no
// pattern attempt produced a verified attitude. On success, the assignment
// covers the 4 pattern stars + the verified 5th star (other stars may still
// be -1; caller / expansion step extends from there).
//
// Logs "[3e] pattern path failed → pyramid fallback" if no candidate
// verifies after the first kPatternAttempts seed centroids — that lets us
// count fallback frequency in production.

// Number of seed centroids we try (brightest N by peak). Originally 4 per
// plan 3e.3. Diagnostics in the 3e.5 follow-up (DEBUG_PATTERN_STATS) showed
// the failure mode on noisy synthetic scenes: when one of the top-4 brightest
// centroids isn't a catalog star (or its nearest detected neighbors aren't a
// subset of its catalog top-8 neighbors), all 4 seeds can fail in lock-step
// even though dimmer-seed alternatives would succeed. Bumping to 10 raises
// the per-frame pattern-path hit rate substantially on Monte Carlo without
// measurably increasing the success-case latency (the path returns on first
// verify, so extra seeds are only exercised on the fallback boundary).
constexpr int kPatternAttempts = 10;
// Per-seed, enumerate C(kPatternKNearest, 3) triples drawn from the seed's
// kPatternKNearest nearest in-radius neighbors. The catalog generator stores
// patterns built from each catalog star's top-8 nearest within FOV/2, so
// kPatternKNearest=8 mirrors the catalog's coverage. We bump to 10 to
// absorb cases where centroid noise / vmag-cutoff shifts which neighbor
// ranks highest.
constexpr int kPatternKNearest = 10;
// Verification pool: where the 5th star is hunted. The seed's 3 nearest
// neighbors are drawn from ALL detected centroids (not just bright ones) —
// the catalog only contains 4-tuples within FOV/2 of the seed, and at wide
// FOVs the brightest stars are sparse enough that a bright-only pool can
// be empty inside that radius. The verify pool is intentionally narrower
// because verification should anchor on real stars (not noise).
constexpr int kVerifyPoolSize = 12;
// Phase 3g.1: kFifthStarVerifyCosTol moved to inlier_expand.h (used both
// inside try_verify_candidate / refine_and_reexpand and below in
// identify_stars's post-pattern FOV-scale loop, so it now lives where the
// helpers live).
// Inlier-expansion acceptance threshold. Same 0.05° gate the brief mandates
// for the post-verify inlier expansion step (cos(0.05°) ≈ 1 - 3.8e-7). Wide
// enough to absorb sub-degree FOV miscal that the TRIAD seed attitude can't
// correct, while still rejecting wrong-star matches (which would land 10×
// further). The follow-up FOV-scale absorption in identify_stars() tightens
// the implied tolerance once the camera focal is rescaled.
constexpr double kInlierExpansionCosTol = 4e-7;
// Pair-consistency cosine tolerance for the 4 pattern stars. Cos noise
// across a 5-20° pattern at our centroid noise level is ~1e-5 — we use 5e-4
// which still rejects gross mismatches but is loose enough to absorb FOV
// calibration drift (1e-3 was the alt60 worst case).
constexpr double kPatternPairCosTol = 5e-4;
// Safety margin on the catalog FOV radius: the catalog stores patterns
// within FOV_bin/2 of each seed. We don't know FOV_bin directly (catalog.h
// doesn't expose it), but the binary picks FOV_bin to bracket the camera
// FOV. For real-image fixtures the camera FOV is ~11° and bin=10° (radius
// 5°); for FOV=20° MC bin=20° (radius 10°). cos(half_FOV) is a safe upper
// bound on the catalog radius — and any quad of stars within that radius is
// either in the catalog or geometrically equivalent to one that is. We pull
// in a small margin (95%) to absorb FOV-bin/camera-FOV mismatch.
constexpr double kPatternRadiusSafetyMargin = 0.95;

template <typename CatVec>
std::vector<int>
identify_stars_pattern(const std::vector<StarCentroid> &image_stars,
                       const PinholeCamera &camera,
                       const std::vector<std::array<double, 3>> &v_cam,
                       const StarDatabase &db,
                       const CatVec &cat_vec,
                       IdentifyStats *stats = nullptr) {
  const int N = static_cast<int>(image_stars.size());
  if (N < 5)
    return {}; // pattern path requires 4 + 1 verification star

  // Pattern catalogs only contain 4-star tuples with all 4 stars within
  // FOV_bin/2 of the brightest. Use the camera half-FOV as a proxy for
  // FOV_bin/2 (main.cpp picks FOV_bin to bracket the camera FOV).
  const double half_fov_rad = camera_half_fov_rad(camera) *
                              kPatternRadiusSafetyMargin;
  const double cos_radius = std::cos(half_fov_rad);

  // Rank all centroids by peak descending. The first kSeedPoolSize are seed
  // candidates; the full ranked list is the neighbor pool (the seed's 3
  // nearest within FOV/2 may be a mag-5 star next to the bright seed, not
  // necessarily another mag-2 star). Verification pool is the top-N
  // brightest so we don't anchor 5th-star verification on faint noise.
  auto ranked = top_n_indices(image_stars, static_cast<int>(image_stars.size()));
  if (static_cast<int>(ranked.size()) < 4) return {};

  // Neighbor pool = ALL detected centroids (peak-ordered, brightest first).
  const std::vector<int> &neighbor_pool = ranked;
  std::vector<int> verify_pool(
      ranked.begin(),
      ranked.begin() + std::min<int>(ranked.size(), kVerifyPoolSize));

  const int seeds = std::min(kPatternAttempts, static_cast<int>(ranked.size()));
  for (int s = 0; s < seeds; ++s) {
    if (stats) stats->pattern_seeds_tried = s + 1;
    int seed_idx = ranked[s];

    // Find the seed's nearest in-radius neighbors. The catalog stores
    // patterns built from each catalog star's top-8 nearest within FOV/2;
    // mirror that on the observation side by enumerating triples of the
    // seed's kPatternKNearest nearest in-radius detected centroids.
    auto knn = k_nearest_within_radius(seed_idx, neighbor_pool, v_cam,
                                       cos_radius, kPatternKNearest);
    if (knn.size() < 3) continue; // not enough neighbors → next seed

    // Enumerate C(knn.size(), 3) triples. With kPatternKNearest=10 this is
    // up to 120 hashes per seed, but most catalog probes return zero or
    // few candidates — the hot path is the hash + tolerant probe, both
    // cheap.
    const int K = static_cast<int>(knn.size());
    for (int a = 0; a < K - 2; ++a) {
      for (int b = a + 1; b < K - 1; ++b) {
        for (int c = b + 1; c < K; ++c) {
          std::array<int, 4> centroid_indices = {seed_idx, knn[a], knn[b],
                                                  knn[c]};

          // Phase 3e.5 (Change 1): noise-robust multi-key probing.
          //
          // The canonical-order key sorts the 6 inter-star angles. Under
          // centroid noise, two near-equal angles can flip rank, shifting the
          // canonical permutation and bumping the key to a distant bucket
          // that ±1 tolerant probing can't reach. The catalog only stores the
          // canonical key for the noise-free geometry, so we enumerate the
          // 2^k orderings produced by flipping each uncertain adjacent-rank
          // pair (where k is the number of uncertain pairs ≤ 5; in practice
          // 1–3). One of those orderings will match the catalog's canonical
          // key. Probe each with find_pattern_tolerant for residual ±1
          // quantization slack.
          //
          // The 24-input-permutation variant tetra3 uses doesn't apply to our
          // signature-based canonical_order: input permutation is invariant
          // by construction here (per-star signatures are geometry-only, and
          // both catalog and query tie-break on the input identifier — HIP
          // and centroid index respectively — so all 24 input permutations
          // produce the same canonical key). The rank-flip mechanism in
          // pattern_keys_noise_robust is what actually addresses the noise
          // failure mode.
          //
          // noise_tol = 2 mrad ≈ 7 arcmin. At synthetic 5-px-noise input on a
          // 5000 px focal, per-pair angular noise is ~1.4e-3 rad; 1.5σ ≈ 2e-3
          // catches the common rank-flip cases without generating an
          // unreasonably large alternate-key set. Real images (lower per-star
          // noise) end up with the same key as the canonical, which is the
          // single-key path.
          std::array<std::array<double, 3>, 4> vecs4{};
          std::array<int, 4> ids4{};
          for (int i = 0; i < 4; ++i) {
            vecs4[i] = v_cam[centroid_indices[i]];
            ids4[i] = centroid_indices[i];
          }
          constexpr double kPatternRankFlipNoiseRad = 2e-3;
          auto raw_keys =
              pattern_keys_noise_robust(vecs4, ids4, kPatternRankFlipNoiseRad);
          std::vector<std::pair<uint64_t, std::array<int, 4>>> probe_keys;
          probe_keys.reserve(raw_keys.size());
          for (const auto &kp : raw_keys) {
            std::array<int, 4> centroid_canonical{};
            for (int i = 0; i < 4; ++i)
              centroid_canonical[i] = centroid_indices[kp.second[i]];
            probe_keys.emplace_back(kp.first, centroid_canonical);
          }
          for (const auto &kp : probe_keys) {
            uint64_t key = kp.first;
            const auto &centroid_canonical = kp.second;
            auto candidates = db.find_pattern_tolerant(key);
            if (candidates.empty()) continue;

            for (const auto &cand : candidates) {
              std::array<int, 4> hips = {cand.hips[0], cand.hips[1],
                                          cand.hips[2], cand.hips[3]};
              double R_triad[3][3] = {{0}};
              auto assignment = try_verify_candidate(
                  centroid_canonical, hips, verify_pool, v_cam, db, cat_vec,
                  kPatternPairCosTol, kFifthStarVerifyCosTol, R_triad);
              if (assignment.empty()) continue;

              // Change 2: tight inlier expansion + QUEST refine. Projects all
              // remaining centroids into inertial using the TRIAD attitude,
              // matches them against the catalog with a 0.05° gate, then
              // refines the attitude via QUEST and re-expands. Bumps the
              // inlier list from 5 (TRIAD on 4 + 5th-star verify) to N
              // (typically 8–25 on a real scene), which lets the downstream
              // QUEST in estimation.cpp converge to ~0.001° accuracy.
              int n_before = 0;
              for (int ii = 0; ii < (int)assignment.size(); ++ii)
                if (assignment[ii] >= 0) ++n_before;
              int n_after = refine_and_reexpand(assignment, v_cam, db, cat_vec,
                                                R_triad, kInlierExpansionCosTol);
              if (std::getenv("STARTRACKER_DEBUG_3E5")) {
                std::fprintf(stderr,
                             "[3e.5] pattern verify: %d -> %d inliers\n",
                             n_before, n_after);
                for (int ii = 0; ii < (int)assignment.size(); ++ii) {
                  if (assignment[ii] >= 0)
                    std::fprintf(stderr, "  centroid %d -> HIP %d\n", ii,
                                 assignment[ii]);
                }
              }
              return assignment;
            }
          }
        }
      }
    }
  }

  return {};
}

// Phase 3g.1 (completion): identify_stars_pyramid moved to pyramid.h (also
// header-only template for the CatVec parameter).

} // namespace

std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance, IdentifyStats *stats_out) {
  // Phase 3g.3: reset out-stats up front so callers don't see leftover state
  // from a prior call if identify_stars early-returns below.
  if (stats_out) {
    *stats_out = IdentifyStats{};
    stats_out->pattern_catalog_loaded = db.has_pattern_catalog();
  }
  int N = image_stars.size();
  if (N < 3)
    return {};

  // === Camera-frame unit vectors ===
  std::vector<std::array<double, 3>> v_cam =
      pyramid_detail::compute_v_cam(image_stars, camera);

  // Cache catalog unit vectors (one per HIP) to avoid repeated map lookups.
  //
  // Phase 3g.2: lambda is passed straight through to the run_pyramid /
  // identify_stars_pattern / etc. helpers as a template parameter, avoiding
  // the std::function indirection (heap-allocated control block + virtual
  // dispatch) that previously sat on every cat_vec invocation in the hot
  // path. The helpers downstream see the lambda's anonymous closure type
  // directly.
  std::unordered_map<int, std::array<double, 3>> cat_vec_cache;
  auto cat_vec = [&](int hip) -> std::array<double, 3> {
    auto it = cat_vec_cache.find(hip);
    if (it != cat_vec_cache.end())
      return it->second;
    CatalogStar s = db.get_star(hip);
    std::array<double, 3> v = {s.x, s.y, s.z};
    cat_vec_cache[hip] = v;
    return v;
  };

  // === Phase 3e.3: pattern-hash path ===
  //
  // Try the pattern path first if the catalog is loaded. Pattern path is
  // O(K * candidates_per_key) where K ≤ 4 seed attempts × 243 probes; each
  // hit is one TRIAD + a partner search. Compared to the pyramid's O(N²)
  // seed enumeration over N≤50 centroids, this is ~3 orders of magnitude
  // faster when it succeeds. Falls back to pyramid on miss.
  if (db.has_pattern_catalog() && N >= 5) {
    auto pattern_assignment =
        identify_stars_pattern(image_stars, camera, v_cam, db, cat_vec,
                               stats_out);
    if (!pattern_assignment.empty()) {
      // The pattern path already performed tight inlier expansion + QUEST
      // refine + re-expand inside identify_stars_pattern (Phase 3e.5 change
      // 2). All inliers in the returned assignment are within 0.05° of their
      // catalog star against the refined attitude, which is tighter than
      // cross_verify's 3*cos_tolerance gate; running cross_verify here would
      // over-reject on miscalibrated cameras.

      // Phase 3e.5: FOV-scale absorption. Mirrors the pyramid path's
      // coarse-refine-reidentify step (3b.0b). The inlier set gives a stable
      // estimate of the focal scale factor s = true_focal / assumed_focal.
      // Real cameras commonly drift 0.1–0.5% from their nominal FOV; the
      // pattern path's TRIAD/QUEST attitude inherits that drift, capping
      // accuracy at ~(s-1) × FOV/2. Rescaling and re-expanding lets the
      // attitude break through that floor and reach centroid-noise limits.
      //
      // Apply iteratively: each pass collapses ~90% of the residual scale
      // error, so 2–3 passes converge to within the noise floor. We stop
      // once the per-pass adjustment falls below the centroid-noise floor
      // (~1e-4 ratio) or we hit the iteration cap.
      CameraModel refined = camera;
      double cumulative_scale = 1.0;
      for (int iter = 0; iter < 3; ++iter) {
        if (stats_out) stats_out->fov_scale_iters = iter + 1;
        double scale = estimate_scale_factor(pattern_assignment, v_cam, db,
                                              cat_vec);
        if (!std::isfinite(scale)) break;
        // Trigger only when the per-pass residual is well above per-star
        // centroid noise (~1e-4 for our focal). Below that we'd just inject
        // bias rather than remove it.
        if (std::abs(scale - 1.0) <= 1e-4) break;

        refined.focal_x *= scale;
        refined.focal_y *= scale;
        cumulative_scale *= scale;
        auto v_cam_refined = pyramid_detail::compute_v_cam(image_stars, refined);

        // Re-expand inliers under the refined camera. Solve TRIAD on the two
        // matched stars with the widest catalog separation, then run the
        // tight-expansion + QUEST-refine loop. The same routine that built
        // the initial inlier set is reused; under a properly scaled camera,
        // it typically converges to all visible cataloged stars.
        std::vector<int> refined_assignment(N, -1);
        int n_matched = 0;
        for (int i = 0; i < N; ++i) {
          if (pattern_assignment[i] < 0) continue;
          refined_assignment[i] = pattern_assignment[i];
          ++n_matched;
        }
        if (n_matched < 2) break;
        // Build a seed TRIAD attitude from the two matched stars with the
        // widest catalog separation (best conditioning).
        int best_i = -1, best_j = -1;
        double smallest_cos = 2.0;
        std::vector<int> matched_idx;
        matched_idx.reserve(n_matched);
        for (int i = 0; i < N; ++i)
          if (refined_assignment[i] >= 0) matched_idx.push_back(i);
        for (size_t ai = 0; ai < matched_idx.size(); ++ai) {
          auto vi = cat_vec(refined_assignment[matched_idx[ai]]);
          for (size_t aj = ai + 1; aj < matched_idx.size(); ++aj) {
            auto vj = cat_vec(refined_assignment[matched_idx[aj]]);
            double cc = vi[0] * vj[0] + vi[1] * vj[1] + vi[2] * vj[2];
            if (cc < smallest_cos) {
              smallest_cos = cc;
              best_i = matched_idx[ai];
              best_j = matched_idx[aj];
            }
          }
        }
        if (best_i < 0 || best_j < 0) break;
        std::array<double, 3> W1 = v_cam_refined[best_i];
        std::array<double, 3> W2 = v_cam_refined[best_j];
        std::array<double, 3> V1 = cat_vec(refined_assignment[best_i]);
        std::array<double, 3> V2 = cat_vec(refined_assignment[best_j]);
        double R_seed[3][3];
        wahba::triad(W1, W2, V1, V2, R_seed);
        refine_and_reexpand(refined_assignment, v_cam_refined, db, cat_vec,
                             R_seed, kFifthStarVerifyCosTol);

        // Count inliers under refined camera; commit to refined if it didn't
        // strictly lose ground.
        int refined_n = 0;
        for (int i = 0; i < N; ++i)
          if (refined_assignment[i] >= 0) ++refined_n;
        if (refined_n < n_matched) break; // refinement removed inliers — back out
        if (std::getenv("STARTRACKER_DEBUG_3E5")) {
          std::fprintf(stderr,
                       "[3e.5] FOV iter %d: s=%.6f (cumulative %.6f), "
                       "%d -> %d inliers\n",
                       iter, scale, cumulative_scale, n_matched, refined_n);
        }
        pattern_assignment = std::move(refined_assignment);
        v_cam = std::move(v_cam_refined);
      }

      int inliers = 0;
      for (int i = 0; i < N; ++i)
        if (pattern_assignment[i] >= 0) ++inliers;
      if (inliers >= 4) {
        std::vector<IdentifiedStar> identified;
        identified.reserve(static_cast<size_t>(inliers));
        for (int i = 0; i < N; ++i) {
          if (pattern_assignment[i] < 0) continue;
          IdentifiedStar is;
          is.image_idx = i;
          is.catalog_hip_id = pattern_assignment[i];
          is.v_cam[0] = v_cam[i][0];
          is.v_cam[1] = v_cam[i][1];
          is.v_cam[2] = v_cam[i][2];
          identified.push_back(is);
        }
        if (stats_out) {
          stats_out->pattern_path_hit = true;
          stats_out->final_inliers = static_cast<int>(identified.size());
        }
        return identified;
      }
      // Fall through to pyramid — pattern path matched but produced too
      // few inliers; treat as a soft fallback.
    }
    std::fprintf(stderr, "[3e] pattern path failed → pyramid fallback\n");
  }

  // === Fallback: pyramid + coarse-refine-reidentify + cross-verify ===
  if (stats_out) stats_out->fallback_to_pyramid = true;
  auto pyramid_result = identify_stars_pyramid(
      image_stars, camera, db, cos_tolerance, std::move(v_cam), cat_vec);
  if (stats_out)
    stats_out->final_inliers = static_cast<int>(pyramid_result.size());
  return pyramid_result;
}
