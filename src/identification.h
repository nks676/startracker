#pragma once
#include "camera_model.h"
#include "catalog.h"
#include "image_processing.h"
// pattern_hash.h exposes pattern_key_canonical / pattern_keys_noise_robust at
// file scope so existing test_identification.cpp callers keep compiling
// without any namespace prefix (Phase 3g.1 file split).
#include "pattern_hash.h"
#include <array>
#include <cstdint>
#include <vector>

struct IdentifiedStar {
  int image_idx;
  int catalog_hip_id;
  double v_cam[3]; // Unit vector in camera frame
};

// Phase 3g.3: optional per-call telemetry, populated by identify_stars when
// a non-null pointer is passed. All counters reset to defaults on entry; the
// fields are independent so the caller can rely on whichever is meaningful
// for the path that actually executed (e.g. fov_scale_iters is 0 if the
// pyramid fallback ran).
struct IdentifyStats {
  // True if the catalog had a pattern hash table loaded for this call.
  bool pattern_catalog_loaded = false;
  // Seed centroids the pattern path tried before either returning a verified
  // assignment or exhausting the pool. 0 if pattern path didn't run.
  int pattern_seeds_tried = 0;
  // True if a pattern-path seed produced a verified assignment that the
  // pipeline returned. False on pattern-path miss or when N < 5.
  bool pattern_path_hit = false;
  // True if the pyramid fallback ran (pattern path failed or N < 5).
  bool fallback_to_pyramid = false;
  // Final identified-star count (== returned vector size).
  int final_inliers = 0;
  // Number of FOV-scale absorption iterations the post-pattern-path loop ran
  // (3g.3: useful for diagnosing camera-calibration drift on Pi captures).
  int fov_scale_iters = 0;
  // Newton iterations the public-path QUEST used (0 if estimate_attitude
  // isn't reached, e.g. identify_stars itself returned an empty list).
  // Populated by estimate_attitude when called with stats; identify_stars
  // does not call estimate_attitude itself.
  int quest_newton_iters = 0;
};

// Identifies stars in the image against the catalog.
//
// `stats_out`, when non-null, receives per-call telemetry (Phase 3g.3). Pass
// `nullptr` (the default) for the no-overhead path used by everything except
// the `--stats` CLI flag.
std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance, IdentifyStats *stats_out = nullptr);
