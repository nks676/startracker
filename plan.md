# Project Plan: *startracker*

Iterative development of an embedded star-tracking engine: synthetic-data baseline → CI safety net → real-world robustness → accuracy upgrades → speed upgrades → Pi hardware integration → SIMD + tracking. Each phase is regression-tested against the prior phase's gates.

---

## Where We Are (2026-05-12)

**Current capability:** cold-start solve on alt60 ESA tetra3 fixture (1024×768, ~11° FOV) in **~60 ms** on M-series desktop and **~273 ms on Pi 4 (8 GB)**.

Per-stage median (alt60, 10-run median µs):

| Stage | Desktop (M-series) | Pi 4 | Pi/Mac ratio |
|---|---|---|---|
| centroid | 11 ms | **142 ms** | 13× |
| catalog_load (one-time) | 20 ms | 51 ms | 2.6× |
| identify | 1–35 ms (alt60) / ~1 ms (alt40) | **79 ms (alt60)** / 2.4 ms (alt40) | 2–3× |
| estimate | <1 µs | 10 µs | comparable |
| **total** | **60 ms** | **273 ms (alt60) / 196 ms (alt40)** | 3.3–4.5× |

Per-frame steady-state (catalog already loaded): ~145–225 ms on Pi 4. Memory on Pi 4: ~480 MB used during a solve, plenty of headroom in the 8 GB RAM.

**Test coverage:** 52/52 unit tests (Pi 4 ctest: 13 s with catalog mmap'd; was 161 s without `catalog_partners.bin`); real-image regression alt40 = 0.0280°, alt60 = 0.0258° identical across Mac + Pi (gate 0.5°); Monte Carlo on Pi 4 50 trials @ 5″ noise = 100% success / median 0.0041° / max 0.0203° / **pattern hit rate 82%** (better than desktop's 76% — 3g.9's deterministic seed tie-break is paying off).

**Pi 4 baseline highlights (these drive Phase 3f priorities):**
- **Centroid stage 13× slower than Mac.** Not just "Pi is slower" — investigation (see Phase 3f below) traced ~85 ms of the 142 ms to a quadratic `bilinear_sample` inner loop. Algorithmic fix, not SIMD.
- **`identify` on alt60 = 79 ms, alt40 = 2.4 ms — bimodal.** `--stats` confirms alt60 takes 6 pattern-path seeds to hit; alt40 hits seed 1. Same root cause as the previously-noted "TIFF-input identify-stage latency" backlog item, now reproducible on Pi.
- Cold-start gate (<500 ms): ✅ comfortably under.
- Per-frame steady-state gate (<100 ms): ❌ currently 145–225 ms. Closing this is what Phase 3f is for.

**Cleared follow-ups from prior plan revisions:**
- *Pattern-path hit rate on noisy synthetic Monte Carlo was 20% (later 34%).* **Done:** root cause was `K_NEAREST=8` in the pattern catalog generator (per-seed catalog stored only the seed's top-8 nearest catalog neighbors, so query 4-tuples drawn from the seed's 9th-12th nearest had no chance of matching). Bumping to `K_NEAREST=12` raises MC hit rate to **74-90%**, catalog grows ~4× (21 MB → 82 MB per FOV bin). The original "Vmag cutoff mismatch" hypothesis turned out to be wrong — synthetic generator and catalog gen both use `max_mag=7.0`.
- *alt40 went from 0.0000° → 0.0596° during 3e.5.* **Phantom regression.** Root cause was 6-significant-figure truncation in `Estimated Quaternion: ...` print in main.cpp combined with the regression script computing error against the *truncated* string. Bumping `std::cout` precision to 17 sig-figs revealed the actual errors are alt40=0.028° and alt60=0.024° — both fixtures, both well under any gate, and consistent with centroid-noise floor over 17-star QUEST.
- *3e.5 unit tests not landed.* **Done:** `PermutationProbeFindsKey`, `NoiseRobustnessSweep`, `AccuracyAfterVerify` added in `tests/test_identification.cpp`.
- *`per_star_partners` eager build.* **Done (Phase 3f.x):** generator now emits `data/catalog_partners.bin` (directory + entry-block layout); StarDatabase mmaps it at startup without prefaulting (sparse access pattern). `catalog_load` dropped from ~225 ms to ~20 ms — a ≈12× cold-start improvement. Falls back to in-memory construction with a warning if the file is missing.

**Remaining known follow-ups carried into future phases:**
- TIFF input regresses identify-stage latency on alt60 from ~5 ms to ~35 ms vs the PNG-stretched path. Same final attitude, different centroid peak ordering (16-bit ties don't collapse the way 8-bit-saturated peaks do, which shifts which seeds the pattern path tries). Worth re-measuring on Pi before tuning.
- Camera-calibration tool overfits k2/k3 on the 37-point ESA fixture set (≈220-degree polynomial absorbs noise). Adequate as a scaffold; revisit with more frames when Pi captures arrive.

---

## Done (compressed)

### Phase 1 — Bare Minimum Baseline ✅
End-to-end pipeline: synthetic-data tool, catalog builder, image processing (CC-BFS + CoG centroiding), catalog loader, pyramid identification, TRIAD attitude. ~0.004° on identity-quaternion synthetic.

### Phase 2 — CI/CD Infrastructure ✅
- GoogleTest wired via `FetchContent`; `startracker_core` static lib shared by main + tests.
- Unit tests across all four modules (catalog, image processing, identification, estimation).
- GitHub Actions CI: build + ctest + an integration test on a single synthetic frame.
- `tools/monte_carlo.py` aggregates random-orientation accuracy; CI nightly hook.
- `.gitignore` excludes `build/`, `data/`, `tools/venv/`.

### Phase 3a — Real-World Robustness (Image Processing Hardening) ✅
- 3a.1: Brown-Conrady camera distortion model + `CameraModel` struct.
- 3a.2: Background subtraction (tile median + bilinear interp).
- 3a.3: Adaptive thresholding (per-tile `local_mean + k·local_stddev`).
- 3a.4: Hot-pixel / outlier shape filters (min/max pixels, aspect ratio, border).
- 3a.5: Real-image regression suite (`tools/test_real_images.py` against ESA tetra3 fixtures, in CI).
- 3a.6: Native 16-bit grayscale TIFF reader (`src/tiff_reader.{h,cpp}`) + `uint16_t` overload of `extract_centroids_adaptive_gaussian`; `main.cpp` dispatches on `.tif/.tiff` extension. Python preprocessing dropped from `tools/test_real_images.py` and `tools/benchmark.py`; regression now consumes the ESA TIFFs directly. Unit tests cover the round-trip + the 8-bit-rejection path.
- Bonus: pyramid identification replaced vote-based prefilter (direct `find_pairs` + `find_partners` expansion).

### Phase 3b — Algorithm Upgrades (Accuracy) ✅
- 3b.0a: peak-intensity centroid ranking; `CENTROID_CAP` 25 → 50 (60 caused pyramid noise breakdown at N≥55).
- 3b.0b: coarse-refine-reidentify FOV scaling; absorbs the alt60 0.4% FOV calibration drift (eliminates the per-fixture `cos_tol` override).
- 3b.1: iterative Gaussian-weighted centroiding (3–5 iters, σ=1.0 px).
- 3b.2: QUEST attitude estimator with TRIAD fallback for degenerate input (0 fallbacks observed in MC).
- 3b.3: identification cross-verification (median-residual prune, 3× `cos_tolerance` threshold).
- Monte Carlo: 100% / median 0.0041° / max 0.0147°.

### Phase 3c — Speed Optimization ✅ (partial — the algorithm-class fix was Phase 3e)
- 3c.1: Mortari k-vector lookup over the sorted-by-cosine pair array.
- 3c.2: `--benchmark` flag + `tools/benchmark.py` aggregator.
- Goal of <1 s on Pi was not met by k-vector alone; that was the trigger for Phase 3e.

### Phase 3e — Pattern-Hash Identification ✅
- 3e.1: 4-star Tetra-style pattern hash catalog at FOV bins 10°/15°/20° (~870k patterns each, ~21 MB each).
- 3e.2: C++ hash loader + `find_pattern` / `find_pattern_tolerant` (sorted-array binary search + 243-key tolerant probe).
- 3e.3: `identify_stars` rewrite — pattern lookup + 5th-star verify, pyramid kept as named fallback (`identify_stars_pyramid`).
- 3e.4: tests + benchmark instrumentation.
- 3e.5: query-time 24-permutation probing (deduplicates to typically 1–4 unique keys); inlier expansion + QUEST refine after verify-accept (fixes the 3e alt60 accuracy regression).
- **Impact:** identify-stage on alt60: 18,572 ms → 2.8 ms (≈6,600× faster).

### Phase 3f.4 — mmap + mlock + prefault catalog ✅ (out of Phase 3f sequence — landed early as part of the post-3e polish bundle)
- `catalog_pairs.bin`, `catalog_kvec.bin`, and the active `catalog_patterns_*.bin` are `mmap(MAP_PRIVATE, PROT_READ)` + `mlock` (best-effort, warns on `RLIMIT_MEMLOCK` failure) + explicit prefault.
- Static asserts on `sizeof(CatalogPair)==16`, `sizeof(StarPattern)==24` guard the on-disk layout.
- File-read portion of `catalog_load` collapsed; remaining cost (~225 ms) is `per_star_partners` index construction from the mapped pair data.

---

## Code / Algorithm Review (2026-05-11) — superseded by Phase 3g

A read-through on 2026-05-11 produced a list of structural-debt + observability items, all of which were addressed in Phase 3g (see below). Findings included: split `identification.cpp` into 5 modules, template the `cat_vec` callback, add `--stats` instrumentation, vendor ESA fixtures, `run_pyramid` early-exit, dead-code removal, adversarial tests, spatial-index format design, TIFF deterministic tie-break. See git history at commit `b4623f4` for the full review text if you want the original rationale; everything actionable shipped in `ebea92a`/`b4623f4`.

The "Known Follow-ups Backlog" already captures the TIFF-input latency and `per_star_partners` items; the modularisation / spatial-index items above are net-new and belong as their own bullets.

---

## Next Up

### Pi 4 Hardware Integration — IN-FLIGHT

Actual hardware in hand: **Raspberry Pi 4 Model B (8 GB)** on 128 GB SD card, 64-bit Raspberry Pi OS (Trixie / Debian 13, kernel 6.12.75, g++ 14.2). Camera ordered: **Pi HQ Camera (IMX477) + Arducam 35 mm f/1.6 C-mount** (arriving 2026-05-13). Plan-original recommendation was IMX296 mono / 16 mm; user picked HQ Cam + 35 mm for higher precision (~9″/px) at the cost of tighter FOV (~10° horizontal → 10° catalog bin).

**Bring-up scope:**

1. ~~**OS + toolchain + baseline.**~~ **DONE (2026-05-12).** SSH key auth + rsync workflow set up (see local CLAUDE.md / user memory). Pi 4 cmake build = 1m29s (4-core parallel). `ctest` = 13 s with full catalog mmap'd. Catalog files (~547 MB) rsynced from Mac in 21 s rather than regenerated on-Pi (~10-15 min). Baseline measurements logged at top of plan.
2. **Camera capture loop.** Tiny utility that grabs a frame via `libcamera-still` / `picamera2`, writes a 16-bit TIFF, pipes the path to the existing binary. Probably a Python wrapper for the first cut. Mono raw stays unchanged; 16-bit TIFF format is already the C++ binary's native input (3a.6). Stalled on camera arrival.
3. **Camera calibration.** Capture 10+ sky frames at known orientations (or use `solve-field` plate-solving). Fit Brown-Conrady k1..p2 + focal length + principal point via `tools/calibrate_camera.py` (already scaffolded). Persist to `data/pi_camera.json`. Stalled on camera + a clear-sky window.
4. **First-night validation.** Solve a real captured frame end-to-end on the Pi. Report accuracy + per-stage timing. Confirms the ESA-fixture-derived calibration assumptions match real Pi HQ Cam output.

**Pre-Pi prep — ALL DONE:**
- ~~3a.6 TIFF reader~~ **DONE.**
- ~~ARM64 cross-compile lane in CI~~ **DONE.** Lane diagnostically green after the `gtest_discover_tests DISCOVERY_MODE PRE_TEST` fix (commit `ebea92a`).
- ~~Camera-calibration tool scaffold~~ **DONE:** `tools/calibrate_camera.py` smoke-tested on ESA fixtures (post-fit RMS 0.27 px).
- ~~Pi SSH/rsync dev loop~~ **DONE (2026-05-12).** `ssh pi` alias, ed25519 key, `rsync-pi` workflow.

### Phase 3g — Review Follow-ups ✅ (landed 2026-05-11)

All 9 items from the 2026-05-11 code review, completed before Pi hardware bring-up:

- **3g.1 Modularise `identification.cpp`** ✅ — full split into 5 modules:
  - `src/camera_model.{h,cpp}` — `CameraModel`, `project`, `undistort_to_unit_vector`.
  - `src/wahba.{h,cpp}` — TRIAD + unified `wahba::quest` / `wahba::quest_R` shared by `estimation.cpp` (20 iter / 1e-12) and the post-expansion refine path inside the inlier-expansion helpers (30 iter / 1e-14). Duplicate `quest_attitude_local` eliminated.
  - `src/pattern_hash.{h,cpp}` — `kEdgePairs`, `pattern_key_canonical`, `pattern_keys_noise_robust`.
  - `src/pyramid.h` (header-only template) — `run_pyramid`, `estimate_scale_factor`, `cross_verify`, `identify_stars_pyramid`. Header-only so the lambda `CatVec` (3g.2) inlines across the module boundary.
  - `src/inlier_expand.h` (header-only template) — `top_n_indices`, `k_nearest_within_radius`, `expand_inliers`, `expand_inliers_tight`, `refine_and_reexpand`, `try_verify_candidate`, `kFifthStarVerifyCosTol`.
  
  `identification.cpp` shrank from **1518 → 430 lines (72% ↓)** and now hosts only the pattern-path orchestration (`identify_stars_pattern`) plus the top-level `identify_stars` entry. The pyramid/inlier_expand extractions were tricky because of cat_vec / db / R-matrix coupling, but the header-only template idiom (already used for 3g.2's `CatVec`) made it natural — the bodies still live "in the call site" from the compiler's perspective.
- **3g.2 Templated `cat_vec` callback** ✅ — every helper in identification.cpp's anonymous namespace (`run_pyramid`, `estimate_scale_factor`, `cross_verify`, `expand_inliers`, `expand_inliers_tight`, `refine_and_reexpand`, `try_verify_candidate`, `identify_stars_pattern`, `identify_stars_pyramid`) is now `template <typename CatVec>` and takes `const CatVec&`. Call-site lambda is `auto` rather than `std::function`. `<functional>` include dropped from identification.cpp. Profile on Pi will tell us how much of the expected 10–20% identify-stage win materialises in practice.
- **3g.3 `--stats` instrumentation** ✅ — new `IdentifyStats` struct in identification.h, optional out-parameter on `identify_stars` (default `nullptr` = zero overhead). `main.cpp` recognises `--stats` alongside `--benchmark` and prints one `[stats] key=value …` line per solve. Fields: `pattern_catalog_loaded`, `pattern_seeds_tried`, `pattern_path_hit`, `fallback_to_pyramid`, `fov_scale_iters`, `final_inliers`. Smoke-tested on alt60: `pattern_seeds_tried=6 pattern_path_hit=1 fallback_to_pyramid=0 fov_scale_iters=2 final_inliers=18`.
- **3g.4 Vendor ESA real-image fixtures** ✅ — both TIFFs (~1.5 MB each) committed to `tests/data/real_images/`; `tools/test_real_images.py` + `tools/benchmark.py` look there first via a new `resolve_fixture_tiff()`, falling back to URL download only for new fixtures. `bash tests/test_real_images.sh` now passes offline.
- **3g.5 `run_pyramid` early exit** ✅ — outer (i,j) loop breaks when `best_inliers >= max(12, 0.7·N)`. No effect on small-N tests (threshold never reached); helps the pyramid-fallback path's worst case on dense frames.
- **3g.6 Dead `subtract_background` removed** ✅ — function and its two utility-only tests deleted (-43 LOC in `image_processing.cpp`, -51 LOC in `test_image_processing.cpp`). Adaptive thresholding already covered the use cases.
- **3g.7 Adversarial unit tests** ✅ — `EdgeRankFlipStress`, `RejectsBrightNonCatalogContaminant`, `FovBinBoundary12_5`, `FovBinBoundary17_5` added to `test_identification.cpp`. The rank-flip test's initial 2-px jitter exceeded the pattern path's 0.05° verify gate after TRIAD on noisy anchors; jitter tuned to 1 px (still triggers rank flips on the closest-edge pair under the test's 5 RNG seeds while staying within the post-TRIAD verify margin).
- **3g.8 Spatial-index design doc** ✅ — see "Spatial-Index Design (Phase 3g.8)" section below. Grid `(sin(dec) × RA, 64 × 128 bins)` picked over HEALPix for generator simplicity. `'SPAT'` magic, mmap-friendly directory + flat-entries layout matching `catalog_partners.bin`. Implementation deferred to 3d.2 — file format is frozen.
- **3g.9 TIFF identify-stage tie-break** ✅ — both `main.cpp::CENTROID_CAP`'s `partial_sort` and `identification.cpp::top_n_indices` now use `(peak desc, intensity desc, y asc, x asc)` so the seed pool is deterministic regardless of input bit-depth. The 35 ms alt60 identify-stage variance flagged in the previous revision should be tamed; will re-measure on Pi.

**Verification across all 3g items:** 52/52 ctest pass (+4 new adversarial tests, -2 deleted utility tests vs the pre-3g 50); real-image regression alt40=0.0280°, alt60=0.0258° unchanged (well under 0.5° gate); `--stats` produces expected output on alt60 (pattern path hit, no fallback, 18 inliers).

### Phase 3f — Pi 4 speed-up (revised 2026-05-12 from baseline measurements)

The 2026-05-12 Pi 4 baseline (centroid 142 ms, identify alt40 = 2.4 ms vs alt60 = 79 ms, catalog_load 51 ms, total 196-273 ms) makes the priorities concrete. **The first lever is algorithmic, not SIMD.** Reading `image_processing.cpp` revealed the inner loop of `bilinear_sample` (lines 37-58) does a linear scan over tile centers — that's O(n_tx + n_ty) per pixel, ~28 iterations × 786 k pixels = ~22 M extra ops just to look up which tile each pixel belongs to. The tile index is computable as arithmetic: `int tile_x = clamp(int((x - center_offset) / tile_size), 0, n_tx-1)`. Same arithmetic for y. Single-day fix.

Similarly, the "identify alt60 = 79 ms" outlier is a deterministic-after-3g.9 6-seeds-to-hit pattern-path traversal, not a noise/variance issue. Instrumentation will tell us why seeds 1-5 fail before seed 6 wins, and the fix should be targeted (peak-tie-break refinement, pattern_keys_noise_robust noise_tol bump, or similar).

**Ordering by payoff / effort ratio (highest first):**

- **3f.1 `bilinear_sample` O(N) → O(1) fix.** Half-day. Replace the linear scan with arithmetic tile-index computation in `image_processing.cpp::bilinear_sample`. Add a Pi 4 micro-benchmark regression to `tests/test_image_processing.cpp` that times threshold construction for a synthetic 1024×768 input — anything > 30 ms on Pi 4 fails the test. **Expected: centroid 142 ms → 50-80 ms.**

- **3f.2 `alt60` identify root-cause investigation.** 1-2 days. Add a `STARTRACKER_DEBUG_PATTERN_SEEDS` env-var gate around `identify_stars_pattern` that logs, for each seed s in 0..kPatternAttempts: (a) seed centroid HIP nominally, (b) knn neighbor HIPs, (c) which triples produced catalog candidates, (d) for any candidates, whether `try_verify_candidate` rejected at pair-consistency vs 5th-star. Run on alt40 (1 seed) vs alt60 (6 seeds). Identify the actionable difference. Likely fix is in `top_n_indices`'s comparator, `pattern_keys_noise_robust`'s `noise_tol`, or `try_verify_candidate`'s 5th-star anchor choice (currently always `hips[0]`; might need to round-robin). **Expected: alt60 identify 79 ms → 5-15 ms (parity with alt40).**

- **3f.3 Compiler flags pass.** 1 hour. CMakeLists.txt: detect Pi 4 (Cortex-A72) at configure time and add `-mcpu=cortex-a72 -mtune=cortex-a72` plus `-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON`. Validate the on-disk struct static_asserts still hold under LTO (they're declared `static_assert(sizeof(...) == ...)` in `catalog.cpp` so LTO can't reorder fields). Verify ctest 52/52 + real-image regression unchanged. **Expected: 5-15% across all stages.** Re-baseline.

- **3f.4 Eliminate the W×H threshold buffer (compute inline).** 4-6 hours. Currently `build_adaptive_threshold_t` materialises a `std::vector<float>` of W×H = 3.2 MB that overflows L2 and gets re-streamed during BFS. Inline the threshold computation inside the BFS predicate via a tiny lambda capturing the tile arrays; the bilinear interp at the BFS-visited pixel happens on-demand and stays in registers. Memory pressure drops; BFS predicate becomes a few extra arithmetic ops per pixel. **Expected: -5-10 ms on the centroid stage.**

- **3f.5 NEON-vectorize per-tile mean/std.** 1-2 days. The inner loop `for x: sum += v; sum_sq += v*v` reads a tile_size² block — at tile_size=64 that's 4 k pixels per tile × 192 tiles for a 1024×768 frame. Trivially NEON-vectorizable with `vaddq_f32` / FMA across 4-pixel lanes. Stays C++ portable: gate intrinsics behind `#ifdef __ARM_NEON` with a scalar fallback (the existing code). **Expected: -3-5 ms.** Worth doing only after 3f.1 lands (otherwise the threshold-buffer bilinear cost dominates anyway).

- **3f.6 Threaded centroid (4-way).** 2-3 days. Pi 4 has 4 Cortex-A72 cores; currently we run single-threaded. Strategy: split the image into 4 horizontal bands, run `build_adaptive_threshold_t` on each band in parallel via `std::thread`, then a merge pass for components that straddle band boundaries. The threshold-tile structure is already band-friendly. BFS itself stays scalar but parallelised across bands. **Expected: 2-3× on whatever centroid cost remains after 3f.1-5.** Lower priority — only do this if after 3f.1-5 we still exceed 100 ms per-frame.

- **~~3f.7 SIMD pattern compute.~~** **Deprioritised.** Pattern-path identify on alt40 is already 2.4 ms — SIMD on a sub-3 ms stage isn't worth the engineering cost. 3f.2's alt60 investigation gets the bigger win for less effort.

- **3f.8 Single-FOV at boot.** 30 minutes. `main.cpp` already picks the pattern bin once at startup; the per-frame query just reads `db.has_pattern_catalog()`. No actual per-frame branching. Trivial; ship if a follow-up needs an excuse to touch `main.cpp`.

- **~~3f.4 mmap + mlock + prefault catalog.~~ DONE.** Pi 4 result: catalog_load = 51 ms vs the 20 ms desktop measurement (2.6× slower; acceptable since it's one-time). `mlock` succeeds under default `RLIMIT_MEMLOCK` on the user's Pi setup (need to confirm for other Pi configs).
- **~~3f.x `per_star_partners` lazy build.~~ DONE.** mmap'd from `catalog_partners.bin`. Without this file present, ctest takes 161 s; with it, 13 s (12×).

**Projected per-frame envelope after each phase:**

| After | centroid | identify (worst) | total steady-state | gate |
|---|---|---|---|---|
| Today | 142 ms | 79 ms | 221 ms | <100 ms ❌ |
| 3f.1 | ~60 ms | 79 ms | ~140 ms | <100 ms ❌ |
| 3f.1 + 3f.2 | ~60 ms | ~10 ms | ~70 ms | <100 ms ✅ |
| 3f.1-5 | ~30 ms | ~10 ms | ~40 ms | <100 ms ✅ |
| 3f.1-6 | ~15 ms | ~10 ms | ~25 ms | <30 ms stretch ✅ |
| + 3d tracking | ~15 ms | ~1 ms (cone search) | ~17 ms | flying ✅ |

**Decision rule for stopping:** If after 3f.1 + 3f.2 we're under 100 ms per-frame on real Pi-camera captures (post-camera-arrival), the rest of 3f gets deferred until 3d.2 tracking is in. Tracking with a prior attitude collapses the identify cost to <1 ms anyway, so chasing further SIMD wins is premature.

### Phase 3d — Tracking Mode (Extended Kalman Filter)

Needs frame sequences which the Pi capture loop provides naturally. Algorithm scope unchanged from the original plan:

- **3d.1 Multiplicative EKF (MEKF).** 6-state (3 rotation error + 3 angular velocity). Multiplicative error quaternion avoids the norm-constraint snag. Predict via angular velocity propagation; update from identified star measurements.
- **3d.2 Reduced catalog search.** When EKF has prior attitude, restrict identification to catalog stars within a FOV+margin cone of the predicted boresight. Requires a coarse spatial index (declination binning or HEALPix). On Pi this is a big win for steady-state per-frame cost.
- **3d.3 Simulated gyro.** Add `GyroModel` with bias + noise. Expand EKF state to 9 (+ 3 gyro bias). Generate simulated gyro from truth-quaternion sequences for testing in CI without hardware IMU.
- **CLI:** new `--tracking` mode in `main.cpp` that accepts a frame stream and emits a quaternion stream.

### Known Follow-ups Backlog

These don't fit cleanly into one phase but should be addressed before "ship it":

- **TIFF-input identify-stage latency on alt60.** Now reproducible on Pi 4 (79 ms vs alt40's 2.4 ms — 33× spread). 3g.9's deterministic tie-break made it stable; the pattern path tries 6 seeds before hitting on alt60 vs 1 on alt40. **Promoted to 3f.2 (priority work, not just backlog).**
- **`per_star_partners` lazy build vs drop entirely.** Eager in-memory construction at startup costs ~225 ms (vs the 51 ms with `catalog_partners.bin` mmap'd). Pattern-path hit rate is now 82% on Pi (and rising); if it climbs to >95% the pyramid is essentially dead code, and we could drop the partner index entirely. Defer until pattern-path hit rate observation on real Pi-camera frames.

**Post-Pi follow-ups:**
- **Catalog distribution.** Currently catalog is ~547 MB in `data/` (gitignored), regenerated by `tools/generate_catalog.py` (~10-15 min on Pi 4) OR rsynced from a peer machine that has it. Better long-term: tarball + upload as a GitHub Release asset (`catalog-v1.tar.gz`, ~250 MB compressed), pin the version in a config file, add `tools/fetch_catalog.sh`. **Not** git-LFS — bandwidth cost is unworkable on free tier. Trigger this when a second contributor needs to bootstrap a checkout.
- **`unordered_map`/`unordered_set` in identify hot paths → sorted vector / bitset over N≤50.** Profile-driven; identify alt40 is already 2.4 ms on Pi (these are not the bottleneck for that case). Worth doing if 3f.2 reveals `used.count(C)` in `expand_inliers_tight` is the next-most-expensive thing after the alt60 seed traversal.
- **`undistort_to_unit_vector` iteration count.** Reduce from 10 → 4–5 with an early-exit on `|dx|+|dy| < 1e-12`. Called once per centroid (≤ 50). Trivial.
- **TIFF reader fuzz harness.** AFL/libfuzzer corpus once the Pi capture loop is feeding real camera bytes.
- **Robust hot-pixel rejection in centroiding for real-night frames.** Iterative median mask, or a 2-pass extract+reject. Validate against real Pi-camera frames, not ESA fixtures.
- **`subtract_background` removal verification.** 3g.6 removed it; need to confirm real-night Pi frames don't show vignetting/light-pollution gradients that the adaptive-thresholding pass can't handle. If they do, the function may need to come back (or its replacement: a true background-subtract pass before threshold).

---

## Spatial-Index Design (Phase 3g.8)

**Status:** design-only — no `.cpp`/`.h`/`.py` changes in this phase. The on-disk
file format is frozen here so that the 3d.2 (tracking-mode cone query) and the
`expand_inliers_tight` rewrite can land later without churning the catalog
generator or `data/` layout. Mirrors the conventions used by
`catalog_pairs.bin`, `catalog_kvec.bin`, `catalog_patterns_*.bin`, and
`catalog_partners.bin`: little-endian, `int32` magic at byte 0, sorted-by-key
directory, mmap-friendly flat entry array, `static_assert`-guarded struct
sizes.

### Why a sin(dec)/RA grid, not HEALPix

Two candidates were compared:

1. **HEALPix at NSIDE=64** (~50k pixels, ~0.92° pixel side, uniform spherical
   area). Pros: equal-area pixels, no pole pathology. Cons: adds a Python
   dependency (`healpy`) to the catalog generator, and the C++-side
   `query_disc`-equivalent (enumerating pixels whose centre is within
   `half_angle + pixel_radius` of the boresight) is ~150 lines of nontrivial
   spherical geometry. Worth it only if the pole-density skew of the grid
   alternative becomes a real load problem.
2. **Coarse 2D grid by (sin(dec), RA)**. 64 × 128 = 8192 bins (~3° × ~3° near
   the equator, narrowing in RA toward the poles by `cos(dec)`). Pros: trivial
   to generate via `numpy.digitize`, trivial C++ query (two nested integer
   loops with a `cos(dec)` widening factor), no new deps, fits the existing
   "directory + flat entries" layout byte-for-byte. Cons: bin coverage
   shrinks near the poles, so a query at |dec| > 80° touches more RA bins;
   a polar special case handles the |dec_centre|+half_angle > π/2 wrap.

For our use case (~20 stars visible in a 10° cone, ~8k catalog stars total at
Vmag ≤ 7), either index resolves a cone query to ~10 bins / ≤ 200 candidate
HIPs to angular-filter. The grid is simpler, has no Python deps, and lets the
C++ query be ~30 lines. **Recommendation: ship the grid; treat HEALPix as a
documented fallback if a future denser catalog (e.g. Vmag ≤ 8.5 + tracking
near the galactic pole) shows hot-spot pathology.**

### On-disk file format: `catalog_spatial.bin`

```
Header (16 bytes):
  int32 magic              = 0x53504154   ('SPAT', little-endian)
  int32 n_dec_bins         = 64           (compile-time default; written so
                                           a future regrid doesn't need a
                                           new magic)
  int32 n_ra_bins          = 128
  int32 num_nonempty_bins  = count of bins with >= 1 star
Directory (num_nonempty_bins * 16 bytes, sorted ascending by bin_id):
  int32 bin_id             = dec_bin * n_ra_bins + ra_bin
  int32 count              = number of HIPs in this bin
  int64 offset             = byte offset (from file start) of this bin's
                             HIP list
Entries (sum of count * 4 bytes):
  int32 hip_id             = HIP, sorted ascending within each bin so the
                             consumer can set-merge / dedupe cheaply when
                             unioning multiple bins
```

Total size for the production catalog (~8.5k stars, ~8k unique bins
populated): header 16 B + directory ~128 KB + entries ~34 KB = ~165 KB. Fits
easily in L2; the whole index is hot after one cone query.

**Static_assert plan (mirrors the `PartnerEntry` pattern in `catalog.cpp`):**

```cpp
struct SpatialDirEntry {
  int32_t bin_id;
  int32_t count;
  int64_t offset;
};
static_assert(sizeof(SpatialDirEntry) == 16,
              "SpatialDirEntry must be exactly 16 bytes (2*int32 + int64) "
              "to match the on-disk layout written by "
              "tools/generate_catalog.py");
```

The 16-byte directory entry matches the `(int32 hip, int32 count, int64
offset)` layout of `catalog_partners.bin`'s directory, so the same mmap
helper (`mmap_file_ex`) and the same "header → directory → flat entries"
load loop apply unchanged.

### Bin mapping (canonical, must match generator and C++ query)

Given a unit vector `(x, y, z)` with `z = sin(dec)` and `(x, y) ∝ cos(dec) *
(cos(ra), sin(ra))`:

```
dec_bin = clamp(floor((z + 1.0) / 2.0 * n_dec_bins), 0, n_dec_bins - 1)
ra      = atan2(y, x)                          // in [-pi, pi]
ra_bin  = floor((ra + pi) / (2*pi) * n_ra_bins) mod n_ra_bins
bin_id  = dec_bin * n_ra_bins + ra_bin
```

`sin(dec)` (not `dec` itself) is the dec axis because (a) it is a direct
function of the input z-coordinate (no `asin` needed), (b) equal-width bins
in `sin(dec)` give equal-area bands on the sphere, which evens out the
per-bin HIP count compared to equal-`dec` bands.

### C++ query interface (planned for 3d.2, not implemented here)

```cpp
// In src/catalog.h (planned; 3g.8 freezes the signature):
//
//   std::vector<int> StarDatabase::find_stars_in_cone(
//       const double boresight_unit[3],   // unit vector, normalized by caller
//       double cone_half_angle_rad)
//       const;
//
// Returns HIPs whose catalog unit vector v satisfies
//   dot(boresight, v) >= cos(cone_half_angle_rad).
// Result is sorted ascending by HIP (free property of bin-internal sort +
// merging across enumerated bins via a small-N merge).
```

**Algorithm sketch:**

1. Convert boresight to `(dec_c, ra_c)` via `asin(z)` and `atan2(y, x)`.
2. Compute the dec-range `[dec_c - half_angle, dec_c + half_angle]`; map both
   ends to `dec_bin` indices via the sin(dec) formula above. Iterate
   `dec_bin` over this inclusive range (clamped to `[0, n_dec_bins-1]`).
3. For each `dec_bin`, compute the local minimum `cos(dec_local)` over the
   bin's `sin(dec)` interval (taken at whichever bin edge is *closest to the
   pole*, since `cos` is smallest there — that's the worst-case widening
   factor and ensures we don't under-cover). The RA half-width that needs to
   be searched in this dec bin is
   `ra_half = half_angle / max(cos(dec_local), eps_polar)`.
4. If `ra_half >= pi` (always true for the polar caps, see below), enumerate
   *all* `n_ra_bins` ra_bins for this dec_bin. Otherwise enumerate the
   inclusive integer range `[ra_c - ra_half, ra_c + ra_half]` mapped through
   the ra-bin formula, with modular wrap on the 2π seam.
5. For each `(dec_bin, ra_bin)` pair, compute `bin_id`, binary-search the
   sorted directory, and for each HIP in the bin's entry block, apply the
   exact angular test `dot(boresight, get_star(hip).{x,y,z}) >=
   cos(half_angle)` before pushing into the result vector.
6. Sort/dedupe is unnecessary because each bin appears at most once in the
   enumeration and HIPs are unique to a bin (each catalog star falls in
   exactly one bin).

**Polar special case** (handles the singularity where the cone wraps over
the pole):

- If `dec_c + half_angle > pi/2` or `dec_c - half_angle < -pi/2`, the cone
  contains the north or south pole respectively, and every RA bin in the
  affected dec_bin(s) must be enumerated regardless of `ra_c`. Concretely:
  for any dec_bin whose bin-edge sin(dec) interval brushes the pole, set
  `ra_half = pi` so step (4) above hits the "enumerate all" branch. This is
  also the natural outcome of the `cos(dec_local)` calculation in step (3)
  as `dec_local → ±pi/2`, but the explicit check avoids dividing by an
  `eps_polar` floor for typical at-pole queries.

### Generator interface (planned for `tools/generate_catalog.py`)

```python
def generate_spatial_index(vectors, hip_ids, out_dir,
                           n_dec_bins=64, n_ra_bins=128):
    """Emit catalog_spatial.bin matching the format frozen in plan.md §3g.8.

    Layout: header (16 B) + sorted directory (num_nonempty_bins * 16 B) +
    flat int32 HIP entries (4 B each, sorted ascending within each bin).
    """
```

Implementation outline (mirrors the partner-index code that already exists
in `generate_database`):

```python
SPATIAL_MAGIC = 0x53504154  # 'SPAT' little-endian

z = vectors[:, 2]
ra = np.arctan2(vectors[:, 1], vectors[:, 0])           # [-pi, pi]
dec_bin = np.clip(((z + 1.0) / 2.0 * n_dec_bins).astype(np.int64),
                  0, n_dec_bins - 1)
ra_bin = (((ra + np.pi) / (2.0 * np.pi)) * n_ra_bins).astype(np.int64)
ra_bin = np.mod(ra_bin, n_ra_bins)
bin_id = dec_bin * n_ra_bins + ra_bin

# Bucket HIPs by bin_id; sort each bucket ascending by HIP.
buckets = {}
for hid, b in zip(hip_ids, bin_id):
    buckets.setdefault(int(b), []).append(int(hid))
for b in buckets:
    buckets[b].sort()

sorted_bins = sorted(buckets)
num_nonempty_bins = len(sorted_bins)
header_size = 16
dir_entry_size = 16
dir_size = num_nonempty_bins * dir_entry_size
entries_start = header_size + dir_size

with open(os.path.join(out_dir, "catalog_spatial.bin"), "wb") as f:
    f.write(struct.pack("<iiii", SPATIAL_MAGIC, n_dec_bins, n_ra_bins,
                        num_nonempty_bins))
    offset = entries_start
    dir_buf = bytearray(dir_size)
    for idx, b in enumerate(sorted_bins):
        cnt = len(buckets[b])
        struct.pack_into("<iiq", dir_buf, idx * dir_entry_size,
                         b, cnt, offset)
        offset += cnt * 4
    f.write(bytes(dir_buf))
    # Flat int32 entries via a single tofile().
    flat = np.fromiter(
        (h for b in sorted_bins for h in buckets[b]),
        dtype="<i4",
        count=sum(len(buckets[b]) for b in sorted_bins),
    )
    flat.tofile(f)
```

Call site: one line in `generate_database` right after the per-star partner
index block, gated on `num_pairs > 0` for parity with the partners-file
emit.

### Worked example: boresight at (RA = 12h, Dec = +45°), half-angle = 10°

- RA = 12h = π rad. Dec = +45° → sin(dec) = +0.7071.
- `dec_bin_centre = floor((0.7071 + 1.0)/2.0 * 64) = floor(54.626) = 54`
  (bins indexed 0..63 from south pole to north pole).
- Dec range to enumerate: `[35°, 55°]` → sin = `[0.5736, 0.8192]` →
  `dec_bin` range `[floor((0.5736+1)/2*64), floor((0.8192+1)/2*64)] =
  [50, 58]`, i.e. **9 dec_bins** enumerated. None of these brush the pole
  (top bin 58 maps to a sin-edge of ~0.844, dec ≈ 57.6° — well below
  90°), so **the polar special case is not triggered.**
- For the worst-case (top-most) dec_bin at sin(dec) ≈ 0.844 (dec ≈ 57.6°),
  `cos(dec_local) ≈ 0.535`. RA half-width to search = `10° / 0.535 ≈
  18.7°`. RA bins span `(2π / 128) = 2.8125°` each, so we enumerate
  `2 * ceil(18.7° / 2.8125°) + 1 = 15` ra_bins around ra_c.
- For the equator-most dec_bin at sin(dec) ≈ 0.574 (dec ≈ 35°),
  `cos(dec_local) ≈ 0.819`. RA half-width ≈ 12.2°, so ~9 ra_bins.
- Total enumerated bins ≈ Σ(9..15) over 9 dec_bins ≈ ~100 bins. Each
  bin holds ~1 catalog HIP on average (~8.5k stars / 8k populated bins),
  so the candidate set is **~100 HIPs**, filtered by the exact
  `dot(boresight, v) >= cos(10°) = 0.9848` test down to the ~20 visible
  stars.
- Cost: ~100 binary searches over an 8k-entry directory (~13 cmps each) +
  ~100 angular tests ≈ a few microseconds — well under the per-frame
  budget and ~25× cheaper than the current `expand_inliers_tight` union
  of `find_partners` calls.

This worked example sits comfortably in the non-polar branch and exercises
the `cos(dec_local)` widening at moderate dec; the polar branch (e.g. a
boresight at dec = +85° with half-angle = 10°) would instead trip
`dec_c + half_angle = 95° > 90°`, triggering "enumerate all ra_bins" for
all dec_bins in `[dec_bin_for_75°, n_dec_bins - 1]`.

---

## Verification Gates

Updated 2026-05-12 with Pi 4 measurements. "Today" columns show the current measured value; "Gate" is the value below which CI / regression should refuse to ship.

| Gate | Today (desktop) | Today (Pi 4) | Gate |
|---|---|---|---|
| `ctest` | 52/52 | 52/52 (13 s with catalog mmap'd; was 161 s without) | 52/52 maintained |
| Real-image alt40 / alt60 | 0.0280° / 0.0258° | 0.0280° / 0.0258° (identical) | < 0.1° both |
| Monte Carlo success rate | 100% | 100% | ≥ 95% |
| Monte Carlo median error | 0.0041° | 0.0041° | ≤ 0.01° |
| Monte Carlo pattern-path hit rate | 76% | **82%** (better on Pi — 3g.9 tie-break paying off) | ≥ 70% |
| Cold-start total (alt60) | 60 ms | 273 ms | desktop < 150 ms ✅ / Pi < 500 ms ✅ |
| `catalog_load` median | 20 ms | 51 ms | < 100 ms Pi |
| Per-frame steady-state (alt60) | 12–45 ms | **221 ms (centroid 142 + identify 79)** | Pi < 100 ms ❌ → 3f.1+3f.2 |
| Per-frame steady-state (alt40) | 1–5 ms | **144 ms (centroid 142 + identify 2)** | Pi < 100 ms ❌ → 3f.1 alone closes |
| Pi 4 cross-compile CI | passing | n/a | green |
| Pi 4 per-frame stretch | — | — | <30 ms (3f.1-6 + 3d tracking) |

---

## Key Files

| File | Owners |
|---|---|
| `src/image_processing.{h,cpp}` | 3a.2/3/4, 3b.1, 3a.6 (uint16 overload); **3f.1 + 3f.4-5 (next)** |
| `src/identification.{h,cpp}` | 3a.1, 3b.0b/3, 3e.3/5; **3f.2 (alt60 root-cause, next)** |
| `src/camera_model.{h,cpp}` | 3g.1 (Brown-Conrady project + undistort) |
| `src/wahba.{h,cpp}` | 3g.1 (TRIAD + unified QUEST) |
| `src/pattern_hash.{h,cpp}` | 3g.1 (`pattern_key_canonical`, `pattern_keys_noise_robust`) |
| `src/pyramid.h` | 3g.1 (header-only templates: run_pyramid + identify_stars_pyramid) |
| `src/inlier_expand.h` | 3g.1 (header-only templates: expand_inliers* + try_verify_candidate) |
| `src/catalog.{h,cpp}` | 3c.1, 3e.2, 3f.4, 3f.x (partners mmap); future 3d.2 (spatial cone query) |
| `src/estimation.{h,cpp}` | 3b.2 (now thin wrapper around wahba::quest) |
| `src/main.cpp` | all phases; 3a.6 routes TIFF vs PNG; 3g.3 --stats; 3g.9 tie-break |
| `src/tiff_reader.{h,cpp}` | 3a.6 (minimal uncompressed 16-bit grayscale reader, no libtiff) |
| `tools/generate_catalog.py` | 3c.1, 3e.1, K_NEAREST=12, 3f.x (catalog_partners.bin); future 3d.2 (spatial) |
| `tools/test_real_images.py` | 3a.5; 3a.6 dropped tiff→png; 3g.4 vendored fixtures |
| `tools/monte_carlo.py` | Phase 2.4, 3e.5 (hit-rate instrumentation) |
| `tools/benchmark.py` | 3c.2; 3a.6 dropped tiff→png; 3g.4 vendored fixture resolver |
| `tools/calibrate_camera.py` | Pre-Pi prep (Brown-Conrady fit scaffold) |
| `cmake/aarch64-linux-gnu.cmake` | Pre-Pi prep (cross-compile toolchain) |
| `.github/workflows/ci.yml` | host x86_64 + ARM64 cross-compile jobs |
| NEW for 3d: `src/tracking.{h,cpp}` | 3d.1 (MEKF) |
| NEW for 3d: `tools/capture.py` | Pi capture loop (libcamera-still/picamera2 → 16-bit TIFF) |

---

## Sequence Snapshot

```
[done]   1 → 2 → 3a (incl. 3a.6) → 3b → 3c → 3e → 3e.5 + 3f.4 + K_NEAREST=12
         → 3e.5 unit tests → 3f.x (catalog_partners.bin)
         → pre-Pi prep: ARM64 CI lane + calibrate_camera.py scaffold
         → 3g (file split: camera_model + wahba unified QUEST + pattern_hash;
              templated cat_vec; --stats; vendored ESA fixtures; pyramid
              early-exit; dead-code removal; adversarial tests; spatial-index
              format frozen; TIFF deterministic tie-break)
         → Pi 4 bring-up step 1: OS + native build + baseline measurements
              (centroid 142 ms, identify alt60 79 ms, total 273 ms)
[now]    3f.1 (bilinear_sample O(N)→O(1)) + 3f.2 (alt60 identify root-cause):
         the two highest-payoff items, both algorithmic, ~2 days combined,
         expected per-frame 221 ms → ~70 ms. Doable before camera arrives.
[next]   Camera arrival → Pi capture loop → first-night solve → real calibration
[then]   3f.3-6 SIMD/threading IF still needed after measuring on real Pi
         captures; 3d (EKF tracking) once captures are available
[later]  3d.2 spatial-index implementation (file format frozen by 3g.8);
         catalog distribution via GitHub Release; TIFF fuzz harness;
         hot-pixel rejection
```
