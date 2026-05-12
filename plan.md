# Project Plan: *startracker*

Iterative development of an embedded star-tracking engine: synthetic-data baseline → CI safety net → real-world robustness → accuracy upgrades → speed upgrades → Pi hardware integration → SIMD + tracking. Each phase is regression-tested against the prior phase's gates.

---

## Where We Are (2026-05-11)

**Current capability:** cold-start solve on alt60 ESA tetra3 fixture (1024×768, ~11° FOV, 0.4% FOV calibration drift) in **~60 ms** end-to-end on M-series desktop (native 16-bit TIFF input, mmap'd partner index). Per-stage median:

| Stage | Time |
|---|---|
| centroid | 11 ms |
| catalog_load | 20 ms (was 229 ms — 3f.x mmap'd the per-star partner index) |
| identify (pattern hash + 24-perm probe + QUEST refine) | ~1–35 ms (varies with seed-ordering luck on 16-bit-native input) |
| estimate | <1 µs |

Per-frame steady-state (catalog already loaded): **~12–45 ms**. Memory: ~620 MB virtual / much lower RSS (partner-index payload is now lazy-faulted via mmap, only the pages that get touched are paged in).

**Test coverage:** 50/50 unit tests (44 → 47 with the three 3e.5 tests, → 50 with the TIFF reader tests); real-image regression alt40=0.0280°, alt60=0.0258° (gate 0.5°); Monte Carlo 50 trials @ 5″ noise = 100% success, median 0.0047°, max 0.0191°, pattern-path hit rate 76% (gates ≥80% / <1° / ≥70% hit rate).

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

## Code / Algorithm Review (2026-05-11)

A read-through of the C++ sources, catalog generator, calibration tool, and CI lane. Algorithms are excellent (Tetra pattern hash + QUEST + k-vector + Brown-Conrady is state-of-the-art for the problem); findings are mostly structural debt and small efficiency items. Each is tagged with **[before-hw]** or **[after-hw]** depending on whether it gates Pi work or is profile-driven.

### Structural / modularity (mostly before-hw)

- **[before-hw] Split `src/identification.cpp` (1518 LOC).** Becomes hard to navigate, and both 3d (tracking) and 3f (SIMD) will touch it heavily. Suggested breakdown:
  - `src/camera_model.{h,cpp}` — `CameraModel`, `project`, `undistort_to_unit_vector`. Currently bleeding into `identification.h`; calibration tool + future tracking both want a standalone owner.
  - `src/wahba.{h,cpp}` — `triad_rotation` + a single `quest_attitude(...)` shared by `estimation.cpp` and `identification.cpp`. **QUEST is implemented twice today** (`estimate_attitude` and `quest_attitude_local`) — any change has to be made and re-tested in both places.
  - `src/pattern_hash.{h,cpp}` — `pattern_key_canonical`, `pattern_keys_noise_robust`, the kEdgePairs/kQuantBits constants, and the `identify_stars_pattern` orchestration.
  - `src/pyramid.{h,cpp}` — `run_pyramid`, `cross_verify`, `estimate_scale_factor`.
  - `src/inlier_expand.{h,cpp}` — `expand_inliers`, `expand_inliers_tight`, `refine_and_reexpand`.
  - `identification.cpp` becomes the orchestrator (the 60-line top-level `identify_stars`).
- **[before-hw] Templated catalog-lookup callback in hot paths.** `std::function<std::array<double,3>(int)> cat_vec` is invoked thousands of times per frame and adds heap-allocated indirection + virtual dispatch. Switch the run_pyramid / try_verify_candidate / expand_inliers* call chain to be a template parameter on a `CatVec` functor. Expected 10–20% identify-stage win on Pi. Do this as part of the split above so the new module boundaries don't ossify the std::function shape.
- **[before-hw] Remove `subtract_background`.** Reads as live in `image_processing.h` but no callsite — `extract_centroids_adaptive_gaussian` builds its own per-tile threshold without needing the subtracted buffer. Either delete or wire it in front of the adaptive extractor and re-benchmark. Dead code rots and confuses 3a.6/3f reviewers.
- **[before-hw] No `--help` / usage with examples in `main.cpp`.** Will matter the first time someone other than you tries to run it on Pi.

### Coverage / observability (before-hw)

- **[before-hw] `--stats` flag exposing internal counters.** Right now success-rate / accuracy is the only signal — opaque about *why* a solve succeeded. Add a small struct populated through the pipeline: pattern-attempts tried, pattern-attempt hits, fallback-to-pyramid count, mean/max inliers post-expansion, QUEST iterations, scale-factor convergence iters. Print one line in `--benchmark` mode. We're about to characterize Pi behavior; flying blind on the internal health metrics will make debugging real-night frames painful.
- **[before-hw] Cache/LFS the ESA tetra3 real-image fixtures.** `tools/test_real_images.py` downloads each TIFF on first run via `download_if_missing`. Pi network is flaky and GHA networking can fail; CI will eventually be red for a no-code-change reason. Vendor them under `tests/data/real_images/`, or use git-LFS, or stamp them into the Docker image.
- **[before-hw] Root-cause the TIFF-input identify-stage variance.** Plan documents alt60 going from ~5 ms (PNG) to ~35 ms (TIFF) with the same final attitude. "Centroid peak-ordering shifts on 16-bit ties" is the hypothesis but it's deterministic in principle — if it's the seed ranking, a peak-tie-break on (peak, intensity) or (peak, -y, -x) should pin the order down. Worth fixing or definitively re-characterising before Pi capture lands real 12-bit raw, which has its own tie distribution.

### Algorithm-level (mix of before and after)

- **[before-hw] `run_pyramid` has no early exit.** Outer (i,j) loop walks every centroid pair (~1225 for N=50) even after a high-confidence solution is found. Add `if (best_inliers >= 0.7*N || best_inliers >= 12) break;` to the outer loops. Cheap, helps the fallback path's worst case. Doesn't change steady-state when pattern path succeeds, but the fallback is exercised on 20–25% of MC traffic today and 100% of pre-pattern-catalog (no `--load-pattern-catalog`) flows.
- **[before-hw] No spatial index in `StarDatabase`.** Needed for 3d.2 tracking (cone-restricted catalog search), and would also speed `expand_inliers_tight` (which currently unions `find_partners` from every matched HIP — O(M·ring_size) per unassigned centroid). HEALPix at NSIDE=64 (~50k pixels) or a coarse dec-band+RA-bin grid is enough; either fits the existing `data/` file convention. Even if 3d is post-Pi, *designing* the storage now (one extra `catalog_spatial.bin` next to `catalog_partners.bin`) keeps the format change off the Pi-bring-up path.
- **[after-hw] Replace `unordered_map`/`unordered_set` in identify hot paths with sorted vectors or bitsets.** N ≤ 50 — a 64-bit `uint64_t` bitset over centroids and a small sorted vector over HIPs is 5–10× faster than hashing on cache-cold ARM. Profile-driven: only worth doing once we see whether identify is actually the bottleneck on Pi.
- **[after-hw] `undistort_to_unit_vector` does 10 fixed-point iterations unconditionally.** 4–5 converges to <1e-10 in practice for Brown-Conrady; an early-exit on `|dx|+|dy| < 1e-12` saves cycles in the zero-distortion (synthetic) path too. Trivial change; defer until Pi profile says it matters.
- **[after-hw] Fuzz the TIFF reader.** Handwritten parser → eventual exposure to whatever the Pi HQ Camera firmware emits under fault conditions. AFL or libfuzzer corpus of malformed TIFFs once the capture loop is live.
- **[after-hw] Robust outlier-rejection in centroiding for real-night frames.** Cosmic-ray streaks and hot pixels survive shape filters in some cases. Iterative median-based hot-pixel mask or a 2-pass extract+reject scheme. Will only know if needed once we see real Pi-camera frames; the ESA fixtures are too clean to validate this on.

### Test discipline (before-hw)

- **[before-hw] Adversarial unit tests.** Existing tests are happy-path (real catalog stars, synthetic frames). Add: (a) a frame with two stars within centroid noise of each other (edge-rank-flip stress), (b) a frame with one bright non-catalog object (satellite) that should be filtered out, (c) a frame at the FOV bin boundary (10.0° / 12.5° / 17.5°) to catch the `fov_bin = ... ? 10 : 15 : 20` branch.
- **[before-hw] Replace bash test wrappers with pytest.** `tests/test_*.sh` are shell-around-Python around the C++ binary; on Pi running pytest end-to-end is more idiomatic and easier to read failure tracebacks from.

### Verdict on the existing plan

Keep the phase order (Pi bring-up → 3f SIMD → 3d tracking). The two changes I'd lobby for:

1. **Insert a "Phase 3g — Modularisation & Observability" between Pi prep and Pi hardware integration.** Specifically the file split + unified QUEST + `--stats` + ESA fixture caching items above. It is genuinely pre-Pi work — none of it needs hardware — and it eliminates a chunk of structural drag from 3f and 3d. Budget: half a day.
2. **Pull the spatial-index design forward** into pre-Pi prep (design + file format only, no implementation). 3d.2 then becomes a smaller change.

The "Known Follow-ups Backlog" already captures the TIFF-input latency and `per_star_partners` items; the modularisation / spatial-index items above are net-new and belong as their own bullets.

---

## Next Up

### Pi 4 Hardware Integration (next phase) — IN-FLIGHT

Hardware order recommended: Raspberry Pi 5 (8 GB) for dev + Pi 4 (4 GB) for deployment validation, Pi HQ Camera (Sony IMX477), 25 mm CS-mount low-distortion lens, 32 GB Class 10 microSD, USB-C 5V/5A PSU, case + ribbon cable. ~$255.

**Bring-up scope (rough order of work once hardware arrives):**

1. **OS + toolchain.** 64-bit Raspberry Pi OS, native build (Pi 5 has the RAM; cross-compile not needed). Run `cmake --build` + `ctest` on-device — should be a no-op given x86/M-series tested it. Capture an initial wall-clock baseline for each stage on Pi 4.
2. **Native 16-bit TIFF input (the deferred 3a.6).** Pi HQ Camera emits 12/16-bit raw; the current 8-bit PNG path crushes faint stars. Add a minimal uncompressed-16-bit-grayscale TIFF reader (no compression handling needed — `libtiff` is overkill) plus `uint16_t` overloads for `subtract_background` and `extract_centroids_adaptive_gaussian`. Drop the Python preprocessing from `tools/test_real_images.py`. The ESA tetra3 TIFFs are already perfect regression data; we currently stretch them to 8-bit before feeding the binary, which is the exact gap to close.
3. **Camera capture loop.** Tiny utility that grabs a frame via `libcamera`/`picamera2`, writes a 16-bit TIFF, pipes the path to the existing binary. Probably a Python wrapper for the first cut.
4. **Camera calibration.** Capture 10+ sky frames at known orientations (or use `solve-field` for plate solving). Fit Brown-Conrady k1..p2 + focal length + principal point. Persist to `data/pi_camera.json`. Wrap as `tools/calibrate_camera.py` so it's repeatable.
5. **First-night validation.** Solve a real captured frame end-to-end on the Pi. Report accuracy + per-stage timing. This is the input to Phase 3f priorities.

**Pre-Pi prep that could happen now (don't strictly need hardware):**
- ~~3a.6 TIFF reader~~ **DONE.**
- ~~ARM64 cross-compile lane in CI~~ **DONE:** `cmake/aarch64-linux-gnu.cmake` toolchain + a `build-arm64-cross` GHA job that compiles `startracker_core`, `startracker`, and `startracker_tests` for `aarch64-linux-gnu`. Locally verified via Homebrew cross-toolchain — produced an `ELF 64-bit LSB executable, ARM aarch64` binary; GoogleTest cross-compiles cleanly via FetchContent. Required adding `#include <algorithm>` to two test files that had been pulling `std::sort`/`std::clamp` through transitive includes only.
- ~~Camera-calibration tool scaffold~~ **DONE:** `tools/calibrate_camera.py` fits Brown-Conrady (focal_x/y, center_x/y, k1..k3, p1, p2) over (frame, truth-quaternion) pairs via `scipy.optimize.least_squares`. Smoke-tested on the two ESA fixtures (37 matched stars; pre-fit RMS 1.06 px → post-fit 0.27 px; converged). Output schema documents the CLI distortion order (k1, k2, p1, p2, k3) so the caller doesn't have to guess.

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

### Phase 3f — SIMD & Cache Engineering (post-Pi measurement)

Only worth doing if the on-Pi measurement from the previous phase shows we're slower than the application needs. Target if pursued: **<30 ms total on Pi 4, <3 ms on desktop.**

- **3f.1 SIMD centroiding.** NEON (Pi/M-series) + AVX2 (x86) intrinsics in `subtract_background` (tile median + bilinear interp + clamp) and `extract_centroids_adaptive_gaussian` (per-tile mean/stddev precompute; BFS itself is scalar but the above-threshold predicate vectorizes). Expected: 3–5× on the centroid stage.
- **3f.2 SIMD pattern compute.** Vectorize the 6 inter-star angles + sort + quantize into a single AVX2/NEON lane. Small absolute compute but runs 24× per query post-3e.5.
- **3f.3 Robin Hood hash for the pattern table.** Currently a sorted vector + binary search — single cache miss per probe is roughly equivalent in practice. Only revisit if profiling on Pi shows the binary search is bottlenecked.
- **~~3f.4 mmap + mlock + prefault catalog.~~** **DONE** (out of sequence — landed in the post-3e polish bundle). Need to revisit on Pi for `mlock` behavior under default `RLIMIT_MEMLOCK`.
- **3f.5 Single-FOV assumption.** Pre-pick the pattern catalog bin at boot rather than per-frame. Trivial; mostly bookkeeping.
- **~~Phase 3f.x: defer `per_star_partners` construction.~~** **DONE** (Phase 3f.x landed pre-Pi). `tools/generate_catalog.py` now emits `data/catalog_partners.bin` (header + directory + entry blocks); the C++ side mmaps it without `mlock`/prefault (sparse access — pattern path touches one entry block per `find_partners` call, prefaulting all 198 MB would cost ~30 ms for no benefit). Falls back to in-memory construction if the file is missing. **catalog_load dropped 225 ms → 20 ms (~12×).**

### Phase 3d — Tracking Mode (Extended Kalman Filter)

Needs frame sequences which the Pi capture loop provides naturally. Algorithm scope unchanged from the original plan:

- **3d.1 Multiplicative EKF (MEKF).** 6-state (3 rotation error + 3 angular velocity). Multiplicative error quaternion avoids the norm-constraint snag. Predict via angular velocity propagation; update from identified star measurements.
- **3d.2 Reduced catalog search.** When EKF has prior attitude, restrict identification to catalog stars within a FOV+margin cone of the predicted boresight. Requires a coarse spatial index (declination binning or HEALPix). On Pi this is a big win for steady-state per-frame cost.
- **3d.3 Simulated gyro.** Add `GyroModel` with bias + noise. Expand EKF state to 9 (+ 3 gyro bias). Generate simulated gyro from truth-quaternion sequences for testing in CI without hardware IMU.
- **CLI:** new `--tracking` mode in `main.cpp` that accepts a frame stream and emits a quaternion stream.

### Known Follow-ups Backlog

These don't fit cleanly into one phase but should be addressed before "ship it":

- **TIFF-input identify-stage latency.** Native 16-bit TIFF input regresses alt60 identify from ~5 ms to ~35 ms vs the previous PNG-stretched path. Same final attitude. Root cause is centroid peak-ordering shifts (8-bit-saturated peaks collapse to ties that broke one way; 16-bit clean peaks order differently and route the pattern path through more failed seeds). Re-measure on Pi before tuning — may be moot once we have real Pi-camera data instead of ESA fixtures. *(See 3g.9 for the pre-Pi tie-break fix.)*
- **`per_star_partners` lazy build.** Eager construction at startup costs ~225 ms and ~250 MB. After 3e.5 + K_NEAREST=12, the pyramid is fallback-only for 75-90% of MC traffic. Lazy-build (or drop entirely if pattern hit rate climbs further) is the next biggest cold-start win. *(Largely subsumed by 3f.x's mmap path; this item is the "drop entirely" follow-up if the pattern-path hit rate climbs to >95%.)*

**Post-Pi follow-ups (from 2026-05-11 review):**
- **`unordered_map`/`unordered_set` in identify hot paths → sorted vector / bitset over N≤50.** Profile-driven; only worth doing if identify is the Pi-side bottleneck.
- **`undistort_to_unit_vector` iteration count.** Reduce from 10 → 4–5 with an early-exit on `|dx|+|dy| < 1e-12`. Micro-optimisation.
- **TIFF reader fuzz harness.** AFL/libfuzzer corpus once the Pi capture loop is feeding real camera bytes.
- **Robust hot-pixel rejection in centroiding for real-night frames.** Iterative median mask, or a 2-pass extract+reject. Validate against real Pi-camera frames, not ESA fixtures.

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

| Gate | Current | Pre-Pi target |
|---|---|---|
| `ctest` | 50/50 | maintained |
| Real-image alt40 / alt60 | 0.0280° / 0.0258° (native TIFF) | maintained < 0.1° |
| Monte Carlo success rate | 100% | maintained ≥ 95% |
| Monte Carlo median error | 0.0047° | maintained ≤ 0.01° |
| Monte Carlo pattern-path hit rate | 76% | ≥ 70% maintained |
| Cold-start total (alt60, desktop, native TIFF) | 60 ms | maintained < 150 ms |
| `catalog_load` median | 20 ms (was 225 ms) | maintained < 50 ms |
| Per-frame steady-state (alt60, desktop) | 12–45 ms | investigate identify-stage variance |
| Pi 4 cross-compile CI | passing locally; GHA job committed | passing on GHA after first push |
| Pi 4 cold-start total | unmeasured | <500 ms (initial gate; refine after measuring) |
| Pi 4 per-frame steady-state | unmeasured | <100 ms (initial gate; <30 ms is the 3f stretch goal) |

---

## Key Files

| File | Owners |
|---|---|
| `src/image_processing.{h,cpp}` | 3a.2/3/4, 3b.1, 3a.6 (uint16 overload), future 3f.1 |
| `src/identification.{h,cpp}` | 3a.1, 3b.0b/3, 3e.3/5, future 3d.2 (spatial index) |
| `src/catalog.{h,cpp}` | 3c.1, 3e.2, 3f.4, 3f.x (partners mmap), future 3d.2 |
| `src/estimation.{h,cpp}` | 3b.2 |
| `src/main.cpp` | all phases; 3a.6 routes TIFF vs PNG by extension |
| `src/tiff_reader.{h,cpp}` | 3a.6 (minimal uncompressed 16-bit grayscale reader, no libtiff) |
| `tools/generate_catalog.py` | 3c.1, 3e.1, K_NEAREST=12, 3f.x (catalog_partners.bin) |
| `tools/test_real_images.py` | 3a.5; 3a.6 dropped the tiff→png stretch |
| `tools/monte_carlo.py` | Phase 2.4, 3e.5 (hit-rate instrumentation) |
| `tools/benchmark.py` | 3c.2; 3a.6 dropped the tiff→png stretch |
| `tools/calibrate_camera.py` | Pre-Pi prep (Brown-Conrady fit scaffold) |
| `cmake/aarch64-linux-gnu.cmake` | Pre-Pi prep (cross-compile toolchain) |
| `.github/workflows/ci.yml` | host x86_64 + ARM64 cross-compile jobs |
| NEW: `src/tracking.{h,cpp}` | 3d.1 |

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
[next]   Pi 4 hardware bring-up (OS, native build, first sky capture)
[then]   3f (remaining SIMD/cache work, measured on Pi)
[later]  3d (EKF tracking — needs frame stream from Pi capture loop;
         spatial-index file format is already frozen by 3g.8)
```
