# Project Plan: *startracker*

Iterative development of an embedded star-tracking engine: synthetic-data baseline → CI safety net → real-world robustness → accuracy upgrades → speed upgrades → Pi hardware integration → SIMD + tracking. Each phase is regression-tested against the prior phase's gates.

---

## Where We Are (2026-05-11)

**Current capability:** cold-start solve on alt60 ESA tetra3 fixture (1024×768, ~11° FOV, 0.4% FOV calibration drift) in **~250 ms** end-to-end on M-series desktop. Per-stage median:

| Stage | Time |
|---|---|
| centroid | 13 ms |
| catalog_load | 236 ms (~225 ms is `per_star_partners` build, not file I/O) |
| identify (pattern hash + 24-perm probe + QUEST refine) | 2.8 ms |
| estimate | <1 µs |

Per-frame steady-state (catalog already loaded): **~16 ms**. Memory: 570 MB RSS.

**Test coverage:** 44/44 unit tests; real-image regression alt40=0.0596°, alt60=0.0000° (gate 0.5°); Monte Carlo 50 trials @ 5″ noise = 100% success, median 0.0047°, max 0.0193° (gate ≥80% / <1°).

**Known follow-ups carried into future phases:**
- Pattern-path hit rate on noisy synthetic Monte Carlo is **20%**; pyramid carries the rest. Root cause: synthetic Vmag~7.5 stars exceed catalog Vmag-7 cutoff, so the 8-brightest set includes non-catalog stars. Real-image fixtures hit pattern path 100%. Fix: either tighten centroid pre-filter or extend catalog to Vmag 7.5.
- alt40 went from 0.0000° → 0.0596° during 3e.5 (inlier expansion now grabs more stars and QUEST averages over them; the "extra" stars on alt40 include some marginal matches). Still well under the 0.5° gate; investigate if it becomes a problem.
- `per_star_partners` index is built eagerly at startup (~225 ms, ~250 MB). Lazy build on first pyramid call would cut both cold-start and idle RAM substantially, since pattern-path covers real-image traffic. Touches the Phase 3f.4 boundary.

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
- Bonus: pyramid identification replaced vote-based prefilter (direct `find_pairs` + `find_partners` expansion).
- **Deferred 3a.6 (native 16-bit TIFF):** still deferred to the Pi hardware integration window.

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
- 3a.6 TIFF reader (portable code; ESA TIFFs are regression data).
- ARM64 cross-compile lane in CI (catches alignment / endianness / missing-intrinsics bugs before flashing an SD card).
- Camera-calibration tool scaffold (the math is hardware-independent; just needs sample images).

### Phase 3f — SIMD & Cache Engineering (post-Pi measurement)

Only worth doing if the on-Pi measurement from the previous phase shows we're slower than the application needs. Target if pursued: **<30 ms total on Pi 4, <3 ms on desktop.**

- **3f.1 SIMD centroiding.** NEON (Pi/M-series) + AVX2 (x86) intrinsics in `subtract_background` (tile median + bilinear interp + clamp) and `extract_centroids_adaptive_gaussian` (per-tile mean/stddev precompute; BFS itself is scalar but the above-threshold predicate vectorizes). Expected: 3–5× on the centroid stage.
- **3f.2 SIMD pattern compute.** Vectorize the 6 inter-star angles + sort + quantize into a single AVX2/NEON lane. Small absolute compute but runs 24× per query post-3e.5.
- **3f.3 Robin Hood hash for the pattern table.** Currently a sorted vector + binary search — single cache miss per probe is roughly equivalent in practice. Only revisit if profiling on Pi shows the binary search is bottlenecked.
- **~~3f.4 mmap + mlock + prefault catalog.~~** **DONE** (out of sequence — landed in the post-3e polish bundle). Need to revisit on Pi for `mlock` behavior under default `RLIMIT_MEMLOCK`.
- **3f.5 Single-FOV assumption.** Pre-pick the pattern catalog bin at boot rather than per-frame. Trivial; mostly bookkeeping.
- **Phase 3f.x (new candidate): defer `per_star_partners` construction.** Currently builds eagerly at startup (~225 ms, ~250 MB). After 3e.5 the pyramid is fallback-only; lazy-build (or drop entirely if the synthetic-noise pattern miss rate gets fixed) is the next biggest cold-start win.

### Phase 3d — Tracking Mode (Extended Kalman Filter)

Needs frame sequences which the Pi capture loop provides naturally. Algorithm scope unchanged from the original plan:

- **3d.1 Multiplicative EKF (MEKF).** 6-state (3 rotation error + 3 angular velocity). Multiplicative error quaternion avoids the norm-constraint snag. Predict via angular velocity propagation; update from identified star measurements.
- **3d.2 Reduced catalog search.** When EKF has prior attitude, restrict identification to catalog stars within a FOV+margin cone of the predicted boresight. Requires a coarse spatial index (declination binning or HEALPix). On Pi this is a big win for steady-state per-frame cost.
- **3d.3 Simulated gyro.** Add `GyroModel` with bias + noise. Expand EKF state to 9 (+ 3 gyro bias). Generate simulated gyro from truth-quaternion sequences for testing in CI without hardware IMU.
- **CLI:** new `--tracking` mode in `main.cpp` that accepts a frame stream and emits a quaternion stream.

### Known Follow-ups Backlog

These don't fit cleanly into one phase but should be addressed before "ship it":

- **Pattern-path noise robustness on Monte Carlo.** Currently 20% hit rate; pyramid carries the rest. Either (a) tighten the centroid pre-filter so we only pass mag-7-and-brighter, or (b) regenerate the pattern catalog with `max_mag=7.5` so the noisy 8-brightest set has full coverage. Option (b) costs ~2× catalog size; option (a) costs noise robustness on real frames where mag-7-and-brighter is the *upper* bound. Decide after Pi data.
- **alt40 accuracy.** Went from 0.0000° → 0.0596° during 3e.5. Investigate which marginal inliers QUEST is averaging.
- **3e.5 unit tests not landed.** Empirical validation (real-image + MC) covers the regression but the three tests requested in the 3e.5 spec (`PermutationProbeFindsKey`, `NoiseRobustnessSweep`, `AccuracyAfterVerify`) were skipped by the subagent. Add them when convenient.

---

## Verification Gates

| Gate | Current | Pre-Pi target |
|---|---|---|
| `ctest` | 44/44 | maintained |
| Real-image alt40 / alt60 | 0.0596° / 0.0000° | both < 0.1° after alt40 investigation |
| Monte Carlo success rate | 100% | maintained ≥ 95% |
| Monte Carlo median error | 0.0047° | maintained ≤ 0.01° |
| Cold-start total (alt60, desktop) | 250 ms | target after 3a.6 + cross-compile CI: still 250 ms (no change expected) |
| Per-frame steady-state (alt60, desktop) | 16 ms | maintained |
| Pi 4 cold-start total | unmeasured | <500 ms (initial gate; refine after measuring) |
| Pi 4 per-frame steady-state | unmeasured | <100 ms (initial gate; <30 ms is the 3f stretch goal) |

---

## Key Files

| File | Owners |
|---|---|
| `src/image_processing.{h,cpp}` | 3a.2/3/4, 3b.1, 3a.6 (TIFF), future 3f.1 |
| `src/identification.{h,cpp}` | 3a.1, 3b.0b/3, 3e.3/5, future 3f.x (per_star_partners defer) |
| `src/catalog.{h,cpp}` | 3c.1, 3e.2, 3f.4, future 3d.2 (spatial index) |
| `src/estimation.{h,cpp}` | 3b.2 |
| `src/main.cpp` | all phases |
| `tools/generate_catalog.py` | 3c.1, 3e.1, possibly future re-extension to Vmag 7.5 |
| `tools/test_real_images.py` | 3a.5, post-3a.6 simplification |
| `tools/monte_carlo.py` | Phase 2.4, 3e.5 (hit-rate instrumentation) |
| `tools/benchmark.py` | 3c.2 |
| NEW: `src/tiff_reader.{h,cpp}` | 3a.6 (during Pi window) |
| NEW: `src/tracking.{h,cpp}` | 3d.1 |
| NEW: `tools/calibrate_camera.py` | Pi window |

---

## Sequence Snapshot

```
[done]   1 → 2 → 3a → 3b → 3c → 3e → 3e.5 + 3f.4 polish
[next]   Pi 4 bring-up + 3a.6 (16-bit TIFF) + camera calibration
[then]   3f (SIMD + remaining cache work, measured on Pi)
[later]  3d (EKF tracking — needs frame stream from Pi capture loop)
```
