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

- **TIFF-input identify-stage latency.** Native 16-bit TIFF input regresses alt60 identify from ~5 ms to ~35 ms vs the previous PNG-stretched path. Same final attitude. Root cause is centroid peak-ordering shifts (8-bit-saturated peaks collapse to ties that broke one way; 16-bit clean peaks order differently and route the pattern path through more failed seeds). Re-measure on Pi before tuning — may be moot once we have real Pi-camera data instead of ESA fixtures.
- **`per_star_partners` lazy build.** Eager construction at startup costs ~225 ms and ~250 MB. After 3e.5 + K_NEAREST=12, the pyramid is fallback-only for 75-90% of MC traffic. Lazy-build (or drop entirely if pattern hit rate climbs further) is the next biggest cold-start win.

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
[next]   Pi 4 hardware bring-up (OS, native build, first sky capture)
[then]   3f (remaining SIMD/cache work, measured on Pi)
[later]  3d (EKF tracking — needs frame stream from Pi capture loop)
```
