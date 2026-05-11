# startracker

An embedded star tracking engine for spacecraft attitude determination. Takes a star field image and outputs a quaternion representing the spacecraft's orientation.

## How It Works

1. **Image processing** — extracts star centroids from a grayscale PNG using connected-component BFS and center-of-gravity centroiding
2. **Star identification** — matches centroids to the Hipparcos catalog using pairwise angular-distance voting
3. **Attitude estimation** — computes a quaternion using the TRIAD algorithm

## Build

```bash
mkdir build && cd build && cmake .. && make -j$(nproc)
```

## Usage

```bash
./build/startracker <image.png> <catalog_stars.bin> <catalog_pairs.bin> [fov_deg]
```

Pass `--benchmark` (any position) to emit per-stage `[bench] stage=... us=...`
timing lines on stdout. `tools/benchmark.py` aggregates these across the
real-image fixtures and writes a CSV with median / p95 per stage.

## Generate Catalog and Test Data

```bash
cd tools
python3 generate_synthetic_data.py --quat 0 0 0 1 --out ../data/test_0
python3 generate_catalog.py
```

## Testing

### Unit Tests
Fast, in-memory tests covering all core algorithms. Runs automatically on every `git commit` (pre-commit hook).

```bash
cd build && ctest --output-on-failure
```

### Full Pipeline Test
End-to-end test: generates a synthetic star image, runs the tracker, verifies accuracy.

```bash
bash tests/test_full_pipeline.sh
```

### Monte Carlo Validation
Runs 50 random orientations and reports aggregate accuracy statistics.

```bash
bash tests/test_monte_carlo.sh
```

Both the full pipeline test and Monte Carlo run automatically on every `git push` via GitHub Actions.

## Accuracy Thresholds

| Test | Metric | Threshold |
|---|---|---|
| Full pipeline test | Angular error (identity quaternion) | < 2.0° |
| Monte Carlo | Per-trial success | error < 1.0° |
| Monte Carlo | Overall success rate | ≥ 80% of trials |

## Project Structure

```
src/                  C++ source (image processing, catalog, identification, estimation)
tools/                Python utilities (catalog builder, synthetic data generator, Monte Carlo)
tests/                Unit tests (GoogleTest) and integration test scripts
data/                 Generated catalog and test images (not tracked in git)
.github/workflows/    CI configuration
```

## Phases

- **Phase 1** — Baseline pipeline (complete)
- **Phase 2** — CI/CD and Monte Carlo validation (complete)
- **Phase 3** — Algorithm optimization (Gaussian centroiding, k-vector search, QUEST)
