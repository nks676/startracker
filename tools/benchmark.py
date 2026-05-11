#!/usr/bin/env python3
"""benchmark.py — micro-benchmark harness for the startracker pipeline.

Runs the compiled `startracker` binary multiple times per real-image fixture
with the `--benchmark` flag and aggregates per-stage timings (centroid,
catalog_load, identify, estimate, total). Median and p95 are printed; raw
results are also written to a CSV.

The fixture set mirrors `tools/test_real_images.py`: each `tests/data/
real_images/*.json` describes a TIFF + truth quaternion, and we reuse the
same TIFF->PNG calibration step so the binary is fed exactly the bytes that
the regression suite uses. By default we discover those fixtures
automatically; `--images` accepts an explicit glob for ad-hoc benchmarking
against arbitrary PNGs.

Usage:
    python tools/benchmark.py [--images <glob>] [--num-runs N]
                              [--csv <path>] [--binary <path>]

CSV columns:
    image, run, stage_centroid_us, stage_catalog_load_us,
    stage_identify_us, stage_estimate_us, stage_total_us

A startracker invocation that fails to identify (or otherwise exits
non-zero) does NOT abort the whole run — the offending row records NaN
for every stage and the script continues.
"""
from __future__ import annotations

import argparse
import csv
import glob as _glob
import json
import math
import re
import subprocess
import sys
from pathlib import Path

# Reuse the TIFF->PNG conversion and download helpers from the real-image
# regression so the benchmark sees exactly the same input bytes as the
# regression suite.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_real_images import (  # type: ignore[import-not-found]
    REPO_ROOT,
    download_if_missing,
    tiff_to_png,
)


_BENCH_RE = re.compile(r"^\[bench\]\s+stage=(\w+)\s+us=(\d+)\s*$")

STAGES = ("centroid", "catalog_load", "identify", "estimate", "total")


def parse_bench_lines(stdout: str) -> dict[str, float]:
    """Extract `[bench] stage=X us=Y` lines from stdout into a dict."""
    out: dict[str, float] = {}
    for line in stdout.splitlines():
        m = _BENCH_RE.match(line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def run_once(
    binary: Path, image: Path, fov_deg: float, cos_tol: float | None
) -> dict[str, float]:
    """Single binary invocation. Returns per-stage µs (NaN on failure)."""
    stars = REPO_ROOT / "data" / "catalog_stars.bin"
    pairs = REPO_ROOT / "data" / "catalog_pairs.bin"
    cmd = [str(binary), str(image), str(stars), str(pairs), f"{fov_deg:.6f}"]
    if cos_tol is not None:
        cmd.append(f"{cos_tol:.3e}")
    cmd.append("--benchmark")
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, check=False, timeout=60
        )
    except subprocess.TimeoutExpired:
        return {s: math.nan for s in STAGES}
    if proc.returncode != 0:
        return {s: math.nan for s in STAGES}
    parsed = parse_bench_lines(proc.stdout)
    return {s: parsed.get(s, math.nan) for s in STAGES}


def prepare_fixture_images(work_dir: Path) -> list[tuple[str, Path, float, float | None]]:
    """Download + convert every JSON fixture, returning the list of inputs
    (label, png_path, fov_deg, cos_tol) the benchmark should iterate."""
    fixtures = sorted(
        (REPO_ROOT / "tests" / "data" / "real_images").glob("*.json")
    )
    work_dir.mkdir(parents=True, exist_ok=True)
    inputs: list[tuple[str, Path, float, float | None]] = []
    for fx in fixtures:
        truth = json.loads(fx.read_text())
        url = truth["source_url"]
        name = url.rsplit("/", 1)[-1]
        tiff_path = work_dir / name
        png_path = tiff_path.with_suffix(".png")
        try:
            download_if_missing(url, tiff_path)
            tiff_to_png(tiff_path, png_path)
        except Exception as e:  # noqa: BLE001 — diagnostics only
            print(f"  WARN: could not prepare {fx.stem}: {e}", file=sys.stderr)
            continue
        inputs.append(
            (fx.stem, png_path, float(truth["fov_horizontal_deg"]),
             truth.get("cos_tol"))
        )
    return inputs


def discover_glob_inputs(pattern: str) -> list[tuple[str, Path, float, float | None]]:
    """For an arbitrary `--images` glob, default to FOV=20° and no cos_tol
    override (matches the binary's default behavior). Caller is responsible
    for ensuring the glob lines up with the catalog's FOV assumption."""
    paths = sorted(Path(p) for p in _glob.glob(pattern))
    return [(p.stem, p.resolve(), 20.0, None) for p in paths]


def percentile(xs: list[float], q: float) -> float:
    """Nearest-rank percentile, NaN-tolerant. q in [0, 100]."""
    clean = [x for x in xs if not math.isnan(x)]
    if not clean:
        return math.nan
    clean.sort()
    k = max(0, min(len(clean) - 1, int(math.ceil(q / 100.0 * len(clean))) - 1))
    return clean[k]


def median(xs: list[float]) -> float:
    clean = sorted(x for x in xs if not math.isnan(x))
    if not clean:
        return math.nan
    n = len(clean)
    if n % 2:
        return clean[n // 2]
    return 0.5 * (clean[n // 2 - 1] + clean[n // 2])


def fmt_us(x: float) -> str:
    return "  NaN" if math.isnan(x) else f"{x:8.0f}"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "build" / "startracker",
        help="Path to the compiled startracker binary.",
    )
    p.add_argument(
        "--images",
        type=str,
        default=None,
        help="Optional glob of PNGs to benchmark. Default: use the real-image "
             "regression fixtures (tests/data/real_images/*.json).",
    )
    p.add_argument(
        "--num-runs",
        type=int,
        default=5,
        help="Number of repeat runs per image (default 5).",
    )
    p.add_argument(
        "--csv",
        type=Path,
        default=REPO_ROOT / "data" / "benchmark_results.csv",
        help="Output CSV path (default data/benchmark_results.csv, gitignored).",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=REPO_ROOT / "data" / "real_images",
        help="Where downloaded TIFFs / generated PNGs are cached.",
    )
    args = p.parse_args()

    if not args.binary.exists():
        print(f"ERROR: binary not found at {args.binary}", file=sys.stderr)
        return 2

    if args.images:
        inputs = discover_glob_inputs(args.images)
        if not inputs:
            print(f"ERROR: glob '{args.images}' matched no files",
                  file=sys.stderr)
            return 2
    else:
        inputs = prepare_fixture_images(args.work_dir)
        if not inputs:
            print("ERROR: no fixtures available to benchmark", file=sys.stderr)
            return 2

    print(f"=== Benchmark: {len(inputs)} image(s), {args.num_runs} run(s) each ===")

    # Each row: (label, run_idx, dict of stage->us)
    rows: list[tuple[str, int, dict[str, float]]] = []
    per_image_per_stage: dict[str, dict[str, list[float]]] = {}
    overall: dict[str, list[float]] = {s: [] for s in STAGES}

    for label, png_path, fov_deg, cos_tol in inputs:
        per_image_per_stage[label] = {s: [] for s in STAGES}
        print(f"\n  {label}  ({png_path.name})")
        for run_idx in range(args.num_runs):
            timings = run_once(args.binary, png_path, fov_deg, cos_tol)
            rows.append((label, run_idx, timings))
            for s in STAGES:
                per_image_per_stage[label][s].append(timings[s])
                overall[s].append(timings[s])
            failed = math.isnan(timings["total"])
            tag = "FAIL " if failed else "     "
            print(
                f"    run {run_idx + 1}/{args.num_runs}  {tag}"
                f"centroid={fmt_us(timings['centroid'])}  "
                f"catalog={fmt_us(timings['catalog_load'])}  "
                f"identify={fmt_us(timings['identify'])}  "
                f"estimate={fmt_us(timings['estimate'])}  "
                f"total={fmt_us(timings['total'])}"
            )

    # CSV output
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow([
            "image", "run",
            "stage_centroid_us", "stage_catalog_load_us",
            "stage_identify_us", "stage_estimate_us", "stage_total_us",
        ])
        for label, run_idx, timings in rows:
            w.writerow([
                label, run_idx,
                timings["centroid"], timings["catalog_load"],
                timings["identify"], timings["estimate"], timings["total"],
            ])
    print(f"\nCSV written to {args.csv}")

    # Aggregate table
    header = f"{'image':<24} {'stage':<14} {'median_us':>11} {'p95_us':>11}"
    print("\n" + "=" * len(header))
    print(header)
    print("-" * len(header))
    for label in per_image_per_stage:
        for s in STAGES:
            xs = per_image_per_stage[label][s]
            print(
                f"{label:<24} {s:<14} {median(xs):>11.0f} "
                f"{percentile(xs, 95):>11.0f}"
            )
    print("-" * len(header))
    for s in STAGES:
        xs = overall[s]
        print(
            f"{'(all)':<24} {s:<14} {median(xs):>11.0f} "
            f"{percentile(xs, 95):>11.0f}"
        )
    print("=" * len(header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
