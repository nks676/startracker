#!/usr/bin/env python3
"""
monte_carlo.py — Monte Carlo accuracy validation for the startracker engine.

Generates random orientations, runs the C++ tracker on each, and reports
aggregate accuracy statistics. Exits with code 1 if success rate < 80%.

Usage:
    python monte_carlo.py [--num-trials N] [--fov DEG] [--res W H] [--noise STD]
"""

import os
import sys
import csv
import argparse
import subprocess
import tempfile
import numpy as np
from scipy.spatial.transform import Rotation

# Import directly from generate_synthetic_data to avoid subprocess overhead
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from generate_synthetic_data import fetch_hipparcos, generate_image

# Thresholds
SUCCESS_THRESHOLD_DEG = 1.0   # Trial is a success if error < 1°
FAILURE_THRESHOLD_DEG = 5.0   # Trial is a failure if error > 5°
MIN_SUCCESS_RATE      = 0.80  # CI fails if success rate < 80%


def random_quaternion():
    """Generate a random unit quaternion with canonical positive-w form."""
    q = np.random.randn(4)
    q /= np.linalg.norm(q)
    if q[3] < 0:
        q = -q
    return q


def run_tracker(binary_path, image_path, catalog_stars, catalog_pairs, fov_deg):
    """
    Run the C++ startracker binary and return (estimated_quat, stars_identified).
    Returns (None, 0) on failure.
    """
    result = subprocess.run(
        [binary_path, image_path, catalog_stars, catalog_pairs, str(fov_deg)],
        capture_output=True, text=True
    )

    stdout = result.stdout
    est_quat = None
    stars_identified = 0

    for line in stdout.splitlines():
        if line.startswith("Estimated Quaternion:"):
            inner = line.split("[")[1].split("]")[0]
            est_quat = np.array([float(x.strip()) for x in inner.split(",")])
        if line.startswith("Identified"):
            try:
                stars_identified = int(line.split()[1])
            except (IndexError, ValueError):
                pass

    return est_quat, stars_identified


def angular_error_deg(q_truth, q_est):
    """Compute angular error in degrees between two quaternions [x,y,z,w]."""
    r_truth = Rotation.from_quat(q_truth)
    r_est   = Rotation.from_quat(q_est)
    r_err   = r_truth * r_est.inv()
    return np.degrees(r_err.magnitude())


def main():
    parser = argparse.ArgumentParser(description="Monte Carlo accuracy validation.")
    parser.add_argument("--num-trials", type=int, default=50)
    parser.add_argument("--fov",        type=float, default=20.0)
    parser.add_argument("--res",        type=int, nargs=2, default=[1024, 1024])
    parser.add_argument("--noise",      type=float, default=5.0)
    args = parser.parse_args()

    repo_root     = os.path.dirname(script_dir)
    binary_path   = os.path.join(repo_root, "build", "startracker")
    catalog_stars = os.path.join(repo_root, "data", "catalog_stars.bin")
    catalog_pairs = os.path.join(repo_root, "data", "catalog_pairs.bin")
    csv_out       = os.path.join(repo_root, "data", "monte_carlo_results.csv")

    # Validate required files exist
    for path, label in [(binary_path, "startracker binary"),
                         (catalog_stars, "catalog_stars.bin"),
                         (catalog_pairs, "catalog_pairs.bin")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found at {path}")
            sys.exit(1)

    # Load Hipparcos catalog once (uses cache if available)
    os.chdir(script_dir)
    print("Loading Hipparcos catalog...")
    catalog = fetch_hipparcos()

    print(f"\nRunning {args.num_trials} Monte Carlo trials "
          f"(FOV={args.fov}°, res={args.res[0]}x{args.res[1]}, noise={args.noise})\n")

    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for i in range(args.num_trials):
            q_truth = random_quaternion()

            generate_image(
                catalog=catalog,
                quat=q_truth.tolist(),
                resolution=args.res,
                fov_deg=args.fov,
                noise_std=args.noise,
                output_dir=tmpdir
            )

            image_path = os.path.join(tmpdir, "synthetic_starfield.png")
            q_est, stars_identified = run_tracker(
                binary_path, image_path, catalog_stars, catalog_pairs, args.fov
            )

            if q_est is None:
                error_deg = None
                status = "failed"
            else:
                q_est /= np.linalg.norm(q_est)
                error_deg = angular_error_deg(q_truth, q_est)
                status = "success" if error_deg < SUCCESS_THRESHOLD_DEG else "failed"

            results.append({
                "trial":            i + 1,
                "truth_quat":       q_truth.tolist(),
                "est_quat":         q_est.tolist() if q_est is not None else None,
                "angular_error_deg": round(error_deg, 6) if error_deg is not None else None,
                "stars_identified": stars_identified,
                "status":           status,
            })

            marker = "." if status == "success" else "F"
            print(f"  Trial {i+1:3d}/{args.num_trials}: "
                  f"{'%.4f°' % error_deg if error_deg is not None else 'NO SOLUTION':>10s}  "
                  f"stars={stars_identified}  {marker}")

    # --- Statistics ---
    successful = [r for r in results if r["status"] == "success"]
    failed     = [r for r in results if r["status"] == "failed"]
    success_rate = len(successful) / args.num_trials

    errors = [r["angular_error_deg"] for r in successful if r["angular_error_deg"] is not None]

    print(f"\n{'='*50}")
    print(f"MONTE CARLO RESULTS  ({args.num_trials} trials)")
    print(f"{'='*50}")
    print(f"Success rate : {success_rate*100:.1f}%  ({len(successful)}/{args.num_trials})"
          f"  [threshold: error < {SUCCESS_THRESHOLD_DEG}°]")

    if errors:
        print(f"Angular error (successful trials):")
        print(f"  Median : {np.median(errors):.4f}°")
        print(f"  Mean   : {np.mean(errors):.4f}°")
        print(f"  95th % : {np.percentile(errors, 95):.4f}°")
        print(f"  Max    : {np.max(errors):.4f}°")

    if failed:
        print(f"\nFailure cases ({len(failed)}):")
        for r in failed:
            err_str = f"{r['angular_error_deg']:.4f}°" if r["angular_error_deg"] is not None else "NO SOLUTION"
            print(f"  Trial {r['trial']:3d}: error={err_str}, stars={r['stars_identified']}, "
                  f"truth={[round(x,4) for x in r['truth_quat']]}")

    # --- Save CSV ---
    os.makedirs(os.path.dirname(csv_out), exist_ok=True)
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["trial", "truth_quat", "est_quat",
                                                "angular_error_deg", "stars_identified", "status"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to {csv_out}")

    # --- CI gate ---
    if success_rate < MIN_SUCCESS_RATE:
        print(f"\nFAILED: Success rate {success_rate*100:.1f}% is below "
              f"minimum {MIN_SUCCESS_RATE*100:.0f}%")
        sys.exit(1)

    print(f"\nPASSED: Success rate {success_rate*100:.1f}% meets minimum "
          f"{MIN_SUCCESS_RATE*100:.0f}%")


if __name__ == "__main__":
    main()
