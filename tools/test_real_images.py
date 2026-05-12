#!/usr/bin/env python3
"""test_real_images.py — regression test for the startracker pipeline on
real camera images. For each truth fixture in tests/data/real_images/*.json,
downloads the referenced TIFF (cached on disk), feeds it directly to the
startracker binary (which now reads 16-bit grayscale TIFF natively), and
asserts the recovered attitude is within attitude_tolerance_deg of the
plate-solved truth quaternion.

Truth quaternions were computed once via nova.astrometry.net and committed;
no astrometry.net access is needed at test time.

Usage:
    python tools/test_real_images.py [--binary PATH] [--fixtures-dir PATH]
                                     [--work-dir PATH] [--keep]
Returns exit code 0 if every fixture passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image  # noqa: F401  (kept for tiff_to_png back-compat helper)


REPO_ROOT = Path(__file__).resolve().parent.parent


def attitude_error_deg(q_est: np.ndarray, q_truth: np.ndarray) -> float:
    """Geodesic angular error between two unit quaternions (xyzw), in degrees."""
    return float(2.0 * np.degrees(np.arccos(min(1.0, abs(q_est @ q_truth)))))


def download_if_missing(url: str, dest: Path) -> None:
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)


def tiff_to_png(tiff_path: Path, png_path: Path) -> None:
    """DEPRECATED but kept as a public helper for any external caller still on
    the 8-bit path. The C++ binary now reads 16-bit TIFFs natively, so the
    regression and benchmark drivers no longer call this — they hand the raw
    TIFF straight to the binary."""
    img = np.asarray(Image.open(tiff_path))
    if img.dtype == np.uint8:
        Image.fromarray(img, mode="L").save(png_path)
        return
    img = img.astype(np.float32)
    bg = float(np.median(img))
    sigma = float(np.std(img))
    scale = 10.0 / max(1.0, sigma)
    out = np.clip((img - bg) * scale + 50.0, 0, 255).astype(np.uint8)
    Image.fromarray(out, mode="L").save(png_path)


# Parses the C++ stdout for the estimated quaternion. Tolerant of whitespace.
_QUAT_RE = re.compile(
    r"Estimated Quaternion:\s*\[([^,]+),\s*([^,]+),\s*([^,]+),\s*([^\]]+)\]"
)


def run_startracker(
    binary: Path, image: Path, fov_deg: float, cos_tol: float | None = None
) -> np.ndarray:
    """Run the binary and parse the estimated quaternion. Raises on failure."""
    stars = REPO_ROOT / "data" / "catalog_stars.bin"
    pairs = REPO_ROOT / "data" / "catalog_pairs.bin"
    cmd = [str(binary), str(image), str(stars), str(pairs), f"{fov_deg:.6f}"]
    if cos_tol is not None:
        cmd.append(f"{cos_tol:.3e}")
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"startracker exited {proc.returncode}\nstdout: {proc.stdout}\n"
            f"stderr: {proc.stderr}"
        )
    m = _QUAT_RE.search(proc.stdout)
    if not m:
        raise RuntimeError(
            "startracker did not print an estimated quaternion.\n"
            f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return np.array([float(x) for x in m.groups()])


def run_fixture(fixture_path: Path, binary: Path, work_dir: Path) -> tuple[bool, str]:
    truth = json.loads(fixture_path.read_text())
    q_truth = np.array(truth["quaternion_xyzw"])
    tol_deg = float(truth.get("attitude_tolerance_deg", 0.5))
    fov_deg = float(truth["fov_horizontal_deg"])
    # Optional per-fixture tolerance override for the catalog pair matcher.
    # Real cameras have FOV-calibration drift (~0.4% on the tetra3 IMX265
    # sample) that pushes correct pair-cos values just outside the synthetic
    # default of 1e-5. The override lets each fixture stay self-contained.
    cos_tol = truth.get("cos_tol")

    url = truth["source_url"]
    name = url.rsplit("/", 1)[-1]
    tiff_path = work_dir / name
    download_if_missing(url, tiff_path)

    # 3a.6: feed the raw TIFF to the binary; it reads 16-bit grayscale TIFF
    # natively via src/tiff_reader.cpp. No PNG intermediate.
    q_est = run_startracker(binary, tiff_path, fov_deg, cos_tol)
    err = attitude_error_deg(q_est, q_truth)
    ok = err <= tol_deg
    return ok, f"{fixture_path.stem}: err={err:.4f}° (tolerance {tol_deg:.2f}°)"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "build" / "startracker",
        help="Path to the compiled startracker binary.",
    )
    p.add_argument(
        "--fixtures-dir",
        type=Path,
        default=REPO_ROOT / "tests" / "data" / "real_images",
        help="Directory of *.json truth fixtures.",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=REPO_ROOT / "data" / "real_images",
        help="Where downloaded TIFFs and generated PNGs are cached.",
    )
    args = p.parse_args()

    if not args.binary.exists():
        print(f"ERROR: binary not found at {args.binary}", file=sys.stderr)
        return 2

    fixtures = sorted(args.fixtures_dir.glob("*.json"))
    if not fixtures:
        print(f"ERROR: no fixtures in {args.fixtures_dir}", file=sys.stderr)
        return 2

    print(f"=== Real-image regression: {len(fixtures)} fixture(s) ===")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    all_ok = True
    for fx in fixtures:
        try:
            ok, msg = run_fixture(fx, args.binary, args.work_dir)
        except Exception as e:
            print(f"  {fx.stem}: ERROR — {e}")
            all_ok = False
            continue
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {msg}")
        if not ok:
            all_ok = False

    print("=" * 50)
    if all_ok:
        print("ALL FIXTURES PASSED")
        return 0
    print("SOME FIXTURES FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
