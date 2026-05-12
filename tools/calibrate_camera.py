#!/usr/bin/env python3
"""calibrate_camera.py — fit a Brown-Conrady camera model (focal_x, focal_y,
center_x, center_y, k1, k2, k3, p1, p2) from a directory of calibration
frames + plate-solved truth quaternions.

For each frame the tool:
  1. Loads the image (TIFF/PNG) and extracts star centroids in Python with
     a simple tile-adaptive median/MAD threshold + connected-component
     moment centroid. The C++ pipeline has a more refined extractor but
     replicating it here keeps the optimizer self-contained and ~free of
     subprocess overhead.
  2. Projects every catalog star into the image plane using the truth
     quaternion (Inertial -> Camera) and a *nominal* pinhole intrinsic
     (focal from --fov, center at image midpoint, zero distortion).
  3. Pairs each projected catalog star with the nearest centroid inside a
     pixel-tolerance gate; rejects ambiguous pairings where the second-
     nearest centroid is also close.

The pooled correspondences across every frame feed a single non-linear
least-squares problem over 9 intrinsic parameters via
`scipy.optimize.least_squares`. Output is written as JSON with the same
field layout main.cpp's CameraModel uses (k1,k2,k3,p1,p2 — note the C++
CLI takes k1 k2 p1 p2 k3 in that argv order, but the JSON is keyed by
name so caller code never has to memorise the order).

Run against the committed ESA tetra3 fixtures as a smoke test:

    python tools/calibrate_camera.py \
        --fixtures-dir tests/data/real_images \
        --output data/esa_tetra3_camera.json

The tetra3 frames carry no real distortion data (the plate-solver fits a
zero-distortion model), so the smoke test exists only to confirm the
pipeline runs end-to-end and the fit converges to small distortion.
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


REPO_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------

def load_catalog_stars(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read data/catalog_stars.bin → (hips, unit_vectors).

    File layout (matches src/catalog.cpp loader):
        int32 num_stars
        num_stars × { int32 hip_id; double x; double y; double z }
    A handful of catalog entries have NaN coordinates (legitimate Hipparcos
    gaps); those rows are dropped silently.
    """
    raw = path.read_bytes()
    (n,) = struct.unpack_from("<i", raw, 0)
    dt = np.dtype([("hip", "<i4"), ("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
    arr = np.frombuffer(raw, dtype=dt, count=n, offset=4)
    vecs = np.stack([arr["x"], arr["y"], arr["z"]], axis=1)
    finite = np.isfinite(vecs).all(axis=1)
    return arr["hip"][finite].copy(), vecs[finite].copy()


# ---------------------------------------------------------------------------
# Image I/O + centroid extraction
# ---------------------------------------------------------------------------

def load_image(path: Path) -> np.ndarray:
    """Return a float32 single-channel image. Handles 8/16-bit TIFF and PNG."""
    img = np.asarray(Image.open(path))
    if img.ndim == 3:
        img = img[..., 0]
    return img.astype(np.float32)


def _tile_bg_sigma(im: np.ndarray, tile: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Per-tile median + MAD background/noise estimate. Returns per-pixel
    arrays of the same shape as `im` (constant within each tile)."""
    h, w = im.shape
    bg = np.empty_like(im)
    sg = np.empty_like(im)
    for y0 in range(0, h, tile):
        y1 = min(y0 + tile, h)
        for x0 in range(0, w, tile):
            x1 = min(x0 + tile, w)
            t = im[y0:y1, x0:x1]
            m = float(np.median(t))
            s = float(1.4826 * np.median(np.abs(t - m)))
            bg[y0:y1, x0:x1] = m
            sg[y0:y1, x0:x1] = max(s, 1.0)
    return bg, sg


def extract_centroids(
    image: np.ndarray, sigma_threshold: float = 4.0, min_pixels: int = 2
) -> np.ndarray:
    """Tile-adaptive threshold → connected components → intensity-weighted
    moments. Returns an (N,2) array of (x, y) pixel coordinates in image
    convention (origin top-left, y increases downward), sorted by peak
    brightness descending."""
    bg, sg = _tile_bg_sigma(image)
    above = (image - bg) > sigma_threshold * sg
    lbl, ncomp = ndimage.label(above)
    if ncomp == 0:
        return np.zeros((0, 2), dtype=np.float64)
    sizes = ndimage.sum_labels(above, lbl, range(1, ncomp + 1))
    keep_idx = np.where(sizes >= min_pixels)[0] + 1
    if keep_idx.size == 0:
        return np.zeros((0, 2), dtype=np.float64)
    weight = np.maximum(image - bg, 0.0)
    coms = ndimage.center_of_mass(weight, lbl, list(keep_idx))
    peaks = ndimage.maximum(weight, lbl, list(keep_idx))
    cents = np.array([(c[1], c[0]) for c in coms], dtype=np.float64)
    order = np.argsort(-np.asarray(peaks))
    return cents[order]


# ---------------------------------------------------------------------------
# Projection (Brown-Conrady forward model, mirrors src/identification.cpp)
# ---------------------------------------------------------------------------

def project_vectorized(
    v_cam: np.ndarray,
    fx: float,
    fy: float,
    cx: float,
    cy: float,
    k1: float,
    k2: float,
    k3: float,
    p1: float,
    p2: float,
) -> np.ndarray:
    """Vectorized port of src/identification.cpp::project. `v_cam` is
    (N,3) with z > 0; returns (N,2) pixel coordinates."""
    x = v_cam[:, 0] / v_cam[:, 2]
    y = v_cam[:, 1] / v_cam[:, 2]
    r2 = x * x + y * y
    radial = 1.0 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2
    x_d = x * radial + 2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x)
    y_d = y * radial + p1 * (r2 + 2.0 * y * y) + 2.0 * p2 * x * y
    px = fx * x_d + cx
    py = fy * y_d + cy
    return np.stack([px, py], axis=1)


# ---------------------------------------------------------------------------
# Per-frame matching
# ---------------------------------------------------------------------------

@dataclass
class FrameMatches:
    name: str
    R_cam_from_inertial: np.ndarray  # (3,3); v_cam = R @ v_inertial
    v_cam: np.ndarray                # (M,3) catalog unit vectors in camera frame
    pix_obs: np.ndarray              # (M,2) observed centroid pixel coords
    hips: np.ndarray                 # (M,) catalog HIP IDs (debug/audit only)


def match_frame(
    truth: dict,
    image_path: Path,
    catalog_hips: np.ndarray,
    catalog_vecs: np.ndarray,
    match_tol_px: float = 5.0,
    ambig_tol_px: float = 8.0,
) -> FrameMatches:
    """Build star↔centroid correspondences for one frame using the truth
    quaternion + nominal intrinsics. `ambig_tol_px` is the distance the
    second-nearest centroid must exceed for the match to be accepted."""
    q = np.asarray(truth["quaternion_xyzw"], dtype=np.float64)
    w = int(truth["image_width_px"])
    h = int(truth["image_height_px"])
    fov = float(truth["fov_horizontal_deg"])
    fx = w / (2.0 * np.tan(np.deg2rad(fov) / 2.0))
    fy = fx
    cx, cy = w / 2.0, h / 2.0

    # v_cam = R @ v_inertial (matches src/estimation.cpp's R = M_W * M_V^T).
    R = Rotation.from_quat(q).as_matrix()
    v_cam_all = catalog_vecs @ R.T
    in_front = v_cam_all[:, 2] > 0.2
    v_sub = v_cam_all[in_front]
    hip_sub = catalog_hips[in_front]

    proj = project_vectorized(v_sub, fx, fy, cx, cy, 0.0, 0.0, 0.0, 0.0, 0.0)
    on_image = (
        (proj[:, 0] > 0) & (proj[:, 0] < w) & (proj[:, 1] > 0) & (proj[:, 1] < h)
    )
    proj = proj[on_image]
    v_sub = v_sub[on_image]
    hip_sub = hip_sub[on_image]

    image = load_image(image_path)
    centroids = extract_centroids(image)

    if len(centroids) == 0 or len(proj) == 0:
        return FrameMatches(
            name=image_path.stem,
            R_cam_from_inertial=R,
            v_cam=np.zeros((0, 3)),
            pix_obs=np.zeros((0, 2)),
            hips=np.zeros((0,), dtype=np.int32),
        )

    # For each projected catalog star, find nearest centroid; accept if
    # within match_tol_px and second-nearest beyond ambig_tol_px.
    matched_v = []
    matched_px = []
    matched_hip = []
    used_centroid = set()
    for i in range(len(proj)):
        d = np.hypot(centroids[:, 0] - proj[i, 0], centroids[:, 1] - proj[i, 1])
        order = np.argsort(d)
        if d[order[0]] >= match_tol_px:
            continue
        if len(order) >= 2 and d[order[1]] < ambig_tol_px:
            continue
        if int(order[0]) in used_centroid:
            continue
        used_centroid.add(int(order[0]))
        matched_v.append(v_sub[i])
        matched_px.append(centroids[order[0]])
        matched_hip.append(int(hip_sub[i]))

    return FrameMatches(
        name=image_path.stem,
        R_cam_from_inertial=R,
        v_cam=np.asarray(matched_v).reshape(-1, 3),
        pix_obs=np.asarray(matched_px).reshape(-1, 2),
        hips=np.asarray(matched_hip, dtype=np.int32),
    )


# ---------------------------------------------------------------------------
# Optimisation
# ---------------------------------------------------------------------------

# 9-parameter vector layout. Centralised so the residual and the result-
# writer agree.
PARAM_NAMES = ("focal_x", "focal_y", "center_x", "center_y",
               "k1", "k2", "k3", "p1", "p2")


def pack_params(model: dict) -> np.ndarray:
    return np.array([model[k] for k in PARAM_NAMES], dtype=np.float64)


def unpack_params(p: np.ndarray) -> dict:
    return {name: float(p[i]) for i, name in enumerate(PARAM_NAMES)}


def residuals(p: np.ndarray, frames: list[FrameMatches]) -> np.ndarray:
    """Stacked (proj - observed) for every matched star across every
    frame. Returns a 1-D array of length 2*sum(M_i)."""
    fx, fy, cx, cy, k1, k2, k3, p1, p2 = p
    out = []
    for fm in frames:
        if len(fm.v_cam) == 0:
            continue
        proj = project_vectorized(fm.v_cam, fx, fy, cx, cy, k1, k2, k3, p1, p2)
        out.append((proj - fm.pix_obs).ravel())
    if not out:
        return np.zeros(0)
    return np.concatenate(out)


def rms_pixels(residual_vec: np.ndarray) -> float:
    if residual_vec.size == 0:
        return float("nan")
    # residual_vec is a flat (dx, dy, dx, dy, ...); RMS of per-star euclidean
    # error is sqrt(mean(dx^2 + dy^2)).
    r = residual_vec.reshape(-1, 2)
    return float(np.sqrt(np.mean(np.sum(r * r, axis=1))))


# ---------------------------------------------------------------------------
# Fixture handling
# ---------------------------------------------------------------------------

def download_if_missing(url: str, dest: Path) -> None:
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)


def resolve_image_for_fixture(truth: dict, fixtures_dir: Path, work_dir: Path) -> Path:
    """Truth fixtures reference `source_url`; the matching image is cached in
    work_dir (mirrors tools/test_real_images.py)."""
    url = truth["source_url"]
    name = url.rsplit("/", 1)[-1]
    dest = work_dir / name
    download_if_missing(url, dest)
    return dest


def initial_guess(truth_frames: list[dict]) -> dict:
    """Use the first frame's width/FOV. The fixtures we care about share a
    single physical camera so width and FOV match across frames; if they
    don't, the optimiser still converges, just from a slightly off start."""
    t0 = truth_frames[0]
    w = float(t0["image_width_px"])
    h = float(t0["image_height_px"])
    fx = w / (2.0 * np.tan(np.deg2rad(float(t0["fov_horizontal_deg"])) / 2.0))
    return {
        "focal_x": fx,
        "focal_y": fx,
        "center_x": w / 2.0,
        "center_y": h / 2.0,
        "k1": 0.0, "k2": 0.0, "k3": 0.0, "p1": 0.0, "p2": 0.0,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def calibrate(
    fixtures_dir: Path,
    catalog_path: Path,
    output_path: Path,
    work_dir: Path,
    match_tol_px: float = 5.0,
) -> dict:
    fixtures = sorted(fixtures_dir.glob("*.json"))
    if not fixtures:
        raise SystemExit(f"No truth fixtures in {fixtures_dir}")

    print(f"=== Calibrating from {len(fixtures)} fixture(s) ===")
    catalog_hips, catalog_vecs = load_catalog_stars(catalog_path)
    print(f"Catalog: {len(catalog_hips)} stars")

    work_dir.mkdir(parents=True, exist_ok=True)
    truths: list[dict] = []
    frames: list[FrameMatches] = []
    for fx_path in fixtures:
        truth = json.loads(fx_path.read_text())
        truths.append(truth)
        img_path = resolve_image_for_fixture(truth, fixtures_dir, work_dir)
        fm = match_frame(
            truth, img_path, catalog_hips, catalog_vecs, match_tol_px=match_tol_px
        )
        frames.append(fm)
        print(f"  {fx_path.stem}: {len(fm.v_cam)} matched star(s)")

    total_matches = sum(len(f.v_cam) for f in frames)
    if total_matches < len(PARAM_NAMES):
        raise SystemExit(
            f"Only {total_matches} total matches across all frames; need at "
            f"least {len(PARAM_NAMES)} to constrain 9 parameters."
        )

    guess = initial_guess(truths)
    p0 = pack_params(guess)
    r0 = residuals(p0, frames)
    rms_pre = rms_pixels(r0)
    print(f"\nInitial residual RMS: {rms_pre:.3f} px over {total_matches} matches")

    result = least_squares(residuals, p0, args=(frames,), method="lm", xtol=1e-10)
    fitted = unpack_params(result.x)
    rms_post = rms_pixels(result.fun)
    print(f"Final residual RMS:   {rms_post:.3f} px (success={result.success})")

    print("\nPer-coefficient deltas (post − pre):")
    for name in PARAM_NAMES:
        before = guess[name]
        after = fitted[name]
        if abs(before) > 1e-9:
            rel = f"  ({(after - before) / before * 100:+.3f}%)"
        else:
            rel = ""
        print(f"  {name:9s}  {before:14.6g} → {after:14.6g}   Δ={after-before:+.6g}{rel}")

    payload = {
        "model": "brown_conrady",
        "image_width_px": int(truths[0]["image_width_px"]),
        "image_height_px": int(truths[0]["image_height_px"]),
        **fitted,
        "metadata": {
            "fixtures": [f.stem for f in fixtures],
            "total_matches": int(total_matches),
            "rms_pre_px": rms_pre,
            "rms_post_px": rms_post,
            "match_tol_px": match_tol_px,
            "least_squares_success": bool(result.success),
            "least_squares_status": int(result.status),
            "least_squares_message": str(result.message),
            "initial_guess": guess,
            "cli_distortion_order_k1_k2_p1_p2_k3": [
                fitted["k1"], fitted["k2"], fitted["p1"], fitted["p2"], fitted["k3"]
            ],
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nWrote {output_path}")
    return payload


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--fixtures-dir",
        type=Path,
        default=REPO_ROOT / "tests" / "data" / "real_images",
        help="Directory of *.json truth fixtures (same schema as tests/data/real_images/).",
    )
    p.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "build" / "startracker",
        help=(
            "Path to the compiled startracker binary. Currently unused (centroids "
            "are extracted in Python), but accepted so the CLI matches the "
            "rest of the tools/ family — and future work can route through it."
        ),
    )
    p.add_argument(
        "--catalog-stars",
        type=Path,
        default=REPO_ROOT / "data" / "catalog_stars.bin",
        help="Path to the binary star catalog (defaults to data/catalog_stars.bin).",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "esa_tetra3_camera.json",
        help="Where to write the fitted camera-model JSON.",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=REPO_ROOT / "data" / "real_images",
        help="Where downloaded fixture images are cached.",
    )
    p.add_argument(
        "--match-tol-px",
        type=float,
        default=5.0,
        help="Max distance (pixels) between projected catalog star and centroid.",
    )
    args = p.parse_args()
    calibrate(
        fixtures_dir=args.fixtures_dir,
        catalog_path=args.catalog_stars,
        output_path=args.output,
        work_dir=args.work_dir,
        match_tol_px=args.match_tol_px,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
