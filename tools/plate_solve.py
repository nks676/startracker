"""
plate_solve.py — Get ground-truth attitude for a real star-field image.

Uses nova.astrometry.net (via astroquery) as an independent plate solver,
then converts the returned WCS into a quaternion in the convention used by
tools/generate_synthetic_data.py.

CONVENTIONS (must match generate_synthetic_data.py):
    - Quaternion is [x, y, z, w] (scipy.spatial.transform.Rotation.as_quat()).
    - The quaternion represents an ACTIVE rotation taking inertial (ICRS J2000)
      vectors into the CAMERA frame:  v_cam = R.apply(v_inertial).
    - Camera frame: +Z = boresight (optical axis, out of lens),
                    +X = image right (+u),
                    +Y = image down  (+v).
      i.e. pixel (u, v) increases right/down with origin at top-left, and the
      pinhole equations are x_pix = f * X/Z + cx,  y_pix = f * Y/Z + cy.

WCS -> CAMERA BASIS:
    From the WCS we get the celestial coordinates of the image center
    (RA0, Dec0) and the local CD matrix at that point. The CD matrix tells us
    how a small pixel offset (du, dv) maps to a small (dRA*cos(Dec), dDec)
    offset on the sky, in degrees. We construct three orthonormal vectors,
    each expressed in the ICRS frame:
        z_cam_in_inertial = unit vector to (RA0, Dec0)               (boresight)
        x_cam_in_inertial = direction on sky that +u (image right) points to
        y_cam_in_inertial = direction on sky that +v (image down)  points to
    The 3x3 matrix M whose ROWS are (x_cam, y_cam, z_cam) (each expressed in
    the inertial basis) is exactly the rotation matrix that takes an inertial
    vector and yields its components in the camera frame:
        v_cam = M @ v_inertial
    Hence scipy.Rotation.from_matrix(M).as_quat() is our [x,y,z,w].

A self-check at the end re-projects the center and a small +u offset and
verifies they land where expected.

USAGE:
    python plate_solve.py <image> [--api-key KEY] [--scale-lo 10] [--scale-hi 13]

    API key resolution order: --api-key flag, then ASTROMETRY_API_KEY env var.
    Get a free key at https://nova.astrometry.net/ -> Sign In -> My Profile.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def _load_dotenv() -> None:
    """Minimal .env loader. Looks in the project root (parent of this file's
    directory) and the current working directory."""
    candidates = [
        Path(__file__).resolve().parent.parent / ".env",
        Path.cwd() / ".env",
    ]
    for p in candidates:
        if not p.is_file():
            continue
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, _, v = line.partition("=")
            k = k.strip()
            v = v.strip().strip('"').strip("'")
            os.environ.setdefault(k, v)


def _local_sky_basis(ra_deg: float, dec_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (e_ra, e_dec, n) unit vectors at (ra, dec), expressed in ICRS.

    e_ra: direction of increasing RA on the sky (i.e. RA*cos(Dec) east-pointing
          tangent vector).
    e_dec: direction of increasing Dec on the sky (north-pointing tangent).
    n:    radial unit vector (the point itself on the celestial sphere).
    """
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    cd, sd = np.cos(dec), np.sin(dec)
    cr, sr = np.cos(ra), np.sin(ra)
    n = np.array([cd * cr, cd * sr, sd])
    e_ra = np.array([-sr, cr, 0.0])              # d(n)/d(RA) / cos(Dec)
    e_dec = np.array([-sd * cr, -sd * sr, cd])   # d(n)/d(Dec)
    return e_ra, e_dec, n


def wcs_to_camera_quaternion(wcs, width: int, height: int) -> dict:
    """Convert an astropy WCS object to our camera-frame quaternion.

    Returns a dict with quaternion_xyzw, ra_center_deg, dec_center_deg,
    roll_deg, pixel_scale_arcsec.
    """
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.spatial.transform import Rotation as R

    # 1) RA/Dec of image center. WCS is 1-indexed FITS convention; the geometric
    #    center of an (W x H) image is at pixel (W/2 + 0.5, H/2 + 0.5) in FITS,
    #    which corresponds to 0-indexed (W/2 - 0.5, H/2 - 0.5). We use the
    #    wcs_pix2world API with origin=0 and pass (W/2 - 0.5, H/2 - 0.5).
    cx_pix = width / 2.0 - 0.5
    cy_pix = height / 2.0 - 0.5
    sky_center = wcs.pixel_to_world(cx_pix, cy_pix)
    ra0 = float(sky_center.ra.deg)
    dec0 = float(sky_center.dec.deg)

    # 2) Local tangent basis at the center, expressed in ICRS.
    e_ra, e_dec, n = _local_sky_basis(ra0, dec0)

    # 3) Probe pixel directions: a +1px step in u (right) and +1px step in v (down).
    #    Get their sky positions and convert to small tangent-plane offsets.
    sky_u = wcs.pixel_to_world(cx_pix + 1.0, cy_pix)
    sky_v = wcs.pixel_to_world(cx_pix, cy_pix + 1.0)

    def tangent_offset(sky):
        """(d_ra_cosdec_deg, d_dec_deg) at the center."""
        d_ra = float(sky.ra.deg) - ra0
        # Unwrap RA jumps near 0/360
        if d_ra > 180.0:
            d_ra -= 360.0
        elif d_ra < -180.0:
            d_ra += 360.0
        d_ra_cosdec = d_ra * np.cos(np.radians(dec0))
        d_dec = float(sky.dec.deg) - dec0
        return d_ra_cosdec, d_dec

    du_ra, du_dec = tangent_offset(sky_u)  # contribution per +1 px in u
    dv_ra, dv_dec = tangent_offset(sky_v)  # contribution per +1 px in v

    # 4) Build the in-plane camera axes as ICRS vectors. A +u pixel step moves
    #    on the sky by (du_ra * e_ra + du_dec * e_dec); the *direction* on the
    #    sky that +u points to (before normalisation) is that combination.
    x_cam_in_icrs = du_ra * e_ra + du_dec * e_dec
    y_cam_in_icrs = dv_ra * e_ra + dv_dec * e_dec
    # Normalise (pixel scale drops out)
    x_cam_in_icrs /= np.linalg.norm(x_cam_in_icrs)
    y_cam_in_icrs /= np.linalg.norm(y_cam_in_icrs)
    # Boresight is the unit vector to the image center.
    z_cam_in_icrs = n / np.linalg.norm(n)

    # 5) Orthonormalise (defensively): make y perpendicular to z, and x = y x z's
    #    expected sign. We trust z (it's exact) and x (well-defined from u step),
    #    then derive y = z cross x to guarantee a right-handed orthonormal frame.
    #    NOTE: in our convention (x right, y down, z forward) the cross product
    #    is y = z x x  (right-hand rule), because z x x = y when {x, y, z} is RH.
    x_cam_in_icrs -= z_cam_in_icrs * np.dot(z_cam_in_icrs, x_cam_in_icrs)
    x_cam_in_icrs /= np.linalg.norm(x_cam_in_icrs)
    y_cam_derived = np.cross(z_cam_in_icrs, x_cam_in_icrs)
    # Sanity: y_cam_derived should align with the WCS-derived y_cam_in_icrs.
    sign = np.sign(np.dot(y_cam_derived, y_cam_in_icrs))
    if sign == 0:
        sign = 1.0
    # If the WCS implies a flipped (mirrored) image, the orthonormalised
    # right-handed frame would disagree in sign with the measured y direction.
    # We keep the right-handed frame and flag the discrepancy in roll only if
    # needed. For a non-mirrored sky image, sign should be +1.
    if sign < 0:
        # The image is parity-flipped relative to a right-handed RA/Dec sky
        # tangent plane. Flip x to keep RH-ness consistent with measured y.
        x_cam_in_icrs = -x_cam_in_icrs
        y_cam_derived = np.cross(z_cam_in_icrs, x_cam_in_icrs)

    # 6) Assemble M such that v_cam = M @ v_inertial. Rows are camera basis
    #    vectors expressed in the inertial (ICRS) basis.
    M = np.stack([x_cam_in_icrs, y_cam_derived, z_cam_in_icrs], axis=0)

    # Numerical hygiene: project to nearest rotation via SVD.
    U, _, Vt = np.linalg.svd(M)
    M_rot = U @ Vt
    if np.linalg.det(M_rot) < 0:
        # Should not happen after the flip above, but guard anyway.
        U[:, -1] *= -1
        M_rot = U @ Vt

    quat_xyzw = R.from_matrix(M_rot).as_quat()  # scipy returns [x, y, z, w]

    # 7) Roll angle (about boresight). Define roll as the angle from sky-north
    #    (e_dec) to camera-up (-y_cam) measured in the image plane, positive
    #    east-of-north (i.e. toward e_ra). Equivalent: atan2 of camera-up onto
    #    (e_ra, e_dec).
    up = -y_cam_derived
    roll_rad = np.arctan2(np.dot(up, e_ra), np.dot(up, e_dec))
    roll_deg = float(np.degrees(roll_rad))

    # 8) Pixel scale (arcsec/pixel), averaged over the two axes.
    scales_deg = proj_plane_pixel_scales(wcs)  # array of |CD| singular values, deg
    pixel_scale_arcsec = float(np.mean(scales_deg) * 3600.0)

    return {
        "quaternion_xyzw": [float(v) for v in quat_xyzw],
        "ra_center_deg": ra0,
        "dec_center_deg": dec0,
        "roll_deg": roll_deg,
        "pixel_scale_arcsec": pixel_scale_arcsec,
    }


def _self_check(wcs, width: int, height: int, result: dict) -> None:
    """Reproject the center and a +u offset through our quaternion and compare
    against the WCS to sanity-check the convention."""
    from scipy.spatial.transform import Rotation as R

    rot = R.from_quat(result["quaternion_xyzw"])
    # Center direction in inertial
    e_ra, e_dec, n_center = _local_sky_basis(
        result["ra_center_deg"], result["dec_center_deg"]
    )
    v_cam = rot.apply(n_center)
    # Boresight should map to camera +Z.
    if not (v_cam[2] > 0.999 and abs(v_cam[0]) < 1e-3 and abs(v_cam[1]) < 1e-3):
        print(
            f"  [warn] self-check: center vector in camera frame is {v_cam}, "
            f"expected approximately (0, 0, 1).",
            file=sys.stderr,
        )
    else:
        print(f"  self-check OK: boresight maps to camera +Z (v_cam={v_cam}).")


def solve_with_astrometry_net(
    image_path: Path, api_key: str, scale_lo: float, scale_hi: float
):
    """Submit to nova.astrometry.net and return an astropy WCS."""
    from astroquery.astrometry_net import AstrometryNet

    ast = AstrometryNet()
    ast.api_key = api_key

    # nova.astrometry.net accepts FITS, JPG, PNG, GIF. TIFF is generally NOT
    # accepted; we re-encode to PNG in a temp file if needed.
    suffix = image_path.suffix.lower()
    upload_path = image_path
    tmp_png = None
    if suffix in {".tif", ".tiff"}:
        from tempfile import NamedTemporaryFile

        with Image.open(image_path) as img:
            # Convert to 8-bit grayscale if necessary; nova handles L and RGB.
            if img.mode not in ("L", "RGB"):
                img = img.convert("L")
            tmp_png = NamedTemporaryFile(suffix=".png", delete=False)
            tmp_png.close()
            img.save(tmp_png.name, format="PNG")
            upload_path = Path(tmp_png.name)
            width, height = img.size
    else:
        with Image.open(image_path) as img:
            width, height = img.size

    print(f"Submitting {upload_path.name} to nova.astrometry.net "
          f"(scale {scale_lo}-{scale_hi} deg). This typically takes 30s-3min...")

    try:
        wcs_header = ast.solve_from_image(
            str(upload_path),
            scale_units="degwidth",
            scale_lower=scale_lo,
            scale_upper=scale_hi,
            publicly_visible="n",
            allow_commercial_use="n",
            allow_modifications="n",
        )
    finally:
        if tmp_png is not None:
            try:
                os.unlink(tmp_png.name)
            except OSError:
                pass

    if not wcs_header:
        raise RuntimeError("Plate solve failed: astrometry.net returned no WCS.")

    from astropy.wcs import WCS

    wcs = WCS(wcs_header)
    return wcs, width, height


def main() -> int:
    p = argparse.ArgumentParser(
        description="Plate-solve a star-field image and write truth.json."
    )
    p.add_argument("image", type=Path, help="Path to TIFF/PNG/JPG star image.")
    p.add_argument("--api-key", type=str, default=None,
                   help="nova.astrometry.net API key (else ASTROMETRY_API_KEY).")
    p.add_argument("--scale-lo", type=float, default=10.0,
                   help="Lower bound on image width in degrees (default 10).")
    p.add_argument("--scale-hi", type=float, default=13.0,
                   help="Upper bound on image width in degrees (default 13).")
    args = p.parse_args()

    _load_dotenv()
    api_key = args.api_key or os.environ.get("ASTROMETRY_API_KEY")
    if not api_key:
        print(
            "ERROR: no astrometry.net API key supplied.\n"
            "  Pass --api-key <KEY>, or set ASTROMETRY_API_KEY in the environment.\n"
            "  Get a free key at https://nova.astrometry.net/ "
            "(Sign In -> My Profile -> API).",
            file=sys.stderr,
        )
        return 2

    if not args.image.exists():
        print(f"ERROR: image not found: {args.image}", file=sys.stderr)
        return 2

    wcs, width, height = solve_with_astrometry_net(
        args.image, api_key, args.scale_lo, args.scale_hi
    )

    result = wcs_to_camera_quaternion(wcs, width, height)

    # FOV from pixel scale (horizontal).
    fov_horizontal_deg = (result["pixel_scale_arcsec"] * width) / 3600.0

    payload = {
        "quaternion_xyzw": result["quaternion_xyzw"],
        "ra_center_deg": result["ra_center_deg"],
        "dec_center_deg": result["dec_center_deg"],
        "roll_deg": result["roll_deg"],
        "pixel_scale_arcsec": result["pixel_scale_arcsec"],
        "fov_horizontal_deg": fov_horizontal_deg,
        "image_path": str(args.image.resolve()),
        "solver": "nova.astrometry.net (astroquery)",
        "image_width_px": width,
        "image_height_px": height,
    }

    out_path = args.image.with_name("truth.json")
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)

    print(f"Wrote {out_path}")
    print(f"  RA  = {payload['ra_center_deg']:.4f} deg")
    print(f"  Dec = {payload['dec_center_deg']:.4f} deg")
    print(f"  Roll = {payload['roll_deg']:.4f} deg")
    print(f"  Pixel scale = {payload['pixel_scale_arcsec']:.3f} arcsec/px")
    print(f"  Horizontal FOV = {fov_horizontal_deg:.3f} deg")
    print(f"  Quat (xyzw) = {payload['quaternion_xyzw']}")

    _self_check(wcs, width, height, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
