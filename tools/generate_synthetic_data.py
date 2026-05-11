import os
import json
import argparse
import numpy as np
from PIL import Image
from scipy.spatial.transform import Rotation as R
from astroquery.vizier import Vizier
import astropy.units as u

def fetch_hipparcos(max_mag=7.0):
    """
    Fetches the Hipparcos catalog from Vizier, keeping stars brighter than max_mag.
    Uses caching to avoid repeated downloads.
    """
    cache_file = "hipparcos_cached.npy"
    if os.path.exists(cache_file):
        print(f"Loading Hipparcos catalog from {cache_file}...")
        return np.load(cache_file, allow_pickle=True).item()
        
    print(f"Downloading Hipparcos catalog (Vmag <= {max_mag})...")
    v = Vizier(columns=['HIP', 'RAICRS', 'DEICRS', 'Vmag'],
               row_limit=-1)
    
    # Query Hipparcos Main Catalog (I/239)
    result = v.query_constraints(catalog='I/239/hip_main', Vmag=f'<={max_mag}')
    if not result:
        raise RuntimeError("Failed to download Hipparcos catalog.")
        
    table = result[0]
    
    # Extract columns
    hip_ids = table['HIP'].data.data
    ra_deg = table['RAICRS'].data.data
    dec_deg = table['DEICRS'].data.data
    vmag = table['Vmag'].data.data
    
    # Pre-compute unit vectors in ICRS
    ra_rad = np.radians(ra_deg)
    dec_rad = np.radians(dec_deg)
    
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    
    vectors = np.stack((x, y, z), axis=-1)
    
    data = {
        'hip': hip_ids,
        'vmag': vmag,
        'vectors': vectors
    }
    np.save(cache_file, data)
    print(f"Loaded {len(hip_ids)} stars and cached to {cache_file}.")
    return data

def generate_image(catalog, quat, resolution, fov_deg, noise_std, output_dir):
    """
    Generates a synthetic star tracker image.
    quat: [x, y, z, w] quaternion transforming INERTIAL to CAMERA frame.
    """
    rot = R.from_quat(quat)
    
    # Calculate focal length in pixels using FOV
    # fov_deg is the horizontal FoV usually, let's assume it's diagonal or width.
    # width_pixels = resolution[0]
    # fov = 2 * arctan( width / (2 * f) )  => f = width / (2 * tan(fov/2))
    width, height = resolution
    f_pixels = width / (2 * np.tan(np.radians(fov_deg) / 2))
    
    cx, cy = width / 2.0, height / 2.0
    
    vectors = catalog['vectors']
    vmag = catalog['vmag']
    hip = catalog['hip']
    
    # Transform inertial to camera
    # If quat transforms inertial IN to camera, then v_cam = R * v_inertial
    cam_vectors = rot.apply(vectors)
    
    # Filter stars behind the camera (Z <= 0)
    # Our camera frame convention: Z is optical axis (+), X is right, Y is down
    front_mask = cam_vectors[:, 2] > 0
    cam_vectors = cam_vectors[front_mask]
    vmag = vmag[front_mask]
    hip = hip[front_mask]
    
    # Pinhole projection
    x_proj = f_pixels * (cam_vectors[:, 0] / cam_vectors[:, 2])
    y_proj = f_pixels * (cam_vectors[:, 1] / cam_vectors[:, 2])
    
    x_pix = x_proj + cx
    y_pix = y_proj + cy
    
    # Filter stars outside image bounds
    margin = 10 # render margin
    in_view_mask = (x_pix > -margin) & (x_pix < width + margin) & \
                   (y_pix > -margin) & (y_pix < height + margin)
    
    x_pix = x_pix[in_view_mask]
    y_pix = y_pix[in_view_mask]
    vmag = vmag[in_view_mask]
    hip = hip[in_view_mask]
    cam_z = cam_vectors[:, 2][in_view_mask]
    
    truth_data = {
        'quaternion_xyzw': quat,
        'matches': []
    }
    
    # Create blank image (float for rendering)
    image = np.zeros((height, width), dtype=np.float32)
    
    # Render 2D Gaussian spots
    # Vmag to intensity:
    # A difference of 5 magnitudes is a factor of 100 in brightness.
    # Let vmag 0 = intensity 1e5 (arbitrary scale)
    base_intensity = 50000.0
    intensities = base_intensity * (10 ** (-0.4 * vmag))
    
    sigma = 1.0 # Gaussian blur sigma in pixels
    
    for i in range(len(x_pix)):
        xp, yp = x_pix[i], y_pix[i]
        intensity = intensities[i]
        
        # bounding box for gaussian
        x0 = max(0, int(xp - 4*sigma))
        x1 = min(width, int(xp + 4*sigma) + 1)
        y0 = max(0, int(yp - 4*sigma))
        y1 = min(height, int(yp + 4*sigma) + 1)
        
        if x0 >= x1 or y0 >= y1:
            continue
            
        yy, xx = np.mgrid[y0:y1, x0:x1]
        gaussian = np.exp(-((xx - xp)**2 + (yy - yp)**2) / (2 * sigma**2))
        
        # Normalize gaussian so sum is 1, then multiply by total intensity
        gaussian = gaussian / (2 * np.pi * sigma**2)
        
        image[y0:y1, x0:x1] += intensity * gaussian
        
        truth_data['matches'].append({
            'hip': int(hip[i]),
            'x': float(xp),
            'y': float(yp),
            'vmag': float(vmag[i])
        })
        
    print(f"Rendered {len(x_pix)} stars in field of view.")
    
    # Add noise
    # Base background level
    background = 50.0
    image += background
    
    # Poisson noise
    image = np.random.poisson(image)
    
    # Add optional Gaussian read noise
    if noise_std > 0:
        image = image + np.random.normal(0, noise_std, image.shape)
        
    image = np.clip(image, 0, 255).astype(np.uint8)
    
    os.makedirs(output_dir, exist_ok=True)
    img_path = os.path.join(output_dir, "synthetic_starfield.png")
    Image.fromarray(image, mode='L').save(img_path)
    
    json_path = os.path.join(output_dir, "truth.json")
    with open(json_path, 'w') as f:
        json.dump(truth_data, f, indent=2)
        
    print(f"Saved image to {img_path} and truth data to {json_path}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic star tracker images.")
    parser.add_argument("--quat", type=float, nargs=4, default=[0, 0, 0, 1], help="Quaternion [x, y, z, w] transforming inertial to camera.")
    parser.add_argument("--res", type=int, nargs=2, default=[1024, 1024], help="Resolution [width, height].")
    parser.add_argument("--fov", type=float, default=20.0, help="Horizontal Field of View (degrees).")
    parser.add_argument("--noise", type=float, default=5.0, help="Standard deviation of Gaussian read noise.")
    parser.add_argument("--out", type=str, default="../data/test_0", help="Output directory.")
    
    args = parser.parse_args()
    
    # Change dir to script dir for caching
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    catalog = fetch_hipparcos()
    
    generate_image(
        catalog=catalog,
        quat=args.quat,
        resolution=args.res,
        fov_deg=args.fov,
        noise_std=args.noise,
        output_dir=os.path.abspath(args.out)
    )
