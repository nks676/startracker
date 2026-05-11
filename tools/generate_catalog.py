import os
import struct
import numpy as np

def generate_database(max_mag=6.0, fov_max=25.0):
    cache_file = "hipparcos_cached.npy"
    if not os.path.exists(cache_file):
        raise RuntimeError("run generate_synthetic_data.py first to cache Hipparcos.")
        
    data = np.load(cache_file, allow_pickle=True).item()
    hip = data['hip']
    vmag = data['vmag']
    vectors = data['vectors']
    
    # Filter by mag
    mask = vmag <= max_mag
    hip = hip[mask]
    vmag = vmag[mask]
    vectors = vectors[mask]
    
    print(f"Generating DB for {len(hip)} stars...")
    
    # Save stars
    out_dir = "../data"
    os.makedirs(out_dir, exist_ok=True)
    
    star_file = os.path.join(out_dir, "catalog_stars.bin")
    # Format: int32 num_stars, then num_stars * (int32 id, double x, double y, double z)
    with open(star_file, 'wb') as f:
        f.write(struct.pack('<i', len(hip)))
        for i in range(len(hip)):
            # id, x, y, z
            f.write(struct.pack('<iddd', int(hip[i]), vectors[i,0], vectors[i,1], vectors[i,2]))
            
    # Pair generation
    # For baseline, we precompute all pairs within max_fov
    # and sort them by angular distance (cosine) to allow quick binary search or brute force
    cos_fov_max = np.cos(np.radians(fov_max))
    
    pairs = []
    
    # Fast pairwise dot products
    # To save memory, we can do it in blocks
    N = len(hip)
    block_size = 1000
    for i in range(0, N, block_size):
        end_i = min(i + block_size, N)
        vi = vectors[i:end_i]
        
        # dot product with all j > i to avoid memory blowup and duplicates
        for j in range(i, N, block_size):
            end_j = min(j + block_size, N)
            vj = vectors[j:end_j]
            
            dots = np.dot(vi, vj.T)
            
            # Diagonal and below where j_idx <= i_idx
            for k in range(dots.shape[0]):
                idx_i = i + k
                # If j block is same as i block, we only want j > i
                start_l = 0
                if i == j:
                    start_l = k + 1
                    
                for l in range(start_l, dots.shape[1]):
                    idx_j = j + l
                    
                    if dots[k, l] >= cos_fov_max:
                        pairs.append((dots[k, l], int(hip[idx_i]), int(hip[idx_j])))
                        
    print(f"Generated {len(pairs)} pairs within {fov_max} degrees.")
    
    # Sort pairs by cosine (descending, meaning ascending angular distance)
    # Actually, sorting them makes binary search by cosine distance very fast.
    pairs.sort(key=lambda x: x[0], reverse=True)
    
    pair_file = os.path.join(out_dir, "catalog_pairs.bin")
    # Format: int32 num_pairs, then num_pairs * (double cos_val, int32 id1, int32 id2)
    with open(pair_file, 'wb') as f:
        f.write(struct.pack('<i', len(pairs)))
        for p in pairs:
            f.write(struct.pack('<dii', p[0], p[1], p[2]))
            
    print(f"Saved database to {out_dir}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    generate_database(max_mag=7.0, fov_max=25.0)
