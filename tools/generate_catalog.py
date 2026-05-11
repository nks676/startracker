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

    # --- Mortari k-vector index ---
    # The pairs array is sorted descending by cos_val. The k-vector replaces
    # the two std::lower_bound calls in find_pairs with O(1) integer indexing.
    #
    # Construction:
    #   y_min = pairs[P-1].cos_val (smallest cos)
    #   y_max = pairs[0].cos_val   (largest cos)
    #   M = num_kvec_bins (we pick 4 * P)
    #   dq = (y_max - y_min) / M
    #   For each i in [0, M]:
    #     y_i = y_min + i * dq
    #     K[i] = largest index j (in the descending pairs array) such that
    #            pairs[j].cos_val >= y_i.
    #     Equivalently: all pairs at index <= K[i] satisfy cos_val >= y_i,
    #                   and all pairs at index >  K[i] satisfy cos_val <  y_i.
    #
    # File format (catalog_kvec.bin):
    #   [int32 num_kvec_bins M]
    #   [double y_min]
    #   [double y_max]
    #   [double dq]
    #   [(M + 1) * int32 K_i]
    num_pairs = len(pairs)
    if num_pairs > 0:
        cos_vals = np.array([p[0] for p in pairs], dtype=np.float64)  # descending
        y_max = float(cos_vals[0])
        y_min = float(cos_vals[-1])
        # Bin count: Mortari's original recipe is M = 4 * P, but for our
        # catalog that produces a ~100MB K array whose pointer-chasing wipes
        # out the algorithmic win. We instead size M so the bin step dq is
        # ~10x smaller than the tightest practical cos_tolerance (~1e-6),
        # giving dq ~ 1e-7. This keeps K small enough to live in L2 cache
        # while still bracketing the query range with negligible slack
        # (a single-bin overshoot adds <1 pair to the filter loop on
        # average). We still cap at 4*P and floor at a few thousand bins
        # for tiny catalogs.
        if y_max <= y_min:
            # Degenerate: a single distinct cos value. Fall back to one bin.
            M = 1
            dq = 1.0
        else:
            target_dq = 1e-7
            M_by_dq = int(np.ceil((y_max - y_min) / target_dq))
            M = min(4 * num_pairs, M_by_dq)
            M = max(M, 4096)
            dq = (y_max - y_min) / M

        # ascending = cos_vals reversed; for a target y, the largest j with
        # pairs[j].cos_val >= y in descending order equals
        # (P - 1) - (first index in ascending where val >= y).
        # We use numpy searchsorted on ascending values.
        asc = cos_vals[::-1]  # ascending
        # For each bin i compute y_i = y_min + i*dq
        i_arr = np.arange(M + 1, dtype=np.int64)
        y_i = y_min + i_arr * dq
        # searchsorted with side='left': first idx in asc such that asc[idx] >= y_i.
        # Number of pairs in descending array with cos_val >= y_i is
        # (P - left_idx). So K[i] (largest such index) = (P - left_idx) - 1.
        left_idx = np.searchsorted(asc, y_i, side='left')
        K = (num_pairs - left_idx).astype(np.int64) - 1
        # Clamp K to [-1, P-1]. K = -1 means "no pairs satisfy cos_val >= y_i"
        # (only possible when y_i > y_max due to float drift). Treat as P-1 ...
        # actually we keep -1 semantics: the C++ side will treat "first" lookup
        # as "K[i_high]+1 is start". To keep things simple, clamp >=0:
        # When all pairs have cos_val >= y_i, left_idx = 0 => K = P-1 (correct).
        # When no pairs have cos_val >= y_i, left_idx = P => K = -1.
        # The C++ side handles both by using K as an inclusive boundary.
        K = np.clip(K, -1, num_pairs - 1).astype(np.int32)

        kvec_file = os.path.join(out_dir, "catalog_kvec.bin")
        with open(kvec_file, 'wb') as f:
            f.write(struct.pack('<i', M))
            f.write(struct.pack('<ddd', y_min, y_max, dq))
            f.write(K.tobytes())
        print(f"Wrote k-vector index: M={M} bins, dq={dq:.6g}")

    print(f"Saved database to {out_dir}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    generate_database(max_mag=7.0, fov_max=25.0)
