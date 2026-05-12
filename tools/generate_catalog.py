import os
import struct
import time
from itertools import combinations

import numpy as np

# --- Phase 3e.1: 4-star pattern hash catalog ---
#
# Canonical-order convention (Wave 2 of Phase 3e depends on this).
#
# Given a 4-star pattern, compute the 6 unordered pairwise angular distances.
# Sort the 6 edges ascending by distance (ties broken by lex order on the
# input-local star pair). For each of the 4 stars, its "edge signature" is
# the sorted 3-tuple of ranks (in 0..5) of the 3 sorted edges it is incident
# to. Sort the 4 stars ascending lexicographically by this signature; break
# remaining ties by the input identifier (HIP for catalog generation, centroid
# index for observed patterns) ascending.
#
# The resulting permutation is the canonical order. `hips[0..3]` in the binary
# stores HIPs in that order so an observer can reproduce the same ordering
# from centroid IDs and map position-by-position.

PATTERN_MAGIC = 0x50415431  # 'PAT1' little-endian
K_NEAREST = 8
QUANT_BITS = 10
QUANT_SCALE = (1 << QUANT_BITS) - 1  # 1023

# Edge enumeration in (i, j) local-index order; same order is used in C++
# helper for tests.
EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def canonical_order_and_key(vecs4, ids4):
    """Compute the canonical permutation and 64-bit key for a 4-star pattern.

    Args:
        vecs4: ndarray shape (4, 3), unit vectors of the 4 stars in input order.
        ids4: length-4 array-like of input identifiers (HIPs).

    Returns:
        (key, ordered_ids) where ordered_ids is a length-4 list of identifiers
        permuted into canonical order.
    """
    # 6 pairwise distances (angular = arccos(dot)). Use angle in radians; the
    # absolute scale cancels after ratio normalization.
    dots = np.clip(vecs4 @ vecs4.T, -1.0, 1.0)
    dists = np.empty(6, dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        dists[e] = float(np.arccos(dots[i, j]))

    # Sort edges by (distance, pair) ascending. Tie-break by pair to make the
    # sort fully deterministic when two distances are bit-equal (very rare).
    edge_order = sorted(range(6), key=lambda e: (dists[e], EDGE_PAIRS[e]))

    # Build per-star edge-rank signature.
    sig = [[], [], [], []]
    for rank, e in enumerate(edge_order):
        i, j = EDGE_PAIRS[e]
        sig[i].append(rank)
        sig[j].append(rank)
    for s in sig:
        s.sort()

    # Sort the 4 local-star indices by (signature, input_id) ascending.
    canonical = sorted(
        range(4), key=lambda s: (tuple(sig[s]), int(ids4[s]))
    )

    # Compute 5 ratios from sorted-by-distance edges and quantize.
    sorted_dists = dists[edge_order]
    largest = sorted_dists[5]
    # largest > 0 in any non-degenerate pattern (no coincident stars).
    ratios = sorted_dists[:5] / largest
    q = np.clip(np.round(ratios * QUANT_SCALE), 0, QUANT_SCALE).astype(np.uint64)

    # Pack 5 × 10-bit into uint64: q[0] in bits 0..9, ..., q[4] in bits 40..49.
    key = (q[4] << 40) | (q[3] << 30) | (q[2] << 20) | (q[1] << 10) | q[0]
    ordered_ids = [int(ids4[c]) for c in canonical]
    return int(key), ordered_ids


def _knn_indices(vectors, cos_threshold, k):
    """For each row of `vectors`, return up to `k` nearest-neighbor row indices
    whose cosine similarity exceeds `cos_threshold` (i.e., within the FOV
    radius). Output shape is (N, k_actual) variable-length as a python list.
    Excludes self.
    """
    N = vectors.shape[0]
    # All-pairs dot product. ~N^2 * 8 bytes; for N=15544 that's ~1.9 GB which
    # is too much. Process in row-blocks.
    block = 1024
    nn_lists = [None] * N
    for start in range(0, N, block):
        end = min(start + block, N)
        sims = vectors[start:end] @ vectors.T  # (block, N)
        # Mask self (diagonal entry per row).
        for r in range(end - start):
            sims[r, start + r] = -np.inf
        # For each row, keep cols above cos_threshold.
        for r in range(end - start):
            row = sims[r]
            cand = np.where(row > cos_threshold)[0]
            if cand.size == 0:
                nn_lists[start + r] = np.empty(0, dtype=np.int64)
                continue
            if cand.size <= k:
                # Sort all by descending sim.
                order = np.argsort(-row[cand])
                nn_lists[start + r] = cand[order]
            else:
                # Top-k by similarity.
                # argpartition for speed, then sort the top-k.
                part = np.argpartition(-row[cand], k - 1)[:k]
                top = cand[part]
                order = np.argsort(-row[top])
                nn_lists[start + r] = top[order]
    return nn_lists


def generate_pattern_catalog(fov_bin_deg, vectors, hip_ids, out_dir):
    """Generate one pattern-hash catalog binary for the given FOV bin.

    File: data/catalog_patterns_{fov_bin_deg}.bin
    Header (5 × int32): magic, fov_bin_deg, K_nearest, quant_bits, num_patterns.
    Body: num_patterns × (uint64 key, int32 hips[4]), sorted ascending by key.
    """
    t0 = time.time()
    half_fov_rad = np.radians(fov_bin_deg / 2.0)
    cos_threshold = np.cos(half_fov_rad)
    N = vectors.shape[0]
    print(
        f"[fov={fov_bin_deg}°] computing k-NN (k={K_NEAREST}) within {fov_bin_deg/2:.1f}° "
        f"radius for {N} stars..."
    )
    nn_lists = _knn_indices(vectors, cos_threshold, K_NEAREST)
    t_knn = time.time()
    print(f"[fov={fov_bin_deg}°] k-NN done in {t_knn - t0:.1f}s")

    # Estimated pattern count for log.
    total_subsets = 0
    for nn in nn_lists:
        if nn.size >= 3:
            n = nn.size
            total_subsets += n * (n - 1) * (n - 2) // 6
    print(f"[fov={fov_bin_deg}°] estimated patterns: {total_subsets}")

    # Allocate arrays up-front.
    keys = np.empty(total_subsets, dtype=np.uint64)
    hips_out = np.empty((total_subsets, 4), dtype=np.int32)

    write_idx = 0
    for s_idx in range(N):
        nn = nn_lists[s_idx]
        if nn.size < 3:
            continue
        s_vec = vectors[s_idx]
        s_hip = int(hip_ids[s_idx])
        nn_vecs = vectors[nn]
        nn_hips = hip_ids[nn]
        nn_count = nn.size
        for a, b, c in combinations(range(nn_count), 3):
            vecs4 = np.stack(
                (s_vec, nn_vecs[a], nn_vecs[b], nn_vecs[c]), axis=0
            )
            ids4 = (s_hip, int(nn_hips[a]), int(nn_hips[b]), int(nn_hips[c]))
            key, ordered = canonical_order_and_key(vecs4, ids4)
            keys[write_idx] = key
            hips_out[write_idx, 0] = ordered[0]
            hips_out[write_idx, 1] = ordered[1]
            hips_out[write_idx, 2] = ordered[2]
            hips_out[write_idx, 3] = ordered[3]
            write_idx += 1
    num_patterns = write_idx
    keys = keys[:num_patterns]
    hips_out = hips_out[:num_patterns]
    t_compute = time.time()
    print(
        f"[fov={fov_bin_deg}°] generated {num_patterns} patterns in "
        f"{t_compute - t_knn:.1f}s"
    )

    # Sort by key for binary-search lookup.
    order = np.argsort(keys, kind="stable")
    keys_sorted = keys[order]
    hips_sorted = hips_out[order]
    t_sort = time.time()

    # Write binary.
    out_path = os.path.join(out_dir, f"catalog_patterns_{fov_bin_deg}.bin")
    with open(out_path, "wb") as f:
        # Header: 5 × int32 little-endian.
        f.write(
            struct.pack(
                "<iiiii",
                PATTERN_MAGIC,
                int(fov_bin_deg),
                K_NEAREST,
                QUANT_BITS,
                num_patterns,
            )
        )
        # Body: pack via numpy structured array for speed.
        rec_dtype = np.dtype([
            ("key", "<u8"),
            ("hips", "<i4", (4,)),
        ])
        rec = np.empty(num_patterns, dtype=rec_dtype)
        rec["key"] = keys_sorted
        rec["hips"] = hips_sorted
        rec.tofile(f)
    t_write = time.time()
    size_mb = os.path.getsize(out_path) / 1e6
    print(
        f"[fov={fov_bin_deg}°] wrote {out_path} ({size_mb:.1f} MB) in "
        f"{t_write - t_sort:.1f}s; total {t_write - t0:.1f}s"
    )
    return out_path, num_patterns, t_write - t0


def generate_database(max_mag=6.0, fov_max=25.0,
                      pattern_fov_bins_deg=(10, 15, 20)):
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

    # --- Phase 3e.1: pattern-hash catalogs (one binary per FOV bin) ---
    for fov_bin_deg in pattern_fov_bins_deg:
        generate_pattern_catalog(fov_bin_deg, vectors, hip, out_dir)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    generate_database(max_mag=7.0, fov_max=25.0)
