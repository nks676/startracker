#pragma once

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

// Phase 3g.1: 4-star pattern-hash key computation, extracted from
// identification.cpp. Pure math on 4-tuples — no catalog, no DB, no camera
// model dependency. Mirrors tools/generate_catalog.py:canonical_order_and_key
// byte-for-byte; the kPatternMagic / kQuantBits constants in catalog.cpp
// guard the on-disk format and must remain consistent with kQuantBits here.

namespace pattern_hash {

// Edge enumeration in (i, j) local-index order. Must match
// tools/generate_catalog.py:EDGE_PAIRS exactly so canonical-order bit-flips
// between generator and runtime are impossible.
constexpr std::array<std::pair<int, int>, 6> kEdgePairs = {{
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
}};
constexpr int kQuantBits = 10;
constexpr int kQuantScale = (1 << kQuantBits) - 1; // 1023

} // namespace pattern_hash

// Computes the canonical-order 64-bit quantized geometric key for a 4-star
// pattern, mirroring tools/generate_catalog.py:canonical_order_and_key.
// `ids4` are the input identifiers (HIPs at catalog gen, centroid indices
// at runtime); they participate only in the lex tie-break, so the resulting
// canonical permutation is geometry-determined modulo ties.
//
// `out_canonical` (4 entries) receives the local indices [0..3] permuted into
// canonical order, so the caller can pair observed[ out_canonical[i] ] with
// hips[i] from a StarPattern record.
//
// Exposed at file scope (not the namespace) for back-compat with the
// pre-3g.1 signature used by test_identification.cpp.
uint64_t pattern_key_canonical(const std::array<std::array<double, 3>, 4> &vecs4,
                               const std::array<int, 4> &ids4,
                               std::array<int, 4> &out_canonical);

// Enumerate the set of pattern keys consistent with the input 4-tuple under
// centroid noise. Edge-rank flips at the canonical sort happen when the gap
// between adjacent sorted distances is below the noise floor; this function
// returns the union of all keys reachable by swapping such uncertain
// adjacent pairs. The first entry is always the noise-free canonical key
// (output of pattern_key_canonical). Deduplicated.
//
// `noise_tol` is in radians; values around the expected per-edge angle noise
// (a few × centroid_noise_rad) are appropriate.
std::vector<std::pair<uint64_t, std::array<int, 4>>>
pattern_keys_noise_robust(const std::array<std::array<double, 3>, 4> &vecs4,
                          const std::array<int, 4> &ids4, double noise_tol);
