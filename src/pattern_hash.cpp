#include "pattern_hash.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>

// Phase 3g.1: extracted from identification.cpp. Pure math on 4-tuples;
// no catalog, DB, or camera-model dependency. Behaviour preserved byte-
// for-byte from the prior identification.cpp definitions.

using pattern_hash::kEdgePairs;
using pattern_hash::kQuantScale;

uint64_t
pattern_key_canonical(const std::array<std::array<double, 3>, 4> &vecs4,
                      const std::array<int, 4> &ids4,
                      std::array<int, 4> &out_canonical) {
  // 1. Six pairwise angular distances (radians, acos of clamped dot).
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    double d = vecs4[i][0] * vecs4[j][0] + vecs4[i][1] * vecs4[j][1] +
               vecs4[i][2] * vecs4[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }

  // 2. Sort the 6 edges by (distance, EDGE_PAIRS) ascending. Tie-break by
  // the EDGE_PAIRS (i, j) tuple — equivalent to Python's
  // sorted(..., key=lambda e: (dists[e], EDGE_PAIRS[e])).
  std::array<int, 6> edge_order = {0, 1, 2, 3, 4, 5};
  std::sort(edge_order.begin(), edge_order.end(), [&](int a, int b) {
    if (dists[a] != dists[b]) return dists[a] < dists[b];
    if (kEdgePairs[a].first != kEdgePairs[b].first)
      return kEdgePairs[a].first < kEdgePairs[b].first;
    return kEdgePairs[a].second < kEdgePairs[b].second;
  });

  // 3. Per-star edge-rank signature: sorted 3-tuple of ranks (in 0..5) of the
  // 3 edges incident to each star.
  std::array<std::array<int, 3>, 4> sig{};
  std::array<int, 4> sig_fill = {0, 0, 0, 0};
  for (int rank = 0; rank < 6; ++rank) {
    int e = edge_order[rank];
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    sig[i][sig_fill[i]++] = rank;
    sig[j][sig_fill[j]++] = rank;
  }
  for (int s = 0; s < 4; ++s)
    std::sort(sig[s].begin(), sig[s].end());

  // 4. Sort the 4 local-star indices ascending by (signature, input_id).
  std::array<int, 4> canonical = {0, 1, 2, 3};
  std::sort(canonical.begin(), canonical.end(), [&](int a, int b) {
    for (int k = 0; k < 3; ++k) {
      if (sig[a][k] != sig[b][k]) return sig[a][k] < sig[b][k];
    }
    return ids4[a] < ids4[b];
  });
  out_canonical = canonical;

  // 5. Quantize the 5 smallest sorted distances normalized by the largest.
  double largest = dists[edge_order[5]];
  if (!(largest > 0.0)) {
    // Degenerate (coincident stars): return 0 key.
    return 0ULL;
  }
  std::array<uint64_t, 5> q{};
  for (int k = 0; k < 5; ++k) {
    double r = dists[edge_order[k]] / largest;
    double v = std::round(r * static_cast<double>(kQuantScale));
    if (v < 0.0) v = 0.0;
    if (v > static_cast<double>(kQuantScale))
      v = static_cast<double>(kQuantScale);
    q[k] = static_cast<uint64_t>(v);
  }

  return (q[4] << 40) | (q[3] << 30) | (q[2] << 20) | (q[1] << 10) | q[0];
}

// Compute the canonical (key, centroid_canonical_order) for a 4-star pattern
// GIVEN a specific edge ordering. Returns {0, ...} on degenerate input.
// Mirrors steps 3–6 of pattern_key_canonical; the only difference is that
// the caller supplies edge_order rather than the sort-determined ordering.
static std::pair<uint64_t, std::array<int, 4>>
pattern_key_with_edge_order(const std::array<double, 6> &dists,
                            const std::array<int, 6> &edge_order,
                            const std::array<int, 4> &ids4) {
  std::array<std::array<int, 3>, 4> sig{};
  std::array<int, 4> sig_fill = {0, 0, 0, 0};
  for (int rank = 0; rank < 6; ++rank) {
    int e = edge_order[rank];
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    sig[i][sig_fill[i]++] = rank;
    sig[j][sig_fill[j]++] = rank;
  }
  for (int s = 0; s < 4; ++s)
    std::sort(sig[s].begin(), sig[s].end());

  std::array<int, 4> canonical = {0, 1, 2, 3};
  std::sort(canonical.begin(), canonical.end(), [&](int a, int b) {
    for (int k = 0; k < 3; ++k) {
      if (sig[a][k] != sig[b][k]) return sig[a][k] < sig[b][k];
    }
    return ids4[a] < ids4[b];
  });

  double largest = dists[edge_order[5]];
  if (!(largest > 0.0)) return {0ULL, canonical};

  std::array<uint64_t, 5> q{};
  for (int k = 0; k < 5; ++k) {
    double r = dists[edge_order[k]] / largest;
    double v = std::round(r * static_cast<double>(kQuantScale));
    if (v < 0.0) v = 0.0;
    if (v > static_cast<double>(kQuantScale))
      v = static_cast<double>(kQuantScale);
    q[k] = static_cast<uint64_t>(v);
  }
  uint64_t key =
      (q[4] << 40) | (q[3] << 30) | (q[2] << 20) | (q[1] << 10) | q[0];
  return {key, canonical};
}

std::vector<std::pair<uint64_t, std::array<int, 4>>>
pattern_keys_noise_robust(const std::array<std::array<double, 3>, 4> &vecs4,
                          const std::array<int, 4> &ids4, double noise_tol) {
  // Compute the 6 pairwise distances (matches step 1 of pattern_key_canonical).
  std::array<double, 6> dists{};
  for (int e = 0; e < 6; ++e) {
    int i = kEdgePairs[e].first;
    int j = kEdgePairs[e].second;
    double d = vecs4[i][0] * vecs4[j][0] + vecs4[i][1] * vecs4[j][1] +
               vecs4[i][2] * vecs4[j][2];
    if (d > 1.0) d = 1.0;
    if (d < -1.0) d = -1.0;
    dists[e] = std::acos(d);
  }

  // Build the base sorted edge order (same tie-break as canonical).
  std::array<int, 6> base_order = {0, 1, 2, 3, 4, 5};
  std::sort(base_order.begin(), base_order.end(), [&](int a, int b) {
    if (dists[a] != dists[b]) return dists[a] < dists[b];
    if (kEdgePairs[a].first != kEdgePairs[b].first)
      return kEdgePairs[a].first < kEdgePairs[b].first;
    return kEdgePairs[a].second < kEdgePairs[b].second;
  });

  // Identify uncertain adjacent pairs (ranks i, i+1) whose distance gap is
  // within noise_tol. At most ~3 in practice for a noisy 4-star scene.
  std::vector<int> uncertain_pos;
  uncertain_pos.reserve(5);
  for (int i = 0; i < 5; ++i) {
    double gap = dists[base_order[i + 1]] - dists[base_order[i]];
    if (gap < noise_tol) uncertain_pos.push_back(i);
  }

  // Enumerate 2^k orderings (k = uncertain count). For each subset of
  // uncertain positions, swap the corresponding adjacent ranks in the order.
  std::vector<std::pair<uint64_t, std::array<int, 4>>> out;
  std::unordered_set<uint64_t> seen;
  const size_t num_subsets = static_cast<size_t>(1) << uncertain_pos.size();
  out.reserve(num_subsets);
  for (size_t mask = 0; mask < num_subsets; ++mask) {
    std::array<int, 6> order = base_order;
    for (size_t bit = 0; bit < uncertain_pos.size(); ++bit) {
      if (mask & (1u << bit)) {
        int pos = uncertain_pos[bit];
        std::swap(order[pos], order[pos + 1]);
      }
    }
    auto kp = pattern_key_with_edge_order(dists, order, ids4);
    if (!seen.insert(kp.first).second) continue;
    out.emplace_back(kp);
  }
  return out;
}
