#include "image_processing.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <vector>

namespace {

// Number of tiles along an axis, rounded up so the last tile covers the
// remaining (possibly partial) strip of pixels.
inline int tile_count(int length, int tile_size) {
  return (length + tile_size - 1) / tile_size;
}

// Tile center in image coordinates. The last tile may be smaller, so we use
// its actual midpoint rather than (tile_index + 0.5) * tile_size to keep
// bilinear interpolation well-behaved at the right/bottom borders.
inline double tile_center(int tile_index, int tile_size, int length) {
  int x0 = tile_index * tile_size;
  int x1 = std::min(x0 + tile_size, length);
  return 0.5 * (x0 + x1 - 1);
}

// Bilinearly interpolate a per-tile scalar field at pixel (x,y). Indices
// outside [0, n_tiles-1] are clamped — pixels near image edges share the
// background estimate of their nearest tile rather than extrapolating.
//
// Phase 3f.1: arithmetic tile-index lookup, was linear scan. For middle
// tiles, tile_center(tx) = tx*tile_size + (tile_size-1)/2, so the largest
// tx with center(tx) <= x is floor((x - (tile_size-1)/2) / tile_size).
// (The last tile may be smaller than tile_size pixels — its center is
// shifted left of the formula's prediction. We handle that by clamping
// the computed tx0 to [0, n_tx-1] and then using the exact tile_center()
// values for the bilinear weights, which absorbs the last-tile offset
// naturally.) Old linear scan cost ~28 iterations per pixel × 786 k pixels
// = ~85 ms on Pi 4 1024×768 — measured as the dominant centroid cost.
template <typename T>
double bilinear_sample(const std::vector<T> &tile_values, int n_tx, int n_ty,
                       int tile_size, int width, int height, double x,
                       double y) {
  const double half_minus = 0.5 * (tile_size - 1);
  int tx0 = static_cast<int>(std::floor((x - half_minus) / tile_size));
  int ty0 = static_cast<int>(std::floor((y - half_minus) / tile_size));
  tx0 = std::clamp(tx0, 0, n_tx - 1);
  ty0 = std::clamp(ty0, 0, n_ty - 1);
  int tx1 = std::min(tx0 + 1, n_tx - 1);
  int ty1 = std::min(ty0 + 1, n_ty - 1);

  // Exact centers (handles the smaller last tile when width % tile_size != 0).
  double ix = tile_center(tx0, tile_size, width);
  double ix1 = tile_center(tx1, tile_size, width);
  double iy = tile_center(ty0, tile_size, height);
  double iy1 = tile_center(ty1, tile_size, height);

  // If the computed tx0 sits to the right of x (pixel is left of the very
  // first tile center, or last-tile center shift overshot), step back one.
  // Guarantees ix <= x for the weight formula below.
  if (ix > x && tx0 > 0) {
    --tx0;
    tx1 = tx0 + 1;
    ix = tile_center(tx0, tile_size, width);
    ix1 = tile_center(tx1, tile_size, width);
  }
  if (iy > y && ty0 > 0) {
    --ty0;
    ty1 = ty0 + 1;
    iy = tile_center(ty0, tile_size, height);
    iy1 = tile_center(ty1, tile_size, height);
  }

  // Compute interpolation weights. Guard against degenerate (single-tile or
  // edge) cases where ix == ix1 or iy == iy1.
  double wx = (ix1 > ix) ? (x - ix) / (ix1 - ix) : 0.0;
  double wy = (iy1 > iy) ? (y - iy) / (iy1 - iy) : 0.0;
  wx = std::clamp(wx, 0.0, 1.0);
  wy = std::clamp(wy, 0.0, 1.0);

  double v00 = static_cast<double>(tile_values[ty0 * n_tx + tx0]);
  double v10 = static_cast<double>(tile_values[ty0 * n_tx + tx1]);
  double v01 = static_cast<double>(tile_values[ty1 * n_tx + tx0]);
  double v11 = static_cast<double>(tile_values[ty1 * n_tx + tx1]);
  double v0 = v00 * (1.0 - wx) + v10 * wx;
  double v1 = v01 * (1.0 - wx) + v11 * wx;
  return v0 * (1.0 - wy) + v1 * wy;
}

// Apply CentroidFilterParams to a candidate component summary. Bounding-box
// dimensions are inclusive (pixel counts), so a 1x10 streak has w=1, h=10.
bool passes_filter(int pixel_count, int bbox_w, int bbox_h, int bbox_x0,
                   int bbox_y0, int bbox_x1, int bbox_y1, int width, int height,
                   const CentroidFilterParams &f) {
  if (pixel_count < f.min_pixels)
    return false;
  if (pixel_count > f.max_pixels)
    return false;
  // Aspect ratio: long axis / short axis. A single-pixel component has w=h=1
  // and trivially passes — but we filter those out via min_pixels=2 anyway.
  int long_axis = std::max(bbox_w, bbox_h);
  int short_axis = std::max(1, std::min(bbox_w, bbox_h));
  if (static_cast<double>(long_axis) / short_axis > f.max_aspect)
    return false;
  // Reject anything whose bounding box touches the border margin. We use the
  // bbox (not the centroid) so that even a small spot whose tail leaks to the
  // edge is rejected; that's where vignetting and decoder artifacts live.
  if (bbox_x0 < f.border_margin || bbox_y0 < f.border_margin)
    return false;
  if (bbox_x1 >= width - f.border_margin ||
      bbox_y1 >= height - f.border_margin)
    return false;
  return true;
}

// Iterative Gaussian-weighted refinement. Initial guess (cx, cy) is the CoG
// already computed from `pixels`. We run a fixed number of iterations of
// intensity * Gaussian-weighted moments, which converges to the true Gaussian
// peak for well-sampled PSFs and is robust to small saturated cores.
//
// sigma = 1.0 px is appropriate for typical star tracker PSFs (~1-1.5 px
// FWHM). The Gaussian weight suppresses pixels far from the current estimate,
// removing the CoG's bias toward asymmetric component pixels (e.g. dim halo
// pixels on one side of the bright core).
struct PixelSample {
  int x;
  int y;
  double i;
};

inline void refine_gaussian(const std::vector<PixelSample> &pixels, double &cx,
                            double &cy) {
  constexpr int kIterations = 4;
  constexpr double kSigma = 1.0;
  const double inv_two_sigma_sq = 1.0 / (2.0 * kSigma * kSigma);
  for (int iter = 0; iter < kIterations; ++iter) {
    double sw_x = 0.0, sw_y = 0.0, sw = 0.0;
    for (const auto &p : pixels) {
      double dx = p.x - cx;
      double dy = p.y - cy;
      double w = std::exp(-(dx * dx + dy * dy) * inv_two_sigma_sq);
      double iw = p.i * w;
      sw_x += p.x * iw;
      sw_y += p.y * iw;
      sw += iw;
    }
    if (sw <= 0.0)
      return; // keep previous estimate; weights collapsed
    cx = sw_x / sw;
    cy = sw_y / sw;
  }
}

// Shared BFS-driven centroid extractor. `is_above_threshold(idx)` decides
// per-pixel membership so the caller can plug in either a constant threshold
// or a precomputed per-pixel adaptive threshold. When `gaussian_refine` is
// true, each component's CoG is refined by intensity * Gaussian-weighted
// moments (see refine_gaussian). `peak` is always populated.
template <typename T, typename Predicate>
std::vector<StarCentroid> extract_centroids_impl_t(const T *image, int width,
                                                    int height,
                                                    Predicate is_above,
                                                    const CentroidFilterParams &f,
                                                    bool gaussian_refine = false) {
  std::vector<StarCentroid> centroids;
  std::vector<uint8_t> visited(static_cast<size_t>(width) * height, 0);

  static const int dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};
  static const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};

  // Reused across components when refinement is on, to avoid reallocating.
  std::vector<PixelSample> pixels;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int idx = y * width + x;
      if (visited[idx] || !is_above(idx))
        continue;

      std::queue<std::pair<int, int>> q;
      q.push({x, y});
      visited[idx] = 1;

      double sum_x = 0, sum_y = 0, sum_i = 0;
      double peak = 0.0;
      int pixel_count = 0;
      int bx0 = x, by0 = y, bx1 = x, by1 = y;
      if (gaussian_refine)
        pixels.clear();

      while (!q.empty()) {
        auto [cx, cy] = q.front();
        q.pop();
        int cidx = cy * width + cx;
        double intensity = static_cast<double>(image[cidx]);
        sum_x += cx * intensity;
        sum_y += cy * intensity;
        sum_i += intensity;
        if (intensity > peak)
          peak = intensity;
        pixel_count++;
        bx0 = std::min(bx0, cx);
        by0 = std::min(by0, cy);
        bx1 = std::max(bx1, cx);
        by1 = std::max(by1, cy);
        if (gaussian_refine)
          pixels.push_back(PixelSample{cx, cy, intensity});

        for (int n = 0; n < 8; ++n) {
          int nx = cx + dx[n];
          int ny = cy + dy[n];
          if (nx < 0 || nx >= width || ny < 0 || ny >= height)
            continue;
          int nidx = ny * width + nx;
          if (!visited[nidx] && is_above(nidx)) {
            visited[nidx] = 1;
            q.push({nx, ny});
          }
        }
      }

      int bbox_w = bx1 - bx0 + 1;
      int bbox_h = by1 - by0 + 1;
      if (sum_i <= 0)
        continue;
      if (!passes_filter(pixel_count, bbox_w, bbox_h, bx0, by0, bx1, by1, width,
                         height, f))
        continue;

      double cx_out = sum_x / sum_i;
      double cy_out = sum_y / sum_i;
      if (gaussian_refine) {
        refine_gaussian(pixels, cx_out, cy_out);
      }

      centroids.push_back(StarCentroid{cx_out, cy_out, sum_i, peak});
    }
  }
  return centroids;
}

} // namespace

namespace {

// Build the per-pixel adaptive threshold buffer used by both the CoG and
// Gaussian-refined adaptive extractors. Factored out of the public function
// bodies so the two share the (non-trivial) preprocessing exactly.
template <typename T>
std::vector<float> build_adaptive_threshold_t(const T *image, int width,
                                              int height, double k_sigma,
                                              int tile_size) {
  if (tile_size < 1)
    tile_size = 1;
  int n_tx = tile_count(width, tile_size);
  int n_ty = tile_count(height, tile_size);

  // Per-tile mean and stddev. We use raw moments for speed (single pass per
  // tile). Stars within a tile bias the mean upward, which is acceptable here
  // because k_sigma=5 puts the threshold well above any bias from a handful
  // of bright pixels in an otherwise-dark tile.
  std::vector<double> mean_tiles(static_cast<size_t>(n_tx) * n_ty, 0.0);
  std::vector<double> std_tiles(static_cast<size_t>(n_tx) * n_ty, 0.0);
  for (int ty = 0; ty < n_ty; ++ty) {
    int y0 = ty * tile_size;
    int y1 = std::min(y0 + tile_size, height);
    for (int tx = 0; tx < n_tx; ++tx) {
      int x0 = tx * tile_size;
      int x1 = std::min(x0 + tile_size, width);
      double sum = 0.0, sum_sq = 0.0;
      int n = 0;
      for (int y = y0; y < y1; ++y) {
        const T *row = image + y * width;
        for (int x = x0; x < x1; ++x) {
          double v = static_cast<double>(row[x]);
          sum += v;
          sum_sq += v * v;
          n++;
        }
      }
      if (n <= 0)
        continue;
      double m = sum / n;
      double var = std::max(0.0, sum_sq / n - m * m);
      // Floor stddev at 1 ADU to avoid div-by-zero / absurdly tight thresholds
      // in pathologically flat tiles (e.g. all zeros). 1 ADU is below the
      // single-bit quantization noise of an 8-bit sensor anyway.
      mean_tiles[ty * n_tx + tx] = m;
      std_tiles[ty * n_tx + tx] = std::max(1.0, std::sqrt(var));
    }
  }

  // Precompute per-pixel threshold. Memory cost is W*H floats (~4 MB on a
  // 1024x1024 image); acceptable, and lets the BFS predicate stay branch-free.
  std::vector<float> thr(static_cast<size_t>(width) * height, 0.f);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double m = bilinear_sample(mean_tiles, n_tx, n_ty, tile_size, width,
                                 height, x, y);
      double s = bilinear_sample(std_tiles, n_tx, n_ty, tile_size, width,
                                 height, x, y);
      thr[y * width + x] = static_cast<float>(m + k_sigma * s);
    }
  }
  return thr;
}

} // namespace

std::vector<StarCentroid>
extract_centroids(const uint8_t *image_data, int width, int height,
                  uint8_t threshold, const CentroidFilterParams &filter) {
  auto is_above = [&](int idx) { return image_data[idx] > threshold; };
  return extract_centroids_impl_t(image_data, width, height, is_above, filter,
                                /*gaussian_refine=*/false);
}

std::vector<StarCentroid>
extract_centroids_gaussian(const uint8_t *image_data, int width, int height,
                           uint8_t threshold,
                           const CentroidFilterParams &filter) {
  auto is_above = [&](int idx) { return image_data[idx] > threshold; };
  return extract_centroids_impl_t(image_data, width, height, is_above, filter,
                                /*gaussian_refine=*/true);
}

std::vector<StarCentroid>
extract_centroids_adaptive(const uint8_t *image, int width, int height,
                           double k_sigma, int tile_size,
                           const CentroidFilterParams &filter) {
  auto thr = build_adaptive_threshold_t(image, width, height, k_sigma, tile_size);
  auto is_above = [&](int idx) {
    return static_cast<double>(image[idx]) > thr[idx];
  };
  return extract_centroids_impl_t(image, width, height, is_above, filter,
                                /*gaussian_refine=*/false);
}

std::vector<StarCentroid> extract_centroids_adaptive_gaussian(
    const uint8_t *image, int width, int height, double k_sigma, int tile_size,
    const CentroidFilterParams &filter) {
  auto thr = build_adaptive_threshold_t(image, width, height, k_sigma, tile_size);
  auto is_above = [&](int idx) {
    return static_cast<double>(image[idx]) > thr[idx];
  };
  return extract_centroids_impl_t(image, width, height, is_above, filter,
                                /*gaussian_refine=*/true);
}

// 16-bit overload — same algorithm, wider pixel type. Pi HQ Camera and the
// ESA tetra3 fixtures emit 12/16-bit raw frames; the 8-bit overload above
// crushes faint stars on those inputs via the 16→8 stretch in
// tools/test_real_images.py. With this overload the binary can consume the
// native TIFF and the adaptive threshold sees the full dynamic range.
std::vector<StarCentroid> extract_centroids_adaptive_gaussian(
    const uint16_t *image, int width, int height, double k_sigma, int tile_size,
    const CentroidFilterParams &filter) {
  auto thr = build_adaptive_threshold_t(image, width, height, k_sigma, tile_size);
  auto is_above = [&](int idx) {
    return static_cast<double>(image[idx]) > thr[idx];
  };
  return extract_centroids_impl_t(image, width, height, is_above, filter,
                                /*gaussian_refine=*/true);
}
