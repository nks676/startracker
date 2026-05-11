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
template <typename T>
double bilinear_sample(const std::vector<T> &tile_values, int n_tx, int n_ty,
                       int tile_size, int width, int height, double x,
                       double y) {
  // Locate the two tile-center rows/cols straddling (x,y).
  double ix = -1.0, iy = -1.0;
  int tx0 = 0, tx1 = 0, ty0 = 0, ty1 = 0;
  // Find tx0 such that tile_center(tx0) <= x < tile_center(tx0+1).
  // Tile centers are monotonically increasing.
  for (int tx = 0; tx < n_tx; ++tx) {
    double c = tile_center(tx, tile_size, width);
    if (c <= x) {
      tx0 = tx;
      ix = c;
    } else {
      break;
    }
  }
  tx1 = std::min(tx0 + 1, n_tx - 1);
  double ix1 = tile_center(tx1, tile_size, width);

  for (int ty = 0; ty < n_ty; ++ty) {
    double c = tile_center(ty, tile_size, height);
    if (c <= y) {
      ty0 = ty;
      iy = c;
    } else {
      break;
    }
  }
  ty1 = std::min(ty0 + 1, n_ty - 1);
  double iy1 = tile_center(ty1, tile_size, height);

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

// Shared BFS-driven centroid extractor. `is_above_threshold(idx)` decides
// per-pixel membership so the caller can plug in either a constant threshold
// or a precomputed per-pixel adaptive threshold.
template <typename Predicate>
std::vector<StarCentroid> extract_centroids_impl(const uint8_t *image,
                                                 int width, int height,
                                                 Predicate is_above,
                                                 const CentroidFilterParams &f) {
  std::vector<StarCentroid> centroids;
  std::vector<uint8_t> visited(static_cast<size_t>(width) * height, 0);

  static const int dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};
  static const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int idx = y * width + x;
      if (visited[idx] || !is_above(idx))
        continue;

      std::queue<std::pair<int, int>> q;
      q.push({x, y});
      visited[idx] = 1;

      double sum_x = 0, sum_y = 0, sum_i = 0;
      int pixel_count = 0;
      int bx0 = x, by0 = y, bx1 = x, by1 = y;

      while (!q.empty()) {
        auto [cx, cy] = q.front();
        q.pop();
        int cidx = cy * width + cx;
        double intensity = static_cast<double>(image[cidx]);
        sum_x += cx * intensity;
        sum_y += cy * intensity;
        sum_i += intensity;
        pixel_count++;
        bx0 = std::min(bx0, cx);
        by0 = std::min(by0, cy);
        bx1 = std::max(bx1, cx);
        by1 = std::max(by1, cy);

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

      centroids.push_back(
          StarCentroid{sum_x / sum_i, sum_y / sum_i, sum_i});
    }
  }
  return centroids;
}

} // namespace

std::vector<uint8_t> subtract_background(const uint8_t *image, int width,
                                         int height, int tile_size) {
  if (tile_size < 1)
    tile_size = 1;
  int n_tx = tile_count(width, tile_size);
  int n_ty = tile_count(height, tile_size);

  // Per-tile median. nth_element is O(N) per tile; total cost O(W*H).
  std::vector<uint8_t> bg_tiles(static_cast<size_t>(n_tx) * n_ty, 0);
  std::vector<uint8_t> buf;
  buf.reserve(static_cast<size_t>(tile_size) * tile_size);
  for (int ty = 0; ty < n_ty; ++ty) {
    int y0 = ty * tile_size;
    int y1 = std::min(y0 + tile_size, height);
    for (int tx = 0; tx < n_tx; ++tx) {
      int x0 = tx * tile_size;
      int x1 = std::min(x0 + tile_size, width);
      buf.clear();
      for (int y = y0; y < y1; ++y) {
        const uint8_t *row = image + y * width;
        buf.insert(buf.end(), row + x0, row + x1);
      }
      if (buf.empty())
        continue;
      auto mid = buf.begin() + buf.size() / 2;
      std::nth_element(buf.begin(), mid, buf.end());
      bg_tiles[ty * n_tx + tx] = *mid;
    }
  }

  std::vector<uint8_t> out(static_cast<size_t>(width) * height, 0);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double bg = bilinear_sample(bg_tiles, n_tx, n_ty, tile_size, width,
                                  height, x, y);
      int v = static_cast<int>(image[y * width + x]) -
              static_cast<int>(std::lround(bg));
      out[y * width + x] = static_cast<uint8_t>(std::clamp(v, 0, 255));
    }
  }
  return out;
}

std::vector<StarCentroid>
extract_centroids(const uint8_t *image_data, int width, int height,
                  uint8_t threshold, const CentroidFilterParams &filter) {
  auto is_above = [&](int idx) { return image_data[idx] > threshold; };
  return extract_centroids_impl(image_data, width, height, is_above, filter);
}

std::vector<StarCentroid>
extract_centroids_adaptive(const uint8_t *image, int width, int height,
                           double k_sigma, int tile_size,
                           const CentroidFilterParams &filter) {
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
        const uint8_t *row = image + y * width;
        for (int x = x0; x < x1; ++x) {
          double v = row[x];
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

  // Precompute per-pixel threshold. Memory cost is W*H doubles (~8 MB on a
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

  auto is_above = [&](int idx) {
    return static_cast<double>(image[idx]) > thr[idx];
  };
  return extract_centroids_impl(image, width, height, is_above, filter);
}
