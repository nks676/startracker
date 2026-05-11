#pragma once
#include <cstdint>
#include <vector>

struct StarCentroid {
  double x;
  double y;
  double intensity;
};

// Shape/size filter parameters applied to connected components by both the
// fixed-threshold and adaptive extractors. Defaults are tuned for typical
// 1024x1024 star tracker images: reject single-pixel hot pixels, oversized
// blobs (satellites, scattered light), streaks, and edge artifacts.
struct CentroidFilterParams {
  int min_pixels = 2;       // > 1 -> reject single hot pixels
  int max_pixels = 200;     // reject extended objects / trails
  double max_aspect = 3.0;  // bbox long/short axis ratio
  int border_margin = 5;    // bbox must be this far from any edge
};

// Subtracts a tile-based median background estimate from an 8-bit image.
// Returns a new buffer with same dimensions. Pixels are clamped to [0,255].
// tile_size: side length of each tile (default 64). Background between tiles
// is bilinearly interpolated from per-tile medians.
std::vector<uint8_t> subtract_background(const uint8_t *image, int width,
                                         int height, int tile_size = 64);

// Extracts centroids from an 8-bit grayscale image using global thresholding
// and CoG. Also applies the shape/size filters in CentroidFilterParams.
std::vector<StarCentroid>
extract_centroids(const uint8_t *image_data, int width, int height,
                  uint8_t threshold,
                  const CentroidFilterParams &filter = CentroidFilterParams{});

// Connected-component centroid extraction with per-tile adaptive thresholding.
// Each tile (tile_size x tile_size) computes its own local mean+stddev; a
// pixel is "above threshold" if it exceeds local_mean + k_sigma * local_stddev.
// Connected components use the same 8-connected BFS as extract_centroids and
// the same CentroidFilterParams shape/size filters are applied.
std::vector<StarCentroid> extract_centroids_adaptive(
    const uint8_t *image, int width, int height, double k_sigma = 5.0,
    int tile_size = 64,
    const CentroidFilterParams &filter = CentroidFilterParams{});
