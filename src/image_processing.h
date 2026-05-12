#pragma once
#include <cstdint>
#include <vector>

struct StarCentroid {
  double x;
  double y;
  double intensity;
  // Peak (max) pixel intensity in the connected component. Used by main.cpp's
  // CENTROID_CAP ranking: peak is robust to saturated blobs that inflate the
  // intensity (sum) metric. Always populated by extract_centroids* APIs.
  double peak = 0.0;
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

// Like extract_centroids but refines each centroid with 3-5 iterations of
// intensity * Gaussian-weighted moments, initialized from the CoG. Sigma is
// fixed at 1.0 px which is appropriate for typical star tracker PSFs.
// Sub-pixel accuracy improves from ~0.1 px (CoG) to ~0.01-0.05 px on
// well-sampled Gaussian PSFs. The `peak` field is populated.
std::vector<StarCentroid> extract_centroids_gaussian(
    const uint8_t *image_data, int width, int height, uint8_t threshold,
    const CentroidFilterParams &filter = CentroidFilterParams{});

// Like extract_centroids_adaptive but with Gaussian-weighted iterative
// refinement after BFS (see extract_centroids_gaussian).
std::vector<StarCentroid> extract_centroids_adaptive_gaussian(
    const uint8_t *image, int width, int height, double k_sigma = 5.0,
    int tile_size = 64,
    const CentroidFilterParams &filter = CentroidFilterParams{});

// 16-bit pixel-type overload. Same algorithm; just consumes the wider raster
// produced by the TIFF reader / Pi HQ Camera path. The float-typed threshold
// buffer is wide enough for the larger dynamic range, so no other tuning is
// needed.
std::vector<StarCentroid> extract_centroids_adaptive_gaussian(
    const uint16_t *image, int width, int height, double k_sigma = 5.0,
    int tile_size = 64,
    const CentroidFilterParams &filter = CentroidFilterParams{});
