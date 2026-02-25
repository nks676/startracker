#pragma once
#include <cstdint>
#include <vector>

struct StarCentroid {
  double x;
  double y;
  double intensity;
};

// Extracts centroids from an 8-bit grayscale image using global thresholding
// and CoG.
std::vector<StarCentroid> extract_centroids(const uint8_t *image_data,
                                            int width, int height,
                                            uint8_t threshold);
