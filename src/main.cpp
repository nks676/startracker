#include <cmath>
#include <iostream>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

#include "catalog.h"
#include "estimation.h"
#include "identification.h"
#include "image_processing.h"

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr
        << "Usage: " << argv[0]
        << " <image.png> <catalog_stars.bin> <catalog_pairs.bin> [fov_deg]\n";
    return 1;
  }

  std::string img_path = argv[1];
  std::string star_path = argv[2];
  std::string pair_path = argv[3];
  double fov_deg = (argc >= 5) ? std::stod(argv[4]) : 20.0;

  // Load Image
  int width, height, channels;
  uint8_t *image_data =
      stbi_load(img_path.c_str(), &width, &height, &channels, 1);
  if (!image_data) {
    std::cerr << "Failed to load image: " << img_path << "\n";
    return 1;
  }

  std::cout << "Loaded image: " << width << "x" << height << "\n";

  // 1. Image Processing
  uint8_t threshold = 100; // Background is ~50, 100 should be safe
  auto centroids = extract_centroids(image_data, width, height, threshold);
  std::cout << "Extracted " << centroids.size() << " centroids.\n";
  stbi_image_free(image_data);

  // 2. Load Catalog
  StarDatabase db(star_path, pair_path);

  // 3. Star Identification
  PinholeCamera camera;
  camera.frame_width = width;
  camera.frame_height = height;
  camera.center_x = width / 2.0;
  camera.center_y = height / 2.0;

  camera.focal_x = width / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  camera.focal_y = camera.focal_x; // square pixels

  // Tolerance is purely angular in cosine space
  // 1e-6 cos tol is ~0.08 degrees. Which is ~4 pixels.
  double cos_tol = 1e-6;

  auto identified = identify_stars(centroids, camera, db, cos_tol);
  std::cout << "Identified " << identified.size() << " stars.\n";
  std::cout << "Identified HIPs: ";
  for (const auto &is : identified) {
    std::cout << is.catalog_hip_id << " ";
  }
  std::cout << "\n";

  if (identified.size() < 2) {
    std::cerr << "Not enough stars identified for attitude estimation.\n";
    return 1;
  }

  // 4. Attitude Estimation
  try {
    Quaternion q = estimate_attitude(identified, db);
    std::cout << "Estimated Quaternion: [" << q.x << ", " << q.y << ", " << q.z
              << ", " << q.w << "]\n";
  } catch (const std::exception &e) {
    std::cerr << "Estimation failed: " << e.what() << "\n";
  }

  return 0;
}
