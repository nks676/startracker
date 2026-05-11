#include <algorithm>
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
  // Accepted argc: 4 (mandatory), 5 (+fov), 6 (+cos_tol), 11 (+all five
  // distortion coefficients). 7..10 means a partial distortion spec, which
  // is ambiguous so we reject it up front.
  if (argc < 4 || (argc > 6 && argc < 11) || argc > 11) {
    std::cerr
        << "Usage: " << argv[0]
        << " <image.png> <catalog_stars.bin> <catalog_pairs.bin> "
           "[fov_deg] [cos_tol] [k1 k2 p1 p2 k3]\n"
        << "  Distortion coefficients are all-or-nothing: pass all five or "
           "none.\n";
    return 1;
  }

  std::string img_path = argv[1];
  std::string star_path = argv[2];
  std::string pair_path = argv[3];
  double fov_deg = (argc >= 5) ? std::stod(argv[4]) : 20.0;
  double cos_tol = (argc >= 6) ? std::stod(argv[5]) : 1e-5;
  double k1 = (argc == 11) ? std::stod(argv[6]) : 0.0;
  double k2 = (argc == 11) ? std::stod(argv[7]) : 0.0;
  double p1 = (argc == 11) ? std::stod(argv[8]) : 0.0;
  double p2 = (argc == 11) ? std::stod(argv[9]) : 0.0;
  double k3 = (argc == 11) ? std::stod(argv[10]) : 0.0;

  // Load Image
  int width, height, channels;
  uint8_t *image_data =
      stbi_load(img_path.c_str(), &width, &height, &channels, 1);
  if (!image_data) {
    std::cerr << "Failed to load image: " << img_path << "\n";
    return 1;
  }

  std::cout << "Loaded image: " << width << "x" << height << "\n";

  // 1. Image Processing — adaptive per-tile thresholding handles varying
  // backgrounds across the frame (vignetting, light pollution, sensor glow).
  auto centroids =
      extract_centroids_adaptive_gaussian(image_data, width, height);
  std::cout << "Extracted " << centroids.size() << " centroids.\n";
  stbi_image_free(image_data);

  // Real cameras commonly see fainter stars than fit in the catalog
  // (catalog cutoff Vmag<7, sensors reach ~Vmag 8-9). Keep only the
  // brightest CENTROID_CAP — non-catalog faint stars can't be identified
  // anyway, and the pyramid's combinatorics scale ~ N^2 per seed-pair check.
  // Rank by peak (max pixel value in the component) rather than intensity
  // (sum): 8-bit saturated blobs inflate sum but cap peak at 255, so peak
  // gives a more honest brightness ordering when stars saturate.
  constexpr size_t CENTROID_CAP = 60;
  if (centroids.size() > CENTROID_CAP) {
    std::partial_sort(centroids.begin(), centroids.begin() + CENTROID_CAP,
                      centroids.end(),
                      [](const StarCentroid &a, const StarCentroid &b) {
                        return a.peak > b.peak;
                      });
    centroids.resize(CENTROID_CAP);
    std::cout << "Kept top " << CENTROID_CAP << " by peak intensity.\n";
  }

  // 2. Load Catalog
  StarDatabase db(star_path, pair_path);

  // 3. Star Identification
  CameraModel camera;
  camera.frame_width = width;
  camera.frame_height = height;
  camera.center_x = width / 2.0;
  camera.center_y = height / 2.0;

  camera.focal_x = width / (2.0 * std::tan(fov_deg * M_PI / 180.0 / 2.0));
  camera.focal_y = camera.focal_x; // square pixels

  camera.k1 = k1;
  camera.k2 = k2;
  camera.k3 = k3;
  camera.p1 = p1;
  camera.p2 = p2;

  auto identified = identify_stars(centroids, camera, db, cos_tol);
  std::cout << "Identified " << identified.size() << " stars.\n";

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
