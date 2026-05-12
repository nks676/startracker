#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

#include "catalog.h"
#include "estimation.h"
#include "identification.h"
#include "image_processing.h"
#include "tiff_reader.h"

int main(int argc, char **argv) {
  // Peel off the optional `--benchmark` flag from any positional slot before
  // the legacy positional argv counting kicks in. Keeps default output byte-
  // identical when the flag is absent.
  bool benchmark = false;
  {
    std::vector<char *> filtered;
    filtered.reserve(argc);
    for (int i = 0; i < argc; ++i) {
      if (std::string(argv[i]) == "--benchmark") {
        benchmark = true;
      } else {
        filtered.push_back(argv[i]);
      }
    }
    argc = static_cast<int>(filtered.size());
    for (int i = 0; i < argc; ++i) {
      argv[i] = filtered[i];
    }
  }

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

  // Load Image. .tif/.tiff route through the native 16-bit reader; everything
  // else (PNG, JPEG, BMP) goes through stb_image at 8-bit. Detection is by
  // file extension — simple and matches the regression fixture naming.
  auto ends_with_ci = [](const std::string &s, const std::string &suffix) {
    if (suffix.size() > s.size()) return false;
    for (size_t i = 0; i < suffix.size(); ++i) {
      char a = s[s.size() - suffix.size() + i];
      char b = suffix[i];
      if (a >= 'A' && a <= 'Z') a = static_cast<char>(a + ('a' - 'A'));
      if (a != b) return false;
    }
    return true;
  };
  const bool is_tiff =
      ends_with_ci(img_path, ".tif") || ends_with_ci(img_path, ".tiff");

  int width = 0, height = 0;
  std::vector<uint16_t> image16;
  uint8_t *image_data = nullptr;
  int channels = 0;
  if (is_tiff) {
    try {
      read_tiff16_grayscale(img_path, image16, width, height);
    } catch (const std::exception &e) {
      std::cerr << "Failed to load TIFF (" << img_path << "): " << e.what()
                << "\n";
      return 1;
    }
  } else {
    image_data = stbi_load(img_path.c_str(), &width, &height, &channels, 1);
    if (!image_data) {
      std::cerr << "Failed to load image: " << img_path << "\n";
      return 1;
    }
  }

  std::cout << "Loaded image: " << width << "x" << height
            << (is_tiff ? " (16-bit TIFF)" : "") << "\n";

  using clk = std::chrono::high_resolution_clock;
  auto t_pipeline_start = clk::now();

  // 1. Image Processing — adaptive per-tile thresholding handles varying
  // backgrounds across the frame (vignetting, light pollution, sensor glow).
  auto t_centroid_start = clk::now();
  auto centroids = is_tiff ? extract_centroids_adaptive_gaussian(
                                 image16.data(), width, height)
                           : extract_centroids_adaptive_gaussian(
                                 image_data, width, height);
  auto t_centroid_end = clk::now();
  std::cout << "Extracted " << centroids.size() << " centroids.\n";
  if (!is_tiff)
    stbi_image_free(image_data);
  if (benchmark) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                  t_centroid_end - t_centroid_start)
                  .count();
    std::cout << "[bench] stage=centroid us=" << us << "\n";
  }

  // Real cameras commonly see fainter stars than fit in the catalog
  // (catalog cutoff Vmag<7, sensors reach ~Vmag 8-9). Keep only the
  // brightest CENTROID_CAP — non-catalog faint stars can't be identified
  // anyway, and the pyramid's combinatorics scale ~ N^2 per seed-pair check.
  // Rank by peak (max pixel value in the component) rather than intensity
  // (sum): 8-bit saturated blobs inflate sum but cap peak at 255, so peak
  // gives a more honest brightness ordering when stars saturate.
  constexpr size_t CENTROID_CAP = 50;
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
  auto t_catalog_start = clk::now();
  StarDatabase db(star_path, pair_path);

  // Phase 3e: load the pattern-hash catalog for the FOV bin closest to the
  // camera's FOV. Pattern catalogs are generated for fixed bins of 10°, 15°,
  // and 20°; we round-robin into the nearest. Missing file → fall through to
  // pyramid path (load_pattern_catalog throws, which we swallow).
  {
    int fov_bin = (fov_deg <= 12.5) ? 10 : (fov_deg <= 17.5) ? 15 : 20;
    std::string pattern_dir;
    auto slash = pair_path.find_last_of("/\\");
    if (slash != std::string::npos) {
      pattern_dir = pair_path.substr(0, slash + 1);
    }
    std::string pattern_path =
        pattern_dir + "catalog_patterns_" + std::to_string(fov_bin) + ".bin";
    try {
      db.load_pattern_catalog(pattern_path);
    } catch (const std::exception &e) {
      std::cerr << "Pattern catalog not loaded (" << e.what()
                << "); identification will use pyramid fallback.\n";
    }
  }
  auto t_catalog_end = clk::now();
  if (benchmark) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                  t_catalog_end - t_catalog_start)
                  .count();
    std::cout << "[bench] stage=catalog_load us=" << us << "\n";
  }

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

  auto t_identify_start = clk::now();
  auto identified = identify_stars(centroids, camera, db, cos_tol);
  auto t_identify_end = clk::now();
  std::cout << "Identified " << identified.size() << " stars.\n";
  if (benchmark) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                  t_identify_end - t_identify_start)
                  .count();
    std::cout << "[bench] stage=identify us=" << us << "\n";
  }

  if (identified.size() < 2) {
    std::cerr << "Not enough stars identified for attitude estimation.\n";
    return 1;
  }

  // 4. Attitude Estimation
  try {
    auto t_estimate_start = clk::now();
    Quaternion q = estimate_attitude(identified, db);
    auto t_estimate_end = clk::now();
    std::cout << std::setprecision(17) << "Estimated Quaternion: [" << q.x
              << ", " << q.y << ", " << q.z << ", " << q.w << "]\n";
    if (benchmark) {
      auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                    t_estimate_end - t_estimate_start)
                    .count();
      std::cout << "[bench] stage=estimate us=" << us << "\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "Estimation failed: " << e.what() << "\n";
  }

  if (benchmark) {
    auto t_pipeline_end = clk::now();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                  t_pipeline_end - t_pipeline_start)
                  .count();
    std::cout << "[bench] stage=total us=" << us << "\n";
  }

  return 0;
}
