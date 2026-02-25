#pragma once
#include "catalog.h"
#include "image_processing.h"
#include <vector>

struct PinholeCamera {
  double frame_width;
  double frame_height;
  double focal_x;
  double focal_y;
  double center_x;
  double center_y;
};

struct IdentifiedStar {
  int image_idx;
  int catalog_hip_id;
  double v_cam[3]; // Unit vector in camera frame
};

// Identifies stars by voting on pairwise angular distances
std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance);
