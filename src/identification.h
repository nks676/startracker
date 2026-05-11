#pragma once
#include "catalog.h"
#include "image_processing.h"
#include <vector>

// Brown-Conrady camera model. Defaults for k1..p2 give pure pinhole behavior,
// so existing callers that only set focal/center fields continue to work.
struct CameraModel {
  double frame_width;
  double frame_height;
  double focal_x;
  double focal_y;
  double center_x;
  double center_y;
  // Brown-Conrady distortion (defaults: zero distortion = pure pinhole)
  double k1 = 0.0; // radial
  double k2 = 0.0;
  double k3 = 0.0;
  double p1 = 0.0; // tangential
  double p2 = 0.0;
};

// Back-compat alias — existing tests and code referring to PinholeCamera still
// work and default-construct to zero distortion.
using PinholeCamera = CameraModel;

// Forward-project a camera-frame unit (or non-unit) vector to pixel
// coordinates, applying Brown-Conrady distortion. Caller is responsible for
// ensuring v[2] > 0 (the point is in front of the camera).
void project(const CameraModel &cam, const double v[3], double &px,
             double &py);

// Inverse map: undistort a pixel coordinate back to a camera-frame unit
// vector. Uses fixed-point iteration (typically converges in ~10 iterations
// to better than 1e-10). If all distortion coefficients are zero this is
// equivalent to the textbook pinhole back-projection.
void undistort_to_unit_vector(const CameraModel &cam, double px, double py,
                              double v_out[3]);

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
