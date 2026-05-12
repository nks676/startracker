#include "camera_model.h"
#include <cmath>

void project(const CameraModel &cam, const double v[3], double &px,
             double &py) {
  // Pinhole step: normalise to z=1 plane.
  const double x = v[0] / v[2];
  const double y = v[1] / v[2];

  const double r2 = x * x + y * y;
  const double radial =
      1.0 + cam.k1 * r2 + cam.k2 * r2 * r2 + cam.k3 * r2 * r2 * r2;
  const double x_d =
      x * radial + 2.0 * cam.p1 * x * y + cam.p2 * (r2 + 2.0 * x * x);
  const double y_d =
      y * radial + cam.p1 * (r2 + 2.0 * y * y) + 2.0 * cam.p2 * x * y;

  px = cam.focal_x * x_d + cam.center_x;
  py = cam.focal_y * y_d + cam.center_y;
}

void undistort_to_unit_vector(const CameraModel &cam, double px, double py,
                              double v_out[3]) {
  // Normalised distorted pixel coordinates.
  const double px_norm = (px - cam.center_x) / cam.focal_x;
  const double py_norm = (py - cam.center_y) / cam.focal_y;

  // Initial guess: assume zero distortion. With k1..p2 = 0 this is the exact
  // answer and the loop body is a no-op (radial == 1, dx = dy = 0), so the
  // result is bitwise identical to the legacy pinhole code path.
  double x = px_norm;
  double y = py_norm;

  for (int iter = 0; iter < 10; ++iter) {
    const double r2 = x * x + y * y;
    const double radial =
        1.0 + cam.k1 * r2 + cam.k2 * r2 * r2 + cam.k3 * r2 * r2 * r2;
    const double dx =
        2.0 * cam.p1 * x * y + cam.p2 * (r2 + 2.0 * x * x);
    const double dy =
        cam.p1 * (r2 + 2.0 * y * y) + 2.0 * cam.p2 * x * y;
    x = (px_norm - dx) / radial;
    y = (py_norm - dy) / radial;
  }

  const double n = std::sqrt(x * x + y * y + 1.0);
  v_out[0] = x / n;
  v_out[1] = y / n;
  v_out[2] = 1.0 / n;
}
