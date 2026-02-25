#pragma once
#include "identification.h"
#include <vector>

struct Quaternion {
  double x, y, z, w;
};

// Estimates the quaternion (Inertial to Camera) using TRIAD.
// Requires at least 2 identified stars.
Quaternion estimate_attitude(const std::vector<IdentifiedStar> &stars,
                             const StarDatabase &db);
