#pragma once

#include <cmath>

#include "../catalog/catalog.h"

#define MAX_FOV_RAD 10.0 * (M_PI / 180.0) // 10 degrees in radians
#define TOLERANCE_RAD 0.01

struct Triangle {
    int star1, star2, star3;
    double a, b, c; // a <= b <= c
};

std::vector<Triangle> catalog_to_triangles(const std::vector<Star>& catalog);

Triangle find_triangle(const Star& s1, const Star& s2, const Star& s3, const std::vector<Triangle>& triangles);