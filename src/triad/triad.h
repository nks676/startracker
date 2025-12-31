#pragma once

#include "../catalog/catalog.h"

#include <vector>

struct Quaternion {
    double w, x, y, z;
};

struct Observation {
    Star body;
    Star inertial;
    double weight;
};

Quaternion compute_attitude(const std::vector<Observation>& obs);