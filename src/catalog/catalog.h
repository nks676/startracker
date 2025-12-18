#pragma once

#include <vector>
#include <string>

struct Star {
    int id;
    double x;
    double y;
    double z;
    double magnitude;
};

std::vector<Star> csv_to_catalog(const std::string& filename);