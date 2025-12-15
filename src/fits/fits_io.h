#pragma once

#include <vector>
#include <string>

struct ImageData
{
    std::vector<double> pixels;
    long width;
    long height;
};

ImageData fits_to_data(const std::string& filename);
