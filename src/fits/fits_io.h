#pragma once

#include <vector>
#include <string>

struct Pixel
{
    long x;
    long y;
    double intensity;
};

struct Cluster
{
    int id;
    std::vector<Pixel> pixels;
};

struct ImageData
{
    std::vector<double> pixels;
    long width;
    long height;

    double intensity_mean;
    double intensity_standard_deviation;

    double intensity_threshold;
    std::vector<bool> pixels_mask; // 1 if star, 0 otherwise

    std::vector<Cluster> clusters;
};

struct UnionFind
{
    std::vector<int> parent;
    
    UnionFind(int size) : parent(size) {
        for (int i = 0; i < size; ++i) {
            parent[i] = i;
        }
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        if (rx != ry) {
            parent[ry] = rx;
        }
    }
};

ImageData fits_to_data(const std::string& filename);
