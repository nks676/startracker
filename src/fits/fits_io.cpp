#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>

#include "fits_io.h"
#include "fitsio.h" // CFITSIO

ImageData fits_to_data(const std::string& filename) {
    ImageData data = {};

    fitsfile *fptr;
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[2] = {1, 1};
    
    if (fits_open_file(&fptr, filename.c_str(), READONLY, &status) ) {
        fits_report_error(stderr, status);
        return data;
    }

    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return data;
    }

    if (naxis != 2) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return data;
    }

    data.width = naxes[0];
    data.height = naxes[1];
    const long num_pixels = data.width * data.height;
    data.pixels.resize(num_pixels); // pixel array
    long fpixel[2] = {1, 1};
    int anynul = 0;

    if (fits_read_pix(fptr, TDOUBLE, fpixel, num_pixels, nullptr, data.pixels.data(), &anynul, &status)){
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return data;
    }

    // mean, stddev, max
    double sum = 0.0;
    double sum_sq = 0.0;
    double max = 0.0;
    for (double val: data.pixels) {
        sum += val;
        sum_sq += val * val;
    
        if (val > max) {
            max = val;
        }
    }
    data.intensity_mean = sum / num_pixels;
    data.intensity_standard_deviation = sqrt((sum_sq / num_pixels) - (data.intensity_mean * data.intensity_mean));
    
    // threshold
    data.intensity_threshold = data.intensity_mean + THRESHOLD_CONSTANT * data.intensity_standard_deviation;
    if (data.intensity_threshold >= max) {
        std::cout << "Warning: Computed threshold exceeds max intensity." << std::endl;
        data.intensity_threshold = data.intensity_mean + 0.8 * (max - data.intensity_mean);
    }

    // mask
    data.pixels_mask.resize(num_pixels);
    for (long i = 0; i < num_pixels; i++) {
        data.pixels_mask[i] = (data.pixels[i] >= data.intensity_threshold);
    }

    // labeling using Union-Find
    std::vector<int> labels(num_pixels, 0);
    UnionFind uf(num_pixels / 2 + 1);
    int label = 1;
    for (long y = 0; y < data.height; ++y) {
        for (long x = 0; x < data.width; ++x) {
            long idx = y * data.width + x;
            if (!data.pixels_mask[idx]) {
                continue;
            }

            int left = (x > 0) ? labels[idx - 1] : 0;
            int top = (y > 0) ? labels[idx - data.width] : 0;

            if (left == 0 && top == 0) {
                labels[idx] = label;
                label++;
            } else if (left != 0 && top == 0) {
                labels[idx] = left;
            } else if (left == 0 && top != 0) {
                labels[idx] = top;
            } else {
                labels[idx] = left;
                if (left != top) 
                {
                    uf.unite(left, top);
                }
            }
        }
    }

    // resolve and cluster
    std::map<int, int> label_to_vector;

    for (long i = 0; i < num_pixels; ++i) {
        if (labels[i] != 0) {
            int root = uf.find(labels[i]);
            labels[i] = root;

            // create new cluster if it doesn't exist
            if (label_to_vector.find(root) == label_to_vector.end()) {
                label_to_vector[root] = data.clusters.size();
                data.clusters.push_back(Cluster{root, {}});
            }

            int cluster_idx = label_to_vector[root];

            long x = i % data.width;
            long y = i / data.width;
            data.clusters[cluster_idx].pixels.push_back(Pixel{x, y, data.pixels[i]});
        }
    }

    // compute centroids and total intensities
    for (Cluster& cluster : data.clusters) {
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_intensity = 0.0;

        for (const Pixel& pixel : cluster.pixels) {
            sum_x += pixel.x * pixel.intensity;
            sum_y += pixel.y * pixel.intensity;
            sum_intensity += pixel.intensity;
        }

        if (sum_intensity > 0) {
            cluster.x_centroid = sum_x / sum_intensity;
            cluster.y_centroid = sum_y / sum_intensity;
            cluster.total_intensity = sum_intensity;
        }
    }

    // keep top 50 brightest clusters
    std::sort(data.clusters.begin(), data.clusters.end(), [](const Cluster& a, const Cluster& b) {
        return a.total_intensity > b.total_intensity;
    });
    if (data.clusters.size() > 50) {
        data.clusters.resize(50);
    }

    // Debug: Print the top 3 to see if they look real
    for(int i=0; i<3 && i < (int)data.clusters.size(); i++) {
        std::cout << "Star " << i << ": Pos(" 
                << data.clusters[i].x_centroid << ", " 
                << data.clusters[i].y_centroid << ") Intensity: " 
                << data.clusters[i].total_intensity << std::endl;
    }

    fits_close_file(fptr, &status);

    return data;
}