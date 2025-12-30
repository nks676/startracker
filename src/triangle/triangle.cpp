#include <vector>
#include <algorithm>
#include <limits>

#include "triangle.h"

std::vector<Triangle> catalog_to_triangles(const std::vector<Star>& catalog) {
    std::vector<Triangle> triangles;
    int n = catalog.size();

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dot_ab = catalog[i].x*catalog[j].x + catalog[i].y*catalog[j].y + catalog[i].z*catalog[j].z;
            dot_ab = std::clamp(dot_ab, -1.0, 1.0);
            double dist_ab = acos(dot_ab);
            if (dist_ab > MAX_FOV_RAD) continue;

            for (int k = j + 1; k < n; ++k) {
                double dot_ac = catalog[i].x*catalog[k].x + catalog[i].y*catalog[k].y + catalog[i].z*catalog[k].z;
                dot_ac = std::clamp(dot_ac, -1.0, 1.0);
                double dist_ac = acos(dot_ac);
                if (dist_ac > MAX_FOV_RAD) continue;

                double dot_bc = catalog[j].x*catalog[k].x + catalog[j].y*catalog[k].y + catalog[j].z*catalog[k].z;
                dot_bc = std::clamp(dot_bc, -1.0, 1.0);
                double dist_bc = acos(dot_bc);
                if (dist_bc > MAX_FOV_RAD) continue;

                Triangle t;
                t.star1 = catalog[i].id;
                t.star2 = catalog[j].id;
                t.star3 = catalog[k].id;
                
                double sides[3] = {dist_ab, dist_ac, dist_bc};
                std::sort(sides, sides + 3);
                t.a = sides[0];
                t.b = sides[1];
                t.c = sides[2];

                triangles.push_back(t);
            }
        }
    }

    std::sort(triangles.begin(), triangles.end(), [](const Triangle& t1, const Triangle& t2) {
        return t1.a < t2.a;
    });

    return triangles;
}

Triangle find_triangle(const Star& s1, const Star& s2, const Star& s3, const std::vector<Triangle>& triangles) {
    double dot_12 = std::clamp(s1.x*s2.x + s1.y*s2.y + s1.z*s2.z, -1.0, 1.0);
    double dot_23 = std::clamp(s2.x*s3.x + s2.y*s3.y + s2.z*s3.z, -1.0, 1.0);
    double dot_31 = std::clamp(s3.x*s1.x + s3.y*s1.y + s3.z*s1.z, -1.0, 1.0);

    std::vector<double> sides = {
        std::acos(dot_12),
        std::acos(dot_23),
        std::acos(dot_31)
    };
    
    std::sort(sides.begin(), sides.end());
    
    double obs_a = sides[0];
    double obs_b = sides[1];
    double obs_c = sides[2];

    Triangle target; 
    target.a = obs_a - TOLERANCE_RAD;

    auto it = std::lower_bound(triangles.begin(), triangles.end(), target, 
        [](const Triangle& t1, const Triangle& t2) {
            return t1.a < t2.a;
        }
    );

    Triangle best = {-1, -1, -1, 0, 0, 0};
    double best_error = std::numeric_limits<double>::infinity();

    for (; it != triangles.end(); ++it) {
        if (it->a > obs_a + TOLERANCE_RAD) break;

        if (std::abs(it->b - obs_b) < TOLERANCE_RAD && 
            std::abs(it->c - obs_c) < TOLERANCE_RAD) {

            double error = std::abs(it->a - obs_a) + std::abs(it->b - obs_b) + std::abs(it->c - obs_c);
            if (error < best_error) {
                best_error = error;
                best = *it;
            }
        }
    }

    return best;
}