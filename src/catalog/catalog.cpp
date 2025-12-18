#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "catalog.h"

std::vector<Star> csv_to_catalog(const std::string& filename) {
    std::vector<Star> catalog;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return catalog;
    }

    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        //  safety check
        if (row.size() < 10) {
            continue;
        }

        try {
            // some stars empty
            if (row[5].empty()) {
                continue;
            }

            // skip faint stars
            double vmag = std::stod(row[5]);
            if (vmag > 6.0) {
                continue; 
            }

            int id = std::stoi(row[1]);
            double ra_rad = std::stod(row[8]) * (M_PI / 180.0);
            double dec_rad = std::stod(row[9]) * (M_PI / 180.0);

            Star s;
            s.id = id;
            s.x = cos(dec_rad) * cos(ra_rad);
            s.y = cos(dec_rad) * sin(ra_rad);
            s.z = sin(dec_rad);
            s.magnitude = vmag;

            catalog.push_back(s);
        } catch(...) {
            continue; // skip bad data
        }
    }

    return catalog;
}
