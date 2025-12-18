#include <iostream>
#include "fits/fits_io.h"
#include "catalog/catalog.h"

int main() {
    std::string m42_filename = "data/m42_40min_red.fits";
    ImageData m42 = fits_to_data(m42_filename);
    std::cout << "Loaded " << m42.clusters.size() << " clusters from " << m42_filename << std::endl;

    std::string m13_filename = "data/m13_test.fits";
    ImageData m13 = fits_to_data(m13_filename);
    std::cout << "Loaded " << m13.clusters.size() << " clusters from " << m13_filename << std::endl;

    std::string starmap_filename = "data/hipparcos-voidmain.csv";
    std::vector<Star> catalog = csv_to_catalog(starmap_filename);
    std::cout << "Loaded " << catalog.size() << " stars from " << starmap_filename << std::endl;

    return 0;
}