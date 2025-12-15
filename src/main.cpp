#include <iostream>
#include "fits/fits_io.h"

int main() {
    std::string filename = "data/m42_40min_red.fits";

    ImageData orion = fits_to_data(filename);

    std::cout << "Image Width: " << orion.width << ", Height: " << orion.height << std::endl;

    return 0;
}