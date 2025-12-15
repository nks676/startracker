#include <iostream>
#include <vector>
#include <string>

#include "fitsio.h"

std::vector<double> fits_to_vector(const std::string& filename, long& width, long& height) {
    fitsfile *fptr;
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[2] = {1, 1};

    if (fits_open_file(&fptr, filename.c_str(), READONLY, &status) )
    {
        fits_report_error(stderr, status);
        return {};
    }

    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status))
    {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return {};
    }

    if (naxis != 2) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return {};
    }

    width = naxes[0];
    height = naxes[1];
    const long num_pixels = width * height;
    std::vector<double> pixels(num_pixels); // pixel array
    long fpixel[2] = {1, 1};
    int anynul = 0;

    if (fits_read_pix(fptr, TDOUBLE, fpixel, num_pixels, nullptr, pixels.data(), &anynul, &status))
    {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return {};
    }

    fits_close_file(fptr, &status);

    return pixels;
}

int main() {
    std::string filename = "m42_40min_red.fits";
    long width = 0;
    long height = 0;

    std::vector<double> image_data = fits_to_vector(filename, width, height);

    std::cout << "Image Width: " << width << ", Height: " << height << std::endl;

    return 0;
}