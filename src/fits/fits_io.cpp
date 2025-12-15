#include <iostream>
#include "fits_io.h"
#include "fitsio.h" // CFITSIO

ImageData fits_to_data(const std::string& filename) {
    ImageData data = {};

    fitsfile *fptr;
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[2] = {1, 1};
    
    if (fits_open_file(&fptr, filename.c_str(), READONLY, &status) )
    {
        fits_report_error(stderr, status);
        return data;
    }

    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status))
    {
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

    if (fits_read_pix(fptr, TDOUBLE, fpixel, num_pixels, nullptr, data.pixels.data(), &anynul, &status))
    {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return data;
    }

    fits_close_file(fptr, &status);

    return data;
}