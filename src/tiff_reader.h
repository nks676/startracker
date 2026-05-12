#pragma once
#include <cstdint>
#include <string>
#include <vector>

// Minimal reader for uncompressed 16-bit grayscale baseline TIFF.
//
// Scope: just enough to consume what the Pi HQ Camera (and the ESA tetra3
// fixtures we use as regression data) emits — single-strip, 16-bit,
// little-endian or big-endian, no compression, single sample-per-pixel. This
// deliberately avoids pulling in libtiff: the format we care about is a
// 320-byte header followed by a flat uint16 raster, and that's all this
// function handles. Anything outside that envelope throws.
//
// On success: `out` is resized to width*height uint16_t values in row-major
// top-to-bottom order, and `width`/`height` are filled. Throws
// std::runtime_error on any I/O failure, format mismatch, or unsupported
// option.
void read_tiff16_grayscale(const std::string &path,
                           std::vector<uint16_t> &out, int &width,
                           int &height);
