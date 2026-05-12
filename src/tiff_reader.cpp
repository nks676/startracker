#include "tiff_reader.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace {

// TIFF tag IDs we care about. Anything else is silently skipped.
constexpr uint16_t kTagImageWidth = 256;
constexpr uint16_t kTagImageLength = 257;
constexpr uint16_t kTagBitsPerSample = 258;
constexpr uint16_t kTagCompression = 259;
constexpr uint16_t kTagPhotometricInterp = 262;
constexpr uint16_t kTagStripOffsets = 273;
constexpr uint16_t kTagSamplesPerPixel = 277;
constexpr uint16_t kTagStripByteCounts = 279;
constexpr uint16_t kTagSampleFormat = 339;

// TIFF data types we need. Each carries an in-IFD-entry size in bytes; values
// that fit in 4 bytes are stored inline in the entry's value field, otherwise
// the field is an offset to the actual data. For our purposes BYTE, SHORT,
// LONG cover every tag we read.
constexpr uint16_t kTypeByte = 1;   // 1 byte
constexpr uint16_t kTypeShort = 3;  // 2 bytes
constexpr uint16_t kTypeLong = 4;   // 4 bytes

inline uint16_t read_u16(const uint8_t *p, bool little) {
  return little ? static_cast<uint16_t>(p[0] | (p[1] << 8))
                : static_cast<uint16_t>((p[0] << 8) | p[1]);
}
inline uint32_t read_u32(const uint8_t *p, bool little) {
  if (little) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
  }
  return (static_cast<uint32_t>(p[0]) << 24) |
         (static_cast<uint32_t>(p[1]) << 16) |
         (static_cast<uint32_t>(p[2]) << 8) | static_cast<uint32_t>(p[3]);
}

// Read the inline value of an IFD entry as an unsigned integer, accepting
// BYTE/SHORT/LONG. The entry value field is 4 bytes; for types smaller than
// 4 bytes the value occupies the leading bytes (endian-corrected).
uint32_t read_entry_inline(const uint8_t *value_bytes, uint16_t typ,
                           bool little) {
  switch (typ) {
  case kTypeByte:
    return value_bytes[0];
  case kTypeShort:
    return read_u16(value_bytes, little);
  case kTypeLong:
    return read_u32(value_bytes, little);
  default:
    throw std::runtime_error("TIFF: unexpected IFD entry type for inline read");
  }
}

} // namespace

void read_tiff16_grayscale(const std::string &path,
                           std::vector<uint16_t> &out, int &width,
                           int &height) {
  std::ifstream fs(path, std::ios::binary);
  if (!fs)
    throw std::runtime_error("TIFF: cannot open " + path);

  // Slurp the whole header region in one go. Most TIFFs put the IFD right
  // after the 8-byte file header; for the ESA fixtures and Pi HQ Cam output
  // the IFD is well under 4 KB, so 64 KB is generous. We re-seek for the
  // pixel payload below.
  std::vector<uint8_t> header(65536);
  fs.read(reinterpret_cast<char *>(header.data()), header.size());
  size_t header_bytes = static_cast<size_t>(fs.gcount());
  if (header_bytes < 8)
    throw std::runtime_error("TIFF: file too short for header");

  // Byte-order: 'II' = little-endian (Intel), 'MM' = big-endian (Motorola).
  bool little;
  if (header[0] == 'I' && header[1] == 'I')
    little = true;
  else if (header[0] == 'M' && header[1] == 'M')
    little = false;
  else
    throw std::runtime_error("TIFF: bad byte-order marker");

  uint16_t magic = read_u16(header.data() + 2, little);
  if (magic != 42)
    throw std::runtime_error("TIFF: bad magic (expected 42, classic TIFF)");

  uint32_t ifd_off = read_u32(header.data() + 4, little);
  if (ifd_off + 2 > header_bytes)
    throw std::runtime_error("TIFF: IFD offset beyond cached header region");

  uint16_t n_entries = read_u16(header.data() + ifd_off, little);
  size_t entries_start = ifd_off + 2;
  if (entries_start + static_cast<size_t>(n_entries) * 12 > header_bytes)
    throw std::runtime_error("TIFF: IFD entries beyond cached header region");

  // Defaults per TIFF 6.0 spec for tags we tolerate as omitted.
  uint32_t f_width = 0, f_height = 0;
  uint32_t f_bps = 1;             // BitsPerSample default 1
  uint32_t f_compression = 1;     // None
  uint32_t f_photometric = 0;     // 0 = WhiteIsZero (we want 1 = BlackIsZero)
  bool have_photometric = false;
  uint32_t f_strip_offset = 0;
  uint32_t f_strip_count = 0;
  uint32_t f_samples_per_pixel = 1;
  uint32_t f_sample_format = 1; // 1 = unsigned int (default)

  // Multi-strip support is intentionally minimal: if a TIFF lists more than
  // one strip, we reject it. ESA fixtures and Pi HQ Cam use single-strip.
  bool multi_strip = false;

  for (int i = 0; i < n_entries; ++i) {
    const uint8_t *e = header.data() + entries_start + i * 12;
    uint16_t tag = read_u16(e, little);
    uint16_t typ = read_u16(e + 2, little);
    uint32_t cnt = read_u32(e + 4, little);
    const uint8_t *val = e + 8; // 4-byte value or offset field

    switch (tag) {
    case kTagImageWidth:
      f_width = read_entry_inline(val, typ, little);
      break;
    case kTagImageLength:
      f_height = read_entry_inline(val, typ, little);
      break;
    case kTagBitsPerSample:
      if (cnt != 1) {
        // For multi-channel images BitsPerSample is an array; we only handle
        // single-channel grayscale, which has count=1.
        throw std::runtime_error(
            "TIFF: BitsPerSample with count != 1 (not grayscale?)");
      }
      f_bps = read_entry_inline(val, typ, little);
      break;
    case kTagCompression:
      f_compression = read_entry_inline(val, typ, little);
      break;
    case kTagPhotometricInterp:
      f_photometric = read_entry_inline(val, typ, little);
      have_photometric = true;
      break;
    case kTagStripOffsets:
      if (cnt != 1) {
        multi_strip = true;
      } else {
        f_strip_offset = read_entry_inline(val, typ, little);
      }
      break;
    case kTagSamplesPerPixel:
      f_samples_per_pixel = read_entry_inline(val, typ, little);
      break;
    case kTagStripByteCounts:
      if (cnt != 1) {
        multi_strip = true;
      } else {
        f_strip_count = read_entry_inline(val, typ, little);
      }
      break;
    case kTagSampleFormat:
      f_sample_format = read_entry_inline(val, typ, little);
      break;
    default:
      break;
    }
  }

  if (f_width == 0 || f_height == 0)
    throw std::runtime_error("TIFF: missing ImageWidth / ImageLength");
  if (f_bps != 16)
    throw std::runtime_error("TIFF: BitsPerSample != 16 (this reader only "
                             "handles 16-bit grayscale)");
  if (f_compression != 1)
    throw std::runtime_error("TIFF: compressed (Compression != 1) unsupported");
  if (f_samples_per_pixel != 1)
    throw std::runtime_error(
        "TIFF: SamplesPerPixel != 1 (only single-channel grayscale)");
  if (multi_strip)
    throw std::runtime_error(
        "TIFF: multi-strip image unsupported (single strip only)");
  if (f_strip_offset == 0 || f_strip_count == 0)
    throw std::runtime_error("TIFF: missing StripOffsets / StripByteCounts");
  if (f_sample_format != 1) // 1 = unsigned int; anything else (signed, float)
    throw std::runtime_error("TIFF: SampleFormat != unsigned int unsupported");
  // PhotometricInterp 0 (WhiteIsZero) is technically valid but inverts the
  // intensity convention. Reject it loudly rather than silently producing an
  // inverted image — Pi HQ Cam emits 1 (BlackIsZero), ESA fixtures do too.
  if (have_photometric && f_photometric != 1)
    throw std::runtime_error(
        "TIFF: PhotometricInterpretation must be 1 (BlackIsZero)");

  const size_t expected_bytes =
      static_cast<size_t>(f_width) * f_height * sizeof(uint16_t);
  if (f_strip_count != expected_bytes)
    throw std::runtime_error("TIFF: StripByteCounts mismatch with WxHx2");

  // Seek to the pixel payload and read raw uint16 raster.
  fs.clear();
  fs.seekg(f_strip_offset, std::ios::beg);
  if (!fs)
    throw std::runtime_error("TIFF: seek to pixel data failed");

  out.assign(static_cast<size_t>(f_width) * f_height, 0);
  fs.read(reinterpret_cast<char *>(out.data()),
          static_cast<std::streamsize>(expected_bytes));
  if (static_cast<size_t>(fs.gcount()) != expected_bytes)
    throw std::runtime_error("TIFF: pixel data short read");

  // Endian-swap if file is big-endian and host is little-endian (or vice
  // versa). We byte-swap unconditionally when the file's byte order doesn't
  // match the host; on macOS/Linux x86-64 and ARM64 the host is little-endian,
  // so only big-endian TIFFs need swapping in practice.
  bool host_little = []() {
    uint16_t v = 1;
    return *reinterpret_cast<const uint8_t *>(&v) == 1;
  }();
  if (little != host_little) {
    for (auto &v : out) {
      v = static_cast<uint16_t>((v >> 8) | (v << 8));
    }
  }

  width = static_cast<int>(f_width);
  height = static_cast<int>(f_height);
}
