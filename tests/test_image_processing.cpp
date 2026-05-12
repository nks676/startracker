#include "image_processing.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <gtest/gtest.h>
#include <random>

// --- 2a. Image Processing Tests ---

// Helper: create a blank (black) image buffer
static std::vector<uint8_t> make_blank(int w, int h) {
  return std::vector<uint8_t>(w * h, 0);
}

// Helper: paint a bright spot (3x3 block) centered at (cx, cy)
static void paint_spot(std::vector<uint8_t> &img, int w, int cx, int cy,
                       uint8_t center_val = 200, uint8_t neighbor_val = 150) {
  for (int dy = -1; dy <= 1; ++dy) {
    for (int dx = -1; dx <= 1; ++dx) {
      int x = cx + dx;
      int y = cy + dy;
      if (x >= 0 && x < w && y >= 0 && y < (int)img.size() / w) {
        img[y * w + x] = (dx == 0 && dy == 0) ? center_val : neighbor_val;
      }
    }
  }
}

TEST(ImageProcessingTest, SingleStar) {
  const int W = 32, H = 32;
  auto img = make_blank(W, H);
  paint_spot(img, W, 16, 16);

  auto centroids = extract_centroids(img.data(), W, H, 100);
  ASSERT_EQ(centroids.size(), 1u);
  EXPECT_NEAR(centroids[0].x, 16.0, 0.5);
  EXPECT_NEAR(centroids[0].y, 16.0, 0.5);
}

TEST(ImageProcessingTest, TwoStars) {
  const int W = 64, H = 64;
  auto img = make_blank(W, H);
  paint_spot(img, W, 16, 16);
  paint_spot(img, W, 48, 48);

  auto centroids = extract_centroids(img.data(), W, H, 100);
  ASSERT_EQ(centroids.size(), 2u);

  // Sort by x to get deterministic order
  if (centroids[0].x > centroids[1].x)
    std::swap(centroids[0], centroids[1]);

  EXPECT_NEAR(centroids[0].x, 16.0, 0.5);
  EXPECT_NEAR(centroids[0].y, 16.0, 0.5);
  EXPECT_NEAR(centroids[1].x, 48.0, 0.5);
  EXPECT_NEAR(centroids[1].y, 48.0, 0.5);
}

TEST(ImageProcessingTest, NoStars) {
  const int W = 32, H = 32;
  auto img = make_blank(W, H);

  auto centroids = extract_centroids(img.data(), W, H, 100);
  EXPECT_EQ(centroids.size(), 0u);
}

TEST(ImageProcessingTest, NoiseRejection) {
  const int W = 64, H = 64;
  std::vector<uint8_t> img(W * H);

  // Fill with random values 0–80 (all below threshold of 100)
  std::mt19937 rng(42);
  std::uniform_int_distribution<int> dist(0, 80);
  for (auto &px : img) {
    px = static_cast<uint8_t>(dist(rng));
  }

  auto centroids = extract_centroids(img.data(), W, H, 100);
  EXPECT_EQ(centroids.size(), 0u);
}

// --- Background subtraction tests ---

TEST(ImageProcessingTest, SubtractBackgroundGradientPreservesStars) {
  const int W = 256, H = 256;
  std::vector<uint8_t> img(W * H, 0);

  // bg(x,y) = 30 + (x/4) — ramps 30..94 left to right.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      img[y * W + x] = static_cast<uint8_t>(30 + (x / 4));
    }
  }

  // Three bright 3x3 spots well inside the border margin (>=5 px from edge).
  struct Pt {
    int x, y;
  };
  Pt truth[3] = {{40, 40}, {128, 128}, {200, 200}};
  for (auto p : truth) {
    paint_spot(img, W, p.x, p.y, 230, 200);
  }

  auto bg_sub = subtract_background(img.data(), W, H, 64);

  // Far from any star, background should be near zero (within a few ADU).
  EXPECT_LT(bg_sub[10 * W + 220], 10);

  auto centroids = extract_centroids(bg_sub.data(), W, H, 50);
  ASSERT_EQ(centroids.size(), 3u);

  // Sort by x so we can compare in deterministic order.
  std::sort(centroids.begin(), centroids.end(),
            [](const StarCentroid &a, const StarCentroid &b) {
              return a.x < b.x;
            });
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(centroids[i].x, truth[i].x, 0.5);
    EXPECT_NEAR(centroids[i].y, truth[i].y, 0.5);
  }
}

TEST(ImageProcessingTest, SubtractBackgroundFlatImageGoesToZero) {
  const int W = 128, H = 128;
  std::vector<uint8_t> img(W * H, 50);
  auto bg_sub = subtract_background(img.data(), W, H, 64);
  // Should be uniformly near 0 (median of 50 is exactly 50 -> result 0).
  for (auto v : bg_sub) {
    EXPECT_LE(static_cast<int>(v), 2);
  }
}

// --- Adaptive thresholding tests ---

TEST(ImageProcessingTest, AdaptiveHandlesSplitBackground) {
  const int W = 128, H = 128;
  std::vector<uint8_t> img(W * H, 0);

  // Left half background 20, right half 120 (sigma ~ 5).
  std::mt19937 rng(7);
  std::normal_distribution<double> nL(20.0, 5.0);
  std::normal_distribution<double> nR(120.0, 5.0);
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      double v = (x < W / 2) ? nL(rng) : nR(rng);
      img[y * W + x] = static_cast<uint8_t>(std::clamp(v, 0.0, 255.0));
    }
  }
  // Paint two spots — one on each half. Center intensity = bg + ~60,
  // neighbors = bg + ~40. With sigma=5, the right-half spot (peak ~180) is
  // well above the right-half local threshold (~120 + 4*5 = 140); the
  // left-half spot (peak ~80) is well above local (~20+4*5=40) but a fixed
  // threshold=100 misses it entirely.
  paint_spot(img, W, 32, 64, 80, 60);
  paint_spot(img, W, 96, 64, 180, 160);

  auto centroids = extract_centroids_adaptive(img.data(), W, H, 4.0, 32);

  // Should find both spots.
  ASSERT_GE(centroids.size(), 2u);

  bool found_left = false, found_right = false;
  for (const auto &c : centroids) {
    if (std::abs(c.x - 32) < 1.5 && std::abs(c.y - 64) < 1.5)
      found_left = true;
    if (std::abs(c.x - 96) < 1.5 && std::abs(c.y - 64) < 1.5)
      found_right = true;
  }
  EXPECT_TRUE(found_left);
  EXPECT_TRUE(found_right);

  // Verify fixed-threshold=100 misses the left-half spot (regression on the
  // motivation for adaptive thresholding).
  auto fixed = extract_centroids(img.data(), W, H, 100);
  bool fixed_found_left = false;
  for (const auto &c : fixed) {
    if (std::abs(c.x - 32) < 2.0 && std::abs(c.y - 64) < 2.0)
      fixed_found_left = true;
  }
  EXPECT_FALSE(fixed_found_left);
}

TEST(ImageProcessingTest, AdaptiveMatchesFixedOnUniformBackground) {
  const int W = 128, H = 128;
  std::vector<uint8_t> img(W * H, 0);
  // Low uniform background, three well-separated spots inset from the edge.
  std::mt19937 rng(11);
  std::uniform_int_distribution<int> bg_dist(10, 20);
  for (auto &px : img)
    px = static_cast<uint8_t>(bg_dist(rng));
  paint_spot(img, W, 30, 30);
  paint_spot(img, W, 64, 80);
  paint_spot(img, W, 100, 30);

  auto fixed = extract_centroids(img.data(), W, H, 100);
  auto adaptive = extract_centroids_adaptive(img.data(), W, H, 5.0, 32);

  ASSERT_EQ(fixed.size(), 3u);
  ASSERT_EQ(adaptive.size(), 3u);

  auto sort_xy = [](std::vector<StarCentroid> &v) {
    std::sort(v.begin(), v.end(),
              [](const StarCentroid &a, const StarCentroid &b) {
                return a.x < b.x;
              });
  };
  sort_xy(fixed);
  sort_xy(adaptive);
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(fixed[i].x, adaptive[i].x, 1.0);
    EXPECT_NEAR(fixed[i].y, adaptive[i].y, 1.0);
  }
}

// --- Shape filter tests (apply to both APIs) ---

TEST(ImageProcessingTest, AdaptiveRejectsSingleHotPixel) {
  const int W = 32, H = 32;
  std::vector<uint8_t> img(W * H, 0);
  img[16 * W + 16] = 255;
  auto centroids = extract_centroids_adaptive(img.data(), W, H, 5.0, 16);
  EXPECT_EQ(centroids.size(), 0u);
}

TEST(ImageProcessingTest, AdaptiveRejectsStreak) {
  const int W = 32, H = 32;
  std::vector<uint8_t> img(W * H, 0);
  // 1 px tall x 10 px wide bright stripe — aspect ratio 10 > 3.
  for (int x = 11; x < 21; ++x)
    img[16 * W + x] = 220;
  auto centroids = extract_centroids_adaptive(img.data(), W, H, 5.0, 16);
  EXPECT_EQ(centroids.size(), 0u);
}

TEST(ImageProcessingTest, AdaptiveRejectsEdgeObject) {
  const int W = 32, H = 32;
  std::vector<uint8_t> img(W * H, 0);
  // 3x3 spot centered at (2, 16) — bbox x0=1 is within the 5px border margin.
  paint_spot(img, W, 2, 16, 220, 180);
  auto centroids = extract_centroids_adaptive(img.data(), W, H, 5.0, 16);
  EXPECT_EQ(centroids.size(), 0u);
}

// --- Gaussian centroiding tests (3b.1) ---

// Helper: paint a Gaussian PSF at sub-pixel (cx, cy) with given sigma and
// amplitude. Writes into a 7x7 neighborhood; pixels outside the image are
// silently dropped.
static void paint_gaussian(std::vector<uint8_t> &img, int w, int h, double cx,
                           double cy, double sigma, double amplitude) {
  int x0 = std::max(0, static_cast<int>(std::floor(cx - 3)));
  int x1 = std::min(w - 1, static_cast<int>(std::ceil(cx + 3)));
  int y0 = std::max(0, static_cast<int>(std::floor(cy - 3)));
  int y1 = std::min(h - 1, static_cast<int>(std::ceil(cy + 3)));
  double inv = 1.0 / (2.0 * sigma * sigma);
  for (int y = y0; y <= y1; ++y) {
    for (int x = x0; x <= x1; ++x) {
      double dx = x - cx;
      double dy = y - cy;
      double v = amplitude * std::exp(-(dx * dx + dy * dy) * inv);
      int pv = static_cast<int>(std::lround(v));
      img[y * w + x] = static_cast<uint8_t>(std::clamp(pv, 0, 255));
    }
  }
}

TEST(ImageProcessingTest, GaussianRefinesSubPixelCentroid) {
  const int W = 64, H = 64;
  const double truth_x = 16.3, truth_y = 16.7;
  auto img = make_blank(W, H);
  // Sigma=1.0 px PSF; amplitude 200 with threshold=70 leaves an asymmetric
  // truncated support around the peak. CoG biases noticeably toward the
  // integer grid; the intensity*Gaussian iteration pulls back to the true
  // center. This is the regime real star tracker pipelines see — the
  // threshold is well above the PSF wings, so the centroid sees only the
  // bright core.
  paint_gaussian(img, W, H, truth_x, truth_y, 1.0, 200.0);

  auto gauss = extract_centroids_gaussian(img.data(), W, H, 70);
  ASSERT_EQ(gauss.size(), 1u);
  double gauss_err = std::hypot(gauss[0].x - truth_x, gauss[0].y - truth_y);
  EXPECT_LT(gauss_err, 0.05)
      << "Gaussian centroid error too large: " << gauss_err;

  auto cog = extract_centroids(img.data(), W, H, 70);
  ASSERT_EQ(cog.size(), 1u);
  double cog_err = std::hypot(cog[0].x - truth_x, cog[0].y - truth_y);
  EXPECT_GT(cog_err, 0.05)
      << "CoG should have noticeable sub-pixel error on a non-integer-center "
         "Gaussian; got "
      << cog_err;

  // Gaussian must beat CoG.
  EXPECT_LT(gauss_err, cog_err);

  // Peak field is populated. The brightest sampled pixel of a sub-pixel-
  // centered Gaussian is below the analytic amplitude (sampling is at the
  // integer-grid nearest neighbor), so allow a comfortable floor.
  EXPECT_GT(gauss[0].peak, 150.0);
  EXPECT_LE(gauss[0].peak, 200.0);
}

TEST(ImageProcessingTest, PeakReflectsBrightestPixelInComponent) {
  const int W = 32, H = 32;
  auto img = make_blank(W, H);

  // 5x5 plateau of 255 starting at (10,10) — saturated core.
  for (int y = 10; y < 15; ++y) {
    for (int x = 10; x < 15; ++x) {
      img[y * W + x] = 255;
    }
  }
  // Adjacent halo pixel at lower intensity to enlarge the component without
  // exceeding the max-pixels filter (5x5 plateau + 1 = 26 pixels, well under
  // the default max_pixels=200).
  img[10 * W + 15] = 200;

  auto centroids = extract_centroids(img.data(), W, H, 100);
  ASSERT_EQ(centroids.size(), 1u);
  // Peak must equal 255 (the saturated core).
  EXPECT_DOUBLE_EQ(centroids[0].peak, 255.0);
  // Intensity (sum) should be well above peak — sanity check the struct.
  EXPECT_GT(centroids[0].intensity, centroids[0].peak);

  // Same expectation for the Gaussian extractor.
  auto gauss = extract_centroids_gaussian(img.data(), W, H, 100);
  ASSERT_EQ(gauss.size(), 1u);
  EXPECT_DOUBLE_EQ(gauss[0].peak, 255.0);

  // And the adaptive variants.
  auto adaptive = extract_centroids_adaptive(img.data(), W, H, 5.0, 16);
  ASSERT_GE(adaptive.size(), 1u);
  EXPECT_DOUBLE_EQ(adaptive[0].peak, 255.0);

  auto adaptive_g =
      extract_centroids_adaptive_gaussian(img.data(), W, H, 5.0, 16);
  ASSERT_GE(adaptive_g.size(), 1u);
  EXPECT_DOUBLE_EQ(adaptive_g[0].peak, 255.0);
}

#include "tiff_reader.h"
#include <cstdio>
#include <fstream>

namespace {

// Helper: write a minimal classic-TIFF (little-endian, single-strip, 16-bit
// grayscale, no compression, BlackIsZero) so the reader has something to
// parse without depending on an external image library. Returns the path of
// the temp file (caller is responsible for std::remove'ing).
std::string write_tmp_tiff16(int width, int height,
                              const std::vector<uint16_t> &pixels) {
  char tpl[] = "/tmp/startracker_tifftest_XXXXXX";
  int fd = mkstemp(tpl);
  std::string path(tpl);
  if (fd >= 0) close(fd);

  std::ofstream out(path, std::ios::binary);

  // 8-byte header: 'II' + magic=42 + IFD0 offset=8 (right after header).
  auto put16 = [&](uint16_t v) {
    out.put(static_cast<char>(v & 0xFF));
    out.put(static_cast<char>((v >> 8) & 0xFF));
  };
  auto put32 = [&](uint32_t v) {
    out.put(static_cast<char>(v & 0xFF));
    out.put(static_cast<char>((v >> 8) & 0xFF));
    out.put(static_cast<char>((v >> 16) & 0xFF));
    out.put(static_cast<char>((v >> 24) & 0xFF));
  };

  out.put('I'); out.put('I');
  put16(42);
  put32(8); // IFD0 offset

  // IFD0: 8 entries, each 12 bytes; then 4 bytes for next-IFD offset (0).
  // Layout: header(8) + ifd_count(2) + 8 entries × 12 + next_off(4) = 8 + 2 + 96 + 4 = 110.
  // Strip data starts at byte 110.
  put16(8); // entry count

  auto put_entry = [&](uint16_t tag, uint16_t typ, uint32_t cnt, uint32_t val) {
    put16(tag); put16(typ); put32(cnt); put32(val);
  };
  // Tags in ascending order per TIFF spec.
  put_entry(256, 4, 1, static_cast<uint32_t>(width));   // ImageWidth (LONG)
  put_entry(257, 4, 1, static_cast<uint32_t>(height));  // ImageLength
  put_entry(258, 3, 1, 16);                              // BitsPerSample (SHORT) — value in low bytes
  put_entry(259, 3, 1, 1);                               // Compression = none
  put_entry(262, 3, 1, 1);                               // PhotometricInterp = BlackIsZero
  put_entry(273, 4, 1, 110);                             // StripOffsets = 110 (right after IFD)
  put_entry(277, 3, 1, 1);                               // SamplesPerPixel = 1
  put_entry(279, 4, 1,
            static_cast<uint32_t>(width) *
                static_cast<uint32_t>(height) * 2); // StripByteCounts

  put32(0); // next IFD offset

  // Pixel data (raw little-endian uint16).
  for (uint16_t v : pixels) {
    out.put(static_cast<char>(v & 0xFF));
    out.put(static_cast<char>((v >> 8) & 0xFF));
  }
  return path;
}

} // namespace

// 3a.6 — round-trip a minimal 16-bit grayscale TIFF through the reader.
// We synthesize the file in-place rather than checking in a binary fixture,
// which would bloat the repo and make the test brittle to any future format
// migrations.
TEST(TiffReaderTest, RoundTripUint16) {
  const int W = 8, H = 6;
  std::vector<uint16_t> pixels(W * H);
  for (int i = 0; i < W * H; ++i) {
    // Pattern with high bytes set so endian handling is exercised.
    pixels[i] = static_cast<uint16_t>(0x1234 + i * 0x111);
  }
  std::string path = write_tmp_tiff16(W, H, pixels);

  std::vector<uint16_t> out;
  int rw = 0, rh = 0;
  ASSERT_NO_THROW(read_tiff16_grayscale(path, out, rw, rh));
  EXPECT_EQ(rw, W);
  EXPECT_EQ(rh, H);
  ASSERT_EQ(out.size(), static_cast<size_t>(W * H));
  for (int i = 0; i < W * H; ++i) {
    EXPECT_EQ(out[i], pixels[i]) << "mismatch at i=" << i;
  }

  std::remove(path.c_str());
}

// 3a.6 — reader rejects an 8-bit TIFF (BitsPerSample = 8) with a clear
// error, instead of silently producing garbage by truncating the strip.
TEST(TiffReaderTest, Rejects8BitTiff) {
  // Easiest path: synthesize a 16-bit TIFF and patch the BitsPerSample value
  // before reading. We can also just write the wrong header from scratch.
  // Do the patch route — fewer bytes touched.
  const int W = 4, H = 2;
  std::vector<uint16_t> pixels(W * H, 0);
  std::string path = write_tmp_tiff16(W, H, pixels);
  // BitsPerSample is the 3rd IFD entry, value byte at file offset
  // header(8) + ifd_count(2) + 2 entries(24) + tag/typ/cnt(8) = 42.
  std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
  f.seekp(42, std::ios::beg);
  f.put(8); f.put(0); // overwrite low two bytes of value field with 8
  f.close();

  std::vector<uint16_t> out;
  int rw = 0, rh = 0;
  EXPECT_THROW(read_tiff16_grayscale(path, out, rw, rh),
               std::runtime_error);
  std::remove(path.c_str());
}

// 3a.6 — end-to-end: feed the 16-bit overload of
// extract_centroids_adaptive_gaussian a TIFF-shaped buffer and confirm it
// finds the bright spot. This exercises the templated path through
// build_adaptive_threshold_t + extract_centroids_impl_t with T=uint16_t.
TEST(ImageProcessingTest, Adaptive16BitFindsBrightSpot) {
  const int W = 64, H = 64;
  std::vector<uint16_t> img(W * H, 1000); // background ~1000 ADU
  // Paint a 3x3 bright spot centered at (32, 32) with much higher ADU than
  // the background. 16-bit headroom: peak ~50000, neighbors ~30000.
  for (int dy = -1; dy <= 1; ++dy) {
    for (int dx = -1; dx <= 1; ++dx) {
      int x = 32 + dx, y = 32 + dy;
      img[y * W + x] =
          (dx == 0 && dy == 0) ? static_cast<uint16_t>(50000)
                                : static_cast<uint16_t>(30000);
    }
  }
  auto centroids =
      extract_centroids_adaptive_gaussian(img.data(), W, H, 5.0, 16);
  ASSERT_GE(centroids.size(), 1u);
  EXPECT_NEAR(centroids[0].x, 32.0, 0.5);
  EXPECT_NEAR(centroids[0].y, 32.0, 0.5);
  EXPECT_DOUBLE_EQ(centroids[0].peak, 50000.0);
}
