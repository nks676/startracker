#include "image_processing.h"
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
