#include "image_processing.h"
#include <queue>
#include <vector>

std::vector<StarCentroid> extract_centroids(const uint8_t *image_data,
                                            int width, int height,
                                            uint8_t threshold) {
  std::vector<StarCentroid> centroids;
  std::vector<bool> visited(width * height, false);

  // Neighborhood offsets (8-connected)
  const int dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};
  const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int idx = y * width + x;

      if (!visited[idx] && image_data[idx] > threshold) {
        // New connected component found
        std::queue<std::pair<int, int>> q;
        q.push({x, y});
        visited[idx] = true;

        double sum_x = 0;
        double sum_y = 0;
        double sum_i = 0;
        int pixel_count = 0;

        while (!q.empty()) {
          auto [cx, cy] = q.front();
          q.pop();

          int cidx = cy * width + cx;
          double intensity = static_cast<double>(image_data[cidx]);

          sum_x += cx * intensity;
          sum_y += cy * intensity;
          sum_i += intensity;
          pixel_count++;

          // Check neighbors
          for (int n = 0; n < 8; ++n) {
            int nx = cx + dx[n];
            int ny = cy + dy[n];

            if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
              int nidx = ny * width + nx;
              if (!visited[nidx] && image_data[nidx] > threshold) {
                visited[nidx] = true;
                q.push({nx, ny});
              }
            }
          }
        }

        // Only keep components that have enough pixels
        if (pixel_count > 1 && sum_i > 0) {
          StarCentroid centroid;
          centroid.x = sum_x / sum_i;
          centroid.y = sum_y / sum_i;
          centroid.intensity = sum_i;
          centroids.push_back(centroid);
        }
      }
    }
  }

  return centroids;
}
