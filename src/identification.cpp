#include "identification.h"
#include <cmath>
#include <map>

std::vector<IdentifiedStar>
identify_stars(const std::vector<StarCentroid> &image_stars,
               const PinholeCamera &camera, const StarDatabase &db,
               double cos_tolerance) {
  int N = image_stars.size();
  if (N < 3)
    return {}; // Need at least 3 stars for reliable identification

  std::vector<std::vector<double>> v_cam(N, std::vector<double>(3));

  // Convert centroids to camera frame unit vectors
  for (int i = 0; i < N; ++i) {
    double vx = (image_stars[i].x - camera.center_x) / camera.focal_x;
    double vy = (image_stars[i].y - camera.center_y) / camera.focal_y;
    double vz = 1.0;

    double norm = std::sqrt(vx * vx + vy * vy + vz * vz);
    v_cam[i][0] = vx / norm;
    v_cam[i][1] = vy / norm;
    v_cam[i][2] = vz / norm;
  }

  // Voting matrix: map[image_idx][catalog_hip] = votes
  std::vector<std::map<int, int>> votes(N);

  // Compute pairwise dot products and vote
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      double dot = v_cam[i][0] * v_cam[j][0] + v_cam[i][1] * v_cam[j][1] +
                   v_cam[i][2] * v_cam[j][2];

      auto matches = db.find_pairs(dot, cos_tolerance);

      for (const auto &match : matches) {
        // Vote both possible assignments: (id1→i, id2→j) and (id1→j, id2→i)
        votes[i][match.id1]++;
        votes[j][match.id2]++;
        votes[i][match.id2]++;
        votes[j][match.id1]++;
      }
    }
  }

  // Determine the highest voted catalog star for each image star
  // Store best_hip and max_votes per image star for duplicate resolution
  std::vector<int> best_hips(N, -1);
  std::vector<int> max_votes_per_star(N, 0);

  for (int i = 0; i < N; ++i) {
    for (auto const &[hip, count] : votes[i]) {
      if (count > max_votes_per_star[i]) {
        max_votes_per_star[i] = count;
        best_hips[i] = hip;
      }
    }
  }

  // Resolve duplicate catalog assignments: if two image stars claim the same
  // catalog HIP, keep only the one with the higher vote count.
  std::map<int, int> hip_to_best_image; // catalog_hip -> image_idx
  for (int i = 0; i < N; ++i) {
    if (max_votes_per_star[i] < 3 || best_hips[i] == -1)
      continue;

    int hip = best_hips[i];
    auto it = hip_to_best_image.find(hip);
    if (it == hip_to_best_image.end()) {
      hip_to_best_image[hip] = i;
    } else {
      // Keep the image star with the higher vote count
      int prev = it->second;
      if (max_votes_per_star[i] > max_votes_per_star[prev]) {
        best_hips[prev] = -1; // discard previous
        hip_to_best_image[hip] = i;
      } else {
        best_hips[i] = -1; // discard current
      }
    }
  }

  // Build result from resolved assignments
  std::vector<IdentifiedStar> identified;
  for (int i = 0; i < N; ++i) {
    if (max_votes_per_star[i] >= 3 && best_hips[i] != -1) {
      IdentifiedStar is;
      is.image_idx = i;
      is.catalog_hip_id = best_hips[i];
      is.v_cam[0] = v_cam[i][0];
      is.v_cam[1] = v_cam[i][1];
      is.v_cam[2] = v_cam[i][2];
      identified.push_back(is);
    }
  }

  return identified;
}
