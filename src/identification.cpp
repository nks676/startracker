#include "identification.h"
#include <cmath>
#include <iostream>
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
        votes[i][match.id1]++;
        votes[i][match.id2]++;
        votes[j][match.id1]++;
        votes[j][match.id2]++;
      }
    }
  }

  std::vector<IdentifiedStar> identified;

  // Determine the highest voted catalog star for each image star
  for (int i = 0; i < N; ++i) {
    int best_hip = -1;
    int max_votes = 0;

    for (auto const &[hip, count] : votes[i]) {
      if (count > max_votes) {
        max_votes = count;
        best_hip = hip;
      }
    }

    // A true star match will have roughly N-1 votes if all lines match.
    // We require at least 3 votes for baseline robustness.
    if (max_votes >= 3 && best_hip != -1) {
      IdentifiedStar is;
      is.image_idx = i;
      is.catalog_hip_id = best_hip;
      is.v_cam[0] = v_cam[i][0];
      is.v_cam[1] = v_cam[i][1];
      is.v_cam[2] = v_cam[i][2];
      identified.push_back(is);
    }
  }

  // TODO: Resolve duplicate catalog assignments (if two image stars claim the
  // same catalog star)

  return identified;
}
