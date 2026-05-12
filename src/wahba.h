#pragma once

#include <array>
#include <vector>

// Phase 3g.1: shared Wahba-problem primitives used by both estimation.cpp
// (public-API attitude estimate) and identification.cpp (TRIAD seed +
// post-inlier-expansion QUEST refine). Previously these were duplicated as
// `estimate_attitude_triad` / `quest_attitude_local`; the duplication caused
// every numerical tweak to be made twice and tested twice.
//
// Convention: all rotations are inertial->camera (so v_camera = R · v_inertial
// and quaternion (qx, qy, qz, qw) rotates inertial unit vectors into camera
// unit vectors).

namespace wahba {

// TRIAD on two camera-frame / inertial-frame vector pairs. Output is the 3x3
// inertial->camera rotation. Identical algebra to the previous private
// triad_rotation; exposed publicly here.
void triad(const std::array<double, 3> &W1,
           const std::array<double, 3> &W2,
           const std::array<double, 3> &V1,
           const std::array<double, 3> &V2,
           double R_out[3][3]);

// QUEST (Shuster & Oh 1981) on M >= 2 camera/inertial unit-vector pairs with
// equal weights 1/M. Writes the optimal quaternion (x, y, z, w) into q_out
// with q_out[3] (w) >= 0 on success. Returns false if Newton-Raphson fails
// to converge or produces non-finite values; q_out is then left unmodified.
//
// `max_iters` and `tol_lambda` are kept caller-tunable because the two
// historical sites used slightly different settings (20 / 1e-12 for the
// public attitude path, 30 / 1e-14 for the tight identification-side
// refinement). Defaults match the public path.
bool quest(const std::vector<std::array<double, 3>> &v_cam,
           const std::vector<std::array<double, 3>> &v_iner,
           double q_out[4],
           int max_iters = 20,
           double tol_lambda = 1e-12);

// Convenience wrapper: run quest() then convert the resulting quaternion to a
// 3x3 rotation matrix (inertial -> camera). Returns whatever quest() returned;
// R_out is unmodified on failure.
bool quest_R(const std::vector<std::array<double, 3>> &v_cam,
             const std::vector<std::array<double, 3>> &v_iner,
             double R_out[3][3],
             int max_iters = 20,
             double tol_lambda = 1e-12);

// Convert an inertial->camera quaternion (x, y, z, w) to a 3x3 rotation
// matrix. Caller is responsible for passing a normalised input; the
// conversion does not re-normalise.
void quat_to_rotation(const double q[4], double R_out[3][3]);

} // namespace wahba
