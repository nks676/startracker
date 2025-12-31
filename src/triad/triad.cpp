#include "triad.h"

#include <cmath>
#include <algorithm>

// --- Helper Functions ---

// Normalize the Star vector (x, y, z)
void normalize(Star& v) {
    double mag = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (mag > 1e-9) {
        v.x /= mag; v.y /= mag; v.z /= mag;
    }
}

// Cross Product (Vector output does not need ID or Magnitude)
Star cross(Star a, Star b) {
    return {
        0, // ID unused
        (double) a.y * b.z - a.z * b.y,
        (double) a.z * b.x - a.x * b.z,
        (double) a.x * b.y - a.y * b.x,
        0.0 // Magnitude unused
    };
}

// --- Main TRIAD Algorithm ---

Quaternion compute_attitude(const std::vector<Observation>& obs) {
    // TRIAD requires at least 2 stars
    if (obs.size() < 2) return {1, 0, 0, 0};

    // 1. Extract vectors from the first two observations
    // We ignore 'magnitude' and 'id' for the math part
    Star r1 = obs[0].inertial;
    Star r2 = obs[1].inertial;
    Star b1 = obs[0].body;
    Star b2 = obs[1].body;

    // 2. Construct the Inertial Reference Triad (V1, V2, V3)
    // V1 = r1 (normalized)
    Star V1 = r1; 
    normalize(V1);

    // V2 = (r1 x r2) / |r1 x r2|  (Unit Normal to the plane)
    Star V2 = cross(r1, r2);
    normalize(V2);

    // V3 = V1 x V2 (Completes the orthonormal basis)
    Star V3 = cross(V1, V2);

    // 3. Construct the Body Measurement Triad (W1, W2, W3)
    // We use the EXACT same logic so the frames match
    Star W1 = b1; 
    normalize(W1);

    Star W2 = cross(b1, b2);
    normalize(W2);

    Star W3 = cross(W1, W2);

    // 4. Calculate Rotation Matrix A = [W] * [V]^T
    // The matrix A transforms Inertial -> Body
    // A = W * V_transpose
    
    // Construct 3x3 matrices from the column vectors
    double W_mat[3][3] = {
        {W1.x, W2.x, W3.x}, 
        {W1.y, W2.y, W3.y}, 
        {W1.z, W2.z, W3.z}
    };
    
    double V_mat[3][3] = {
        {V1.x, V2.x, V3.x}, 
        {V1.y, V2.y, V3.y}, 
        {V1.z, V2.z, V3.z}
    };

    // Perform Multiplication: A = W_mat * Transpose(V_mat)
    double A[3][3] = {{0}};
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            for(int k=0; k<3; k++) {
                // Notice V_mat[j][k] is used instead of V_mat[k][j] because we want the transpose
                A[i][j] += W_mat[i][k] * V_mat[j][k];
            }
        }
    }

    // 5. Convert Rotation Matrix A to Quaternion
    // We use "Stanley's Method" (checking the trace) to avoid dividing by zero 
    // This handles all rotation angles safely.
    double tr = A[0][0] + A[1][1] + A[2][2];
    Quaternion q;

    if (tr > 0) {
        double S = std::sqrt(tr + 1.0) * 2; // S=4*qw
        q.w = 0.25 * S;
        q.x = (A[2][1] - A[1][2]) / S;
        q.y = (A[0][2] - A[2][0]) / S;
        q.z = (A[1][0] - A[0][1]) / S;
    } else if ((A[0][0] > A[1][1]) && (A[0][0] > A[2][2])) {
        double S = std::sqrt(1.0 + A[0][0] - A[1][1] - A[2][2]) * 2; // S=4*qx
        q.w = (A[2][1] - A[1][2]) / S;
        q.x = 0.25 * S;
        q.y = (A[0][1] + A[1][0]) / S;
        q.z = (A[0][2] + A[2][0]) / S;
    } else if (A[1][1] > A[2][2]) {
        double S = std::sqrt(1.0 + A[1][1] - A[0][0] - A[2][2]) * 2; // S=4*qy
        q.w = (A[0][2] - A[2][0]) / S;
        q.x = (A[0][1] + A[1][0]) / S;
        q.y = 0.25 * S;
        q.z = (A[1][2] + A[2][1]) / S;
    } else {
        double S = std::sqrt(1.0 + A[2][2] - A[0][0] - A[1][1]) * 2; // S=4*qz
        q.w = (A[1][0] - A[0][1]) / S;
        q.x = (A[0][2] + A[2][0]) / S;
        q.y = (A[1][2] + A[2][1]) / S;
        q.z = 0.25 * S;
    }

    return q;
}