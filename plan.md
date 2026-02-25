# Project Plan: *startracker*

This plan outlines the iterative development approach for the highly robust embedded star tracking engine. We will start with a simple, accurate baseline, establish rigorous testing, and progressively upgrade the algorithms for performance.

## Phase 1: Bare Minimum Baseline (Slow but Accurate)
The goal is to get an end-to-end pipeline working with the simplest possible algorithms to ensure accuracy before optimizing for speed.

*   **Image Processing:** Standard Center of Gravity (CoG) for centroiding instead of Fast Gaussian Fitting.
*   **Identification:** Simple brute-force pattern matching or basic subgraph matching, rather than optimized $k$-vector database lookups.
*   **Attitude Estimation:** Basic TRIAD or Singular Value Decomposition (SVD) algorithm to solve Wahba's problem, instead of QUEST.

## Phase 2: End-to-End CI/CD Infrastructure
Before making the system complex, we need an automated safety net to ensure continued accuracy as we optimize.

*   Set up automated unit testing for the C++ core.
*   Build a CI/CD pipeline (e.g., GitHub Actions) to compile the code and run tests on every commit.
*   Implement baseline Monte Carlo simulations in Python to validate the accuracy of the slow-but-accurate algorithms against synthetic star fields.

## Phase 3: Algorithm Optimization (Faster & More Robust)
With CI/CD in place, we will iteratively replace the basic algorithms with state-of-the-art, high-performance versions, validating against our tests along the way.

*   **Centroiding Upgrade:** Replace standard CoG with Fast Gaussian Fitting (FGF) for fast sub-pixel precision.
*   **Search & Match Upgrade:** Implement $k$-vector $O(1)$ search lookups and Geometric Voting for robust "Lost-in-Space" identification.
*   **Estimation Upgrade:** Replace the TRIAD/SVD solver entirely with the highly optimized QUEST (QUaternion ESTimator) algorithm.
