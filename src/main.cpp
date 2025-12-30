#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <cassert>

#include "fits/fits_io.h"
#include "catalog/catalog.h"
#include "triangle/triangle.h"

// ---------------------------------------------------------
// Test helpers for robustness checks
// ---------------------------------------------------------
static void normalize(Star& s) {
    double mag = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
    if (mag == 0.0) return;
    s.x /= mag;
    s.y /= mag;
    s.z /= mag;
}

static Star random_star(int id, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Star s{id, dist(rng), dist(rng), dist(rng), 0.0};
    normalize(s);
    return s;
}

static Star perturb_star(const Star& original, double noise_rad, std::mt19937& rng) {
    std::normal_distribution<double> noise_dist(0.0, noise_rad);
    Star noise{0, noise_dist(rng), noise_dist(rng), noise_dist(rng), 0.0};
    Star s{original.id, original.x + noise.x, original.y + noise.y, original.z + noise.z, original.magnitude};
    normalize(s);
    return s;
}

static void test_gaussian_noise(const std::vector<Triangle>& db, const std::vector<Star>& catalog, std::mt19937& rng) {
    std::cout << "\n[TEST] Gaussian Noise Robustness..." << std::endl;
    if (catalog.size() < 3) {
        std::cout << "  SKIP: Not enough stars in test catalog." << std::endl;
        return;
    }

    Star s1 = catalog[0];
    Star s2 = catalog[1];
    Star s3 = catalog[2];

    const double noise_level = 0.0007; // ~2 pixels in a typical sensor model
    Triangle match = find_triangle(
        perturb_star(s1, noise_level, rng),
        perturb_star(s2, noise_level, rng),
        perturb_star(s3, noise_level, rng),
        db);

    if (match.star1 != -1) {
        std::cout << "  PASS: Match found despite noise (" << noise_level << " rad)." << std::endl;
        std::cout << "  Matched IDs: " << match.star1 << ", " << match.star2 << ", " << match.star3 << std::endl;
    } else {
        std::cout << "  FAIL: Algorithm lost the star due to noise!" << std::endl;
    }
}

static void test_false_star(const std::vector<Triangle>& db, const std::vector<Star>& catalog, std::mt19937& rng) {
    std::cout << "\n[TEST] False Star (Space Junk) Rejection..." << std::endl;
    if (catalog.size() < 2) {
        std::cout << "  SKIP: Not enough stars in test catalog." << std::endl;
        return;
    }

    Star s1 = catalog[0];
    Star s2 = catalog[1];
    Star junk = random_star(9999, rng);

    Triangle match = find_triangle(s1, s2, junk, db);
    if (match.star1 == -1) {
        std::cout << "  PASS: Correctly rejected invalid triangle." << std::endl;
    } else {
        std::cout << "  FAIL: False positive! Matched junk to IDs: "
                  << match.star1 << ", " << match.star2 << ", " << match.star3 << std::endl;
    }
}

static void test_magnitude_dropout(const std::vector<Triangle>& db, const std::vector<Star>& catalog) {
    std::cout << "\n[TEST] Magnitude Dropout (Alternative Triangle)..." << std::endl;
    if (catalog.size() < 4) {
        std::cout << "  SKIP: Not enough stars in test catalog." << std::endl;
        return;
    }

    Star s1 = catalog[0];
    Star s2 = catalog[1];
    Star s4 = catalog[3];

    Triangle match = find_triangle(s1, s2, s4, db);
    if (match.star1 != -1) {
        std::cout << "  PASS: Successfully identified backup triangle (A, B, D)." << std::endl;
    } else {
        std::cout << "  FAIL: Could not identify backup triangle." << std::endl;
    }
}

static std::vector<Star> build_test_catalog() {
    std::vector<Star> catalog;
    catalog.push_back({101, 0.99, 0.10, 0.05, 0.0});
    catalog.push_back({102, 0.98, 0.15, 0.06, 0.0});
    catalog.push_back({103, 0.985, 0.12, 0.09, 0.0});
    catalog.push_back({104, 0.97, 0.11, 0.04, 0.0});
    catalog.push_back({105, 0.975, 0.18, 0.07, 0.0});

    for (auto& s : catalog) {
        normalize(s);
    }
    return catalog;
}

static std::vector<Triangle> build_test_db(const std::vector<Star>& catalog) {
    std::vector<Triangle> db;
    int n = static_cast<int>(catalog.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                double d1 = std::acos(std::clamp(catalog[i].x * catalog[j].x + catalog[i].y * catalog[j].y + catalog[i].z * catalog[j].z, -1.0, 1.0));
                double d2 = std::acos(std::clamp(catalog[i].x * catalog[k].x + catalog[i].y * catalog[k].y + catalog[i].z * catalog[k].z, -1.0, 1.0));
                double d3 = std::acos(std::clamp(catalog[j].x * catalog[k].x + catalog[j].y * catalog[k].y + catalog[j].z * catalog[k].z, -1.0, 1.0));

                std::vector<double> s = {d1, d2, d3};
                std::sort(s.begin(), s.end());
                db.push_back({catalog[i].id, catalog[j].id, catalog[k].id, s[0], s[1], s[2]});
            }
        }
    }

    std::sort(db.begin(), db.end(), [](const Triangle& t1, const Triangle& t2) {
        return t1.a < t2.a;
    });

    return db;
}

int main() {
    std::string m42_filename = "data/m42_40min_red.fits";
    ImageData m42 = fits_to_data(m42_filename);
    std::cout << "Loaded " << m42.clusters.size() << " clusters from " << m42_filename << std::endl;

    std::string m13_filename = "data/m13_test.fits";
    ImageData m13 = fits_to_data(m13_filename);
    std::cout << "Loaded " << m13.clusters.size() << " clusters from " << m13_filename << std::endl;

    std::string starmap_filename = "data/hipparcos-voidmain.csv";
    std::vector<Star> catalog = csv_to_catalog(starmap_filename);
    std::cout << "Loaded " << catalog.size() << " stars from " << starmap_filename << std::endl;

    std::vector<Triangle> triangles = catalog_to_triangles(catalog);
    std::cout << "Generated " << triangles.size() << " triangles from catalog" << std::endl;

    // ---------------------------------------------------------
    // PHASE 1.2: SYNTHETIC STRESS TEST
    // ---------------------------------------------------------
    std::cout << "\n--- STARTING SYNTHETIC STRESS TEST ---\n";

    if (triangles.empty()) {
        std::cerr << "Error: No triangles generated. Cannot run test." << std::endl;
        return -1;
    }

    int passed = 0;
    int failed = 0;
    int total_tests = 0;

    // We test every 1000th triangle to save time
    // (Testing all 1.1 million would take a few minutes)
    for (size_t i = 0; i < triangles.size(); i += 1000) {
        total_tests++;
        Triangle expected = triangles[i];

        // 1. Reconstruct the 3 stars from the IDs
        // (In a real scenario, we don't know the IDs yet, just the positions)
        Star s1, s2, s3;
        
        // Simple linear search to find the stars by ID
        // (Optimization: A hash map would be faster, but this is fine for testing)
        int found_count = 0;
        for(const auto& s : catalog) {
            if (s.id == expected.star1) { s1 = s; found_count++; }
            else if (s.id == expected.star2) { s2 = s; found_count++; }
            else if (s.id == expected.star3) { s3 = s; found_count++; }
            if (found_count == 3) break;
        }

        // 2. Run the Matcher
        // We pass the 3 star objects. The matcher measures them and searches the DB.
        Triangle result = find_triangle(s1, s2, s3, triangles);

        // 3. Verify Result
        bool match_found = (result.star1 != -1);
        
        // Check if the IDs match the expected ones (ignoring order)
        bool ids_match = false;
        if (match_found) {
            std::vector<int> exp_ids = {expected.star1, expected.star2, expected.star3};
            std::vector<int> res_ids = {result.star1, result.star2, result.star3};
            std::sort(exp_ids.begin(), exp_ids.end());
            std::sort(res_ids.begin(), res_ids.end());
            
            if (exp_ids == res_ids) ids_match = true;
        }

        if (match_found && ids_match) {
            passed++;
        } else {
            failed++;
        }
        
        // Progress Indicator (dots)
        if (total_tests % 100 == 0) std::cout << "." << std::flush;
    }

    std::cout << "\n\n--- RESULTS ---\n";
    std::cout << "Total Tests: " << total_tests << std::endl;
    std::cout << "Passed:      " << passed << std::endl;
    std::cout << "Failed:      " << failed << std::endl;
    std::cout << "Accuracy:    " << (double)passed / total_tests * 100.0 << "%" << std::endl;

    // ---------------------------------------------------------
    // Additional robustness tests on a synthetic mini-catalog
    // ---------------------------------------------------------
    std::mt19937 rng(42);
    std::vector<Star> test_catalog = build_test_catalog();
    std::vector<Triangle> test_db = build_test_db(test_catalog);
    std::cout << "\nBuilt synthetic test DB with " << test_db.size() << " triangles." << std::endl;

    test_gaussian_noise(test_db, test_catalog, rng);
    test_false_star(test_db, test_catalog, rng);
    test_magnitude_dropout(test_db, test_catalog);

    return 0;
} 