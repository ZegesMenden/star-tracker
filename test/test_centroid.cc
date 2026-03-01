// Unit tests for CenterOfGravityCentroid (centroid.hpp / centroid.cc)

#include <catch.hpp>
#include "centroid.hpp"
#include <cstring>
#include <cmath>
#include <vector>

// Helper: create a blank image
static std::vector<uint8_t> blankImage(size_t w, size_t h) {
    return std::vector<uint8_t>(w * h, 0);
}

// Helper: draw a filled circle (approximate) into an image
static void drawStar(std::vector<uint8_t> &img, size_t w, size_t h,
                     int cx, int cy, int radius, uint8_t brightness) {
    for (int dy = -radius; dy <= radius; dy++) {
        for (int dx = -radius; dx <= radius; dx++) {
            if (dx * dx + dy * dy <= radius * radius) {
                int x = cx + dx;
                int y = cy + dy;
                if (x >= 0 && x < (int)w && y >= 0 && y < (int)h) {
                    img[y * w + x] = brightness;
                }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Blank image → no centroids
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: blank image yields no centroids", "[centroid]") {
    auto img = blankImage(64, 64);
    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), 64, 64);
    CHECK(result.empty());
}

// ═══════════════════════════════════════════════════════════════════════════
//  Uniform image → no centroids (everything below threshold)
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: uniform low-intensity image yields no centroids", "[centroid]") {
    auto img = blankImage(64, 64);
    std::memset(img.data(), 10, img.size());  // uniform low brightness
    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), 64, 64);
    CHECK(result.empty());
}

// ═══════════════════════════════════════════════════════════════════════════
//  Single bright star
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: single bright star detected", "[centroid]") {
    size_t w = 64, h = 64;
    auto img = blankImage(w, h);
    drawStar(img, w, h, 32, 32, 3, 200);

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    REQUIRE(result.size() == 1);
    CHECK(result[0].x == Approx(32.5).margin(1.0));
    CHECK(result[0].y == Approx(32.5).margin(1.0));
    CHECK(result[0].numPixels > 0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Multiple stars
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: two well-separated stars", "[centroid]") {
    size_t w = 128, h = 128;
    auto img = blankImage(w, h);
    drawStar(img, w, h, 30, 30, 3, 255);
    drawStar(img, w, h, 90, 90, 3, 255);

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    REQUIRE(result.size() == 2);

    // Sort by x position for deterministic checks
    if (result[0].x > result[1].x) std::swap(result[0], result[1]);

    CHECK(result[0].x == Approx(30.5).margin(1.5));
    CHECK(result[0].y == Approx(30.5).margin(1.5));
    CHECK(result[1].x == Approx(90.5).margin(1.5));
    CHECK(result[1].y == Approx(90.5).margin(1.5));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Star on border is rejected
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: star touching the image border is excluded", "[centroid]") {
    size_t w = 64, h = 64;
    auto img = blankImage(w, h);
    // Place a star right on the left edge
    drawStar(img, w, h, 0, 32, 3, 255);

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    // Should be rejected because it touches x=0
    CHECK(result.empty());
}

// ═══════════════════════════════════════════════════════════════════════════
//  Centroid accuracy with asymmetric brightness
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: weighted position with asymmetric brightness", "[centroid]") {
    size_t w = 64, h = 64;
    auto img = blankImage(w, h);

    // Manually place a 3×1 horizontal bar with non-uniform brightness:
    // pixel at (29, 30) = 100, (30, 30) = 200, (31, 30) = 100
    img[30 * w + 29] = 100;
    img[30 * w + 30] = 200;
    img[30 * w + 31] = 100;

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    // The weighted CoG should be centered at x=30 (heavier pixel there)
    REQUIRE(result.size() == 1);
    CHECK(result[0].x == Approx(30.5).margin(0.5));
    CHECK(result[0].y == Approx(30.5).margin(0.5));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Radius calculation
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: radius reflects bounding box", "[centroid]") {
    size_t w = 64, h = 64;
    auto img = blankImage(w, h);

    // 5-pixel horizontal line
    for (int x = 28; x <= 32; x++) {
        img[30 * w + x] = 250;
    }

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    REQUIRE(result.size() == 1);
    // radiusX = (32-28+1)/2 = 2.5, radiusY = 1/2 = 0.5
    CHECK(result[0].radiusX == Approx(2.5));
    CHECK(result[0].radiusY == Approx(0.5));
    CHECK(result[0].numPixels == 5);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Single maximum-brightness pixel
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Centroid: single bright pixel", "[centroid]") {
    size_t w = 32, h = 32;
    auto img = blankImage(w, h);
    img[16 * w + 16] = 255;

    CenterOfGravityCentroid algo;
    auto result = algo.compute(img.data(), w, h);

    REQUIRE(result.size() == 1);
    CHECK(result[0].x == Approx(16.5).margin(0.01));
    CHECK(result[0].y == Approx(16.5).margin(0.01));
    CHECK(result[0].numPixels == 1);
}
