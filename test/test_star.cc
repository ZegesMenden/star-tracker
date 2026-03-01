// Unit tests for CatalogStar, Star, and StarIdentifier types

#include <catch.hpp>
#include "star.hpp"
#include "test_precision.hpp"
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════
//  CatalogStar
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("CatalogStar from ra/de produces unit spatial vector", "[catalog]") {
    CatalogStar s(1.0, 0.5, 300, 42);
    CHECK(s.spatial.magnitude() == Approx(1.0).margin(TIGHT_EPS));
    CHECK(s.magnitude == 300);
    CHECK(s.name == 42);
}

TEST_CASE("CatalogStar from ra=0 de=0 points along x-axis", "[catalog]") {
    CatalogStar s(0.0, 0.0, 100, 1);
    CHECK(s.spatial.x == Approx(1.0));
    CHECK(s.spatial.y == Approx(0.0).margin(NEAR_ZERO));
    CHECK(s.spatial.z == Approx(0.0).margin(NEAR_ZERO));
}

TEST_CASE("CatalogStar from Vec3 stores spatial directly", "[catalog]") {
    Vec3 v{0.577, 0.577, 0.577};
    CatalogStar s(v, 200, 7);
    CHECK(s.spatial.x == Approx(v.x));
    CHECK(s.spatial.y == Approx(v.y));
    CHECK(s.spatial.z == Approx(v.z));
}

TEST_CASE("Angular distance between CatalogStars is consistent", "[catalog]") {
    // Two catalog stars at known angular separation
    CatalogStar s1(0.0, 0.0, 100, 1);    // along x
    CatalogStar s2(M_PI / 4, 0.0, 100, 2); // 45° around equator
    stfloat d = angleUnit(s1.spatial, s2.spatial);
    CHECK(d == Approx(M_PI / 4.0).margin(LOOSE_EPS));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Star (centroid)
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Star default constructor", "[star]") {
    Star s;
    CHECK(s.position.x == Approx(0.0));
    CHECK(s.position.y == Approx(0.0));
    CHECK(s.radiusX == Approx(0.0));
    CHECK(s.radiusY == Approx(0.0));
    CHECK(s.magnitude == 0);
}

TEST_CASE("Star 5-arg constructor", "[star]") {
    Star s(10.5, 20.5, 3.0, 4.0, 150);
    CHECK(s.position.x == Approx(10.5));
    CHECK(s.position.y == Approx(20.5));
    CHECK(s.radiusX == Approx(3.0));
    CHECK(s.radiusY == Approx(4.0));
    CHECK(s.magnitude == 150);
}

TEST_CASE("Star 3-arg constructor sets radiusX==radiusY", "[star]") {
    Star s(5.0, 6.0, 2.0);
    CHECK(s.radiusX == Approx(s.radiusY));
    CHECK(s.magnitude == 0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  StarIdentifier
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("StarIdentifier stores indices", "[starid]") {
    StarIdentifier id(3, 42, 0.95);
    CHECK(id.starIndex == 3);
    CHECK(id.catalogIndex == 42);
    CHECK(id.weight == Approx(0.95));
}

TEST_CASE("StarIdentifier default weight", "[starid]") {
    StarIdentifier id(0, 1);
    CHECK(id.weight == Approx(1.0));
}

TEST_CASE("StarIdentifier equality", "[starid]") {
    StarIdentifier a(1, 2);
    StarIdentifier b(1, 2, 0.5);  // weight doesn't affect equality
    StarIdentifier c(1, 3);
    CHECK(a == b);
    CHECK(!(a == c));
}
