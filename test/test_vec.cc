// Unit tests for Vec2, Vec3, and angle utilities (star.hpp / star.cc)

#include <catch.hpp>
#include "star.hpp"
#include "test_precision.hpp"
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════
//  Vec2 tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Vec2 magnitude", "[vec2]") {
    Vec2 v{3.0, 4.0};
    CHECK(v.magnitude() == Approx(5.0));
    CHECK(v.magnitudeSq() == Approx(25.0));
}

TEST_CASE("Vec2 zero magnitude", "[vec2]") {
    Vec2 v{0.0, 0.0};
    CHECK(v.magnitude() == Approx(0.0));
    CHECK(v.magnitudeSq() == Approx(0.0));
}

TEST_CASE("Vec2 normalize", "[vec2]") {
    Vec2 v{3.0, 4.0};
    Vec2 n = v.normalize();
    CHECK(n.x == Approx(0.6));
    CHECK(n.y == Approx(0.8));
    CHECK(n.magnitude() == Approx(1.0));
}

TEST_CASE("Vec2 dot product", "[vec2]") {
    Vec2 a{1.0, 2.0};
    Vec2 b{3.0, 4.0};
    CHECK((a * b) == Approx(11.0));  // 1*3 + 2*4
}

TEST_CASE("Vec2 perpendicular dot product is zero", "[vec2]") {
    Vec2 a{1.0, 0.0};
    Vec2 b{0.0, 1.0};
    CHECK((a * b) == Approx(0.0));
}

TEST_CASE("Vec2 scalar multiply", "[vec2]") {
    Vec2 v{2.0, 3.0};
    Vec2 r = v * 5.0;
    CHECK(r.x == Approx(10.0));
    CHECK(r.y == Approx(15.0));
}

TEST_CASE("Vec2 addition", "[vec2]") {
    Vec2 a{1.0, 2.0};
    Vec2 b{3.0, 4.0};
    Vec2 c = a + b;
    CHECK(c.x == Approx(4.0));
    CHECK(c.y == Approx(6.0));
}

TEST_CASE("Vec2 subtraction", "[vec2]") {
    Vec2 a{5.0, 7.0};
    Vec2 b{3.0, 2.0};
    Vec2 c = a - b;
    CHECK(c.x == Approx(2.0));
    CHECK(c.y == Approx(5.0));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Vec3 tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Vec3 magnitude", "[vec3]") {
    Vec3 v{1.0, 2.0, 2.0};
    CHECK(v.magnitude() == Approx(3.0));
    CHECK(v.magnitudeSq() == Approx(9.0));
}

TEST_CASE("Vec3 normalize", "[vec3]") {
    Vec3 v{0.0, 0.0, 5.0};
    Vec3 n = v.normalize();
    CHECK(n.x == Approx(0.0));
    CHECK(n.y == Approx(0.0));
    CHECK(n.z == Approx(1.0));
    CHECK(n.magnitude() == Approx(1.0));
}

TEST_CASE("Vec3 dot product", "[vec3]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    CHECK((a * b) == Approx(32.0));  // 1*4 + 2*5 + 3*6
}

TEST_CASE("Vec3 cross product", "[vec3]") {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    Vec3 c = a.cross(b);
    CHECK(c.x == Approx(0.0));
    CHECK(c.y == Approx(0.0));
    CHECK(c.z == Approx(1.0));  // i × j = k
}

TEST_CASE("Vec3 cross product anticommutative", "[vec3]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 ab = a.cross(b);
    Vec3 ba = b.cross(a);
    CHECK(ab.x == Approx(-ba.x));
    CHECK(ab.y == Approx(-ba.y));
    CHECK(ab.z == Approx(-ba.z));
}

TEST_CASE("Vec3 cross product perpendicular to both inputs", "[vec3]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 c = a.cross(b);
    CHECK((a * c) == Approx(0.0).margin(TIGHT_EPS));
    CHECK((b * c) == Approx(0.0).margin(TIGHT_EPS));
}

TEST_CASE("Vec3 scalar multiply", "[vec3]") {
    Vec3 v{1.0, 2.0, 3.0};
    Vec3 r = v * 2.0;
    CHECK(r.x == Approx(2.0));
    CHECK(r.y == Approx(4.0));
    CHECK(r.z == Approx(6.0));
}

TEST_CASE("Vec3 addition and subtraction", "[vec3]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    Vec3 sum = a + b;
    Vec3 diff = a - b;
    CHECK(sum.x == Approx(5.0));
    CHECK(sum.y == Approx(7.0));
    CHECK(sum.z == Approx(9.0));
    CHECK(diff.x == Approx(-3.0));
    CHECK(diff.y == Approx(-3.0));
    CHECK(diff.z == Approx(-3.0));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Angle utilities
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("angle() between identical vectors is zero", "[angle]") {
    Vec3 a{1.0, 0.0, 0.0};
    CHECK(angle(a, a) == Approx(0.0).margin(LOOSE_EPS));
}

TEST_CASE("angle() between orthogonal vectors is pi/2", "[angle]") {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    CHECK(angle(a, b) == Approx(M_PI / 2.0).margin(LOOSE_EPS));
}

TEST_CASE("angle() between opposite vectors is pi", "[angle]") {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{-1.0, 0.0, 0.0};
    CHECK(angle(a, b) == Approx(M_PI).margin(LOOSE_EPS));
}

TEST_CASE("angle() with unnormalized vectors", "[angle]") {
    Vec3 a{10.0, 0.0, 0.0};
    Vec3 b{0.0, 0.0, 7.0};
    CHECK(angle(a, b) == Approx(M_PI / 2.0).margin(LOOSE_EPS));
}

TEST_CASE("angleUnit() agrees with angle() for unit vectors", "[angle]") {
    Vec3 a = Vec3{1.0, 1.0, 0.0}.normalize();
    Vec3 b = Vec3{0.0, 1.0, 1.0}.normalize();
    CHECK(angleUnit(a, b) == Approx(angle(a, b)).margin(LOOSE_EPS));
}

TEST_CASE("angle() is symmetric", "[angle]") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    CHECK(angle(a, b) == Approx(angle(b, a)));
}

// ═══════════════════════════════════════════════════════════════════════════
//  sphericalToSpatial
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("sphericalToSpatial at origin gives (1,0,0)", "[spherical]") {
    Vec3 v = sphericalToSpatial(0.0, 0.0);
    CHECK(v.x == Approx(1.0));
    CHECK(v.y == Approx(0.0).margin(NEAR_ZERO));
    CHECK(v.z == Approx(0.0).margin(NEAR_ZERO));
}

TEST_CASE("sphericalToSpatial at ra=pi/2 gives (0,1,0)", "[spherical]") {
    Vec3 v = sphericalToSpatial(M_PI / 2.0, 0.0);
    CHECK(v.x == Approx(0.0).margin(NEAR_ZERO));
    CHECK(v.y == Approx(1.0));
    CHECK(v.z == Approx(0.0).margin(NEAR_ZERO));
}

TEST_CASE("sphericalToSpatial at de=pi/2 gives north pole", "[spherical]") {
    Vec3 v = sphericalToSpatial(0.0, M_PI / 2.0);
    CHECK(v.x == Approx(0.0).margin(NEAR_ZERO));
    CHECK(v.y == Approx(0.0).margin(NEAR_ZERO));
    CHECK(v.z == Approx(1.0));
}

TEST_CASE("sphericalToSpatial produces unit vectors", "[spherical]") {
    stfloat ra = 1.23;
    stfloat de = 0.45;
    Vec3 v = sphericalToSpatial(ra, de);
    CHECK(v.magnitude() == Approx(1.0).margin(TIGHT_EPS));
}

TEST_CASE("sphericalToSpatial preserves angular distance", "[spherical]") {
    // Two stars 90° apart on the equator
    Vec3 a = sphericalToSpatial(0.0, 0.0);
    Vec3 b = sphericalToSpatial(M_PI / 2.0, 0.0);
    CHECK(angleUnit(a, b) == Approx(M_PI / 2.0).margin(LOOSE_EPS));
}
