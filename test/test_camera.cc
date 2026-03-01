// Unit tests for Camera (camera.hpp / camera.cc)

#include <catch.hpp>
#include "camera.hpp"
#include "test_precision.hpp"
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════
//  Projection round-trip
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: pixel -> spatial -> pixel round-trip", "[camera]") {
    Camera cam(100.0, 512, 1024);

    stfloat px = GENERATE(142.0, 90.0, 256.0, 400.0);
    stfloat py = GENERATE(18.0, 512.0, 100.0, 800.0);

    Vec3 spatial = cam.cameraToSpatial({px, py});
    Vec2 result  = cam.spatialToCamera(spatial);

    CHECK(result.x == Approx(px).margin(1e-6));
    CHECK(result.y == Approx(py).margin(1e-6));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Center pixel → boresight
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: center pixel maps to boresight (y=0, z=0)", "[camera]") {
    Camera cam(100.0, 512, 1024);

    Vec3 spatial = cam.cameraToSpatial({256.0, 512.0});  // center
    CHECK(spatial.y == Approx(0.0).margin(TIGHT_EPS));
    CHECK(spatial.z == Approx(0.0).margin(TIGHT_EPS));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Known explicit projection
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: explicit pixel -> spatial direction", "[camera]") {
    Camera cam(100.0, 512, 1024);

    // Pixel at (300, 728)
    Vec3 spatial = cam.cameraToSpatial({300.0, 728.0}).normalize();

    // Expected (ish): x positive, y negative (pixel right of center), z negative (pixel below center)
    CHECK(spatial.x > 0.0);
    CHECK(spatial.y < 0.0);
    CHECK(spatial.z < 0.0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Angle between pixels
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: angle across full horizontal FOV", "[camera]") {
    Camera cam(128.0, 256, 256);

    Vec3 s1 = cam.cameraToSpatial({0.0, 128.0});
    Vec3 s2 = cam.cameraToSpatial({256.0, 128.0});
    stfloat a = angle(s1, s2);
    CHECK(a == Approx(M_PI / 2.0).margin(1e-6));  // 90° with focal=128, width=256
}

TEST_CASE("Camera: angle across full vertical FOV", "[camera]") {
    Camera cam(128.0, 256, 256);

    Vec3 s1 = cam.cameraToSpatial({128.0, 0.0});
    Vec3 s2 = cam.cameraToSpatial({128.0, 256.0});
    stfloat a = angle(s1, s2);
    CHECK(a == Approx(M_PI / 2.0).margin(1e-6));
}

TEST_CASE("Camera: diagonal angle larger than horizontal", "[camera]") {
    Camera cam(200.0, 512, 512);
    Vec3 topLeft  = cam.cameraToSpatial({0.0, 0.0});
    Vec3 botRight = cam.cameraToSpatial({512.0, 512.0});
    Vec3 midLeft  = cam.cameraToSpatial({0.0, 256.0});
    Vec3 midRight = cam.cameraToSpatial({512.0, 256.0});

    stfloat diag = angle(topLeft, botRight);
    stfloat horiz = angle(midLeft, midRight);
    CHECK(diag > horiz);
}

// ═══════════════════════════════════════════════════════════════════════════
//  inSensor
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: inSensor boundary check", "[camera]") {
    Camera cam(100.0, 640, 480);

    CHECK(cam.inSensor({0.0, 0.0}));
    CHECK(cam.inSensor({320.0, 240.0}));
    CHECK(cam.inSensor({640.0, 480.0}));
    CHECK(!cam.inSensor({-1.0, 0.0}));
    CHECK(!cam.inSensor({0.0, -1.0}));
    CHECK(!cam.inSensor({641.0, 0.0}));
    CHECK(!cam.inSensor({0.0, 481.0}));
}

// ═══════════════════════════════════════════════════════════════════════════
//  FOV calculation
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: fov() is consistent with fovToFocalLength()", "[camera]") {
    Camera cam(200.0, 1024, 768);

    stfloat computedFov = cam.fov();
    stfloat computedFL  = fovToFocalLength(computedFov, 1024.0);
    CHECK(computedFL == Approx(200.0).margin(1e-6));
}

TEST_CASE("fovToFocalLength: 90° FOV on 256px sensor gives 128px focal", "[camera]") {
    stfloat fl = fovToFocalLength(M_PI / 2.0, 256.0);
    CHECK(fl == Approx(128.0).margin(1e-6));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Custom principal point
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Camera: custom principal point affects projection", "[camera]") {
    Camera cam1(100.0, 512, 512);                     // center at (256, 256)
    Camera cam2(100.0, 200.0, 200.0, 512, 512);       // center at (200, 200)

    Vec3 s1 = cam1.cameraToSpatial({256.0, 256.0});   // boresight for cam1
    Vec3 s2 = cam2.cameraToSpatial({200.0, 200.0});   // boresight for cam2

    // Both boresights should map to pure x direction
    CHECK(s1.y == Approx(0.0).margin(TIGHT_EPS));
    CHECK(s1.z == Approx(0.0).margin(TIGHT_EPS));
    CHECK(s2.y == Approx(0.0).margin(TIGHT_EPS));
    CHECK(s2.z == Approx(0.0).margin(TIGHT_EPS));
}
