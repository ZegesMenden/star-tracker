// Unit tests for Mat3, Quaternion, Attitude, and the QUEST / TRIAD algorithms

#include <catch.hpp>
#include "attitude.hpp"
#include "test_precision.hpp"
#include <cmath>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
//  Mat3 tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Mat3: identity trace is 3", "[mat3]") {
    CHECK(kIdentityMat3.trace() == Approx(3.0));
}

TEST_CASE("Mat3: identity determinant is 1", "[mat3]") {
    CHECK(kIdentityMat3.det() == Approx(1.0));
}

TEST_CASE("Mat3: identity * vec = vec", "[mat3]") {
    Vec3 v{1.0, 2.0, 3.0};
    Vec3 r = kIdentityMat3 * v;
    CHECK(r.x == Approx(v.x));
    CHECK(r.y == Approx(v.y));
    CHECK(r.z == Approx(v.z));
}

TEST_CASE("Mat3: transpose of identity is identity", "[mat3]") {
    Mat3 t = kIdentityMat3.transpose();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            CHECK(t.at(i,j) == Approx(kIdentityMat3.at(i,j)));
}

TEST_CASE("Mat3: transpose reverses indices", "[mat3]") {
    Mat3 m = {1,2,3, 4,5,6, 7,8,9};
    Mat3 t = m.transpose();
    CHECK(t.at(0,1) == Approx(m.at(1,0)));
    CHECK(t.at(0,2) == Approx(m.at(2,0)));
    CHECK(t.at(1,2) == Approx(m.at(2,1)));
}

TEST_CASE("Mat3: determinant of known matrix", "[mat3]") {
    Mat3 m = {1,2,3, 0,1,4, 5,6,0};
    // det = 1*(0-24) - 2*(0-20) + 3*(0-5) = -24+40-15 = 1
    CHECK(m.det() == Approx(1.0));
}

TEST_CASE("Mat3: inverse round-trip", "[mat3]") {
    Mat3 m = {1,2,3, 0,1,4, 5,6,0};
    Mat3 inv = m.inverse();
    Mat3 product = m * inv;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            CHECK(product.at(i,j) == Approx(i == j ? 1.0 : 0.0).margin(LOOSE_EPS));
}

TEST_CASE("Mat3: addition", "[mat3]") {
    Mat3 a = {1,0,0, 0,1,0, 0,0,1};
    Mat3 b = {0,1,0, 1,0,0, 0,0,1};
    Mat3 c = a + b;
    CHECK(c.at(0,0) == Approx(1.0));
    CHECK(c.at(0,1) == Approx(1.0));
    CHECK(c.at(1,0) == Approx(1.0));
    CHECK(c.at(2,2) == Approx(2.0));
}

TEST_CASE("Mat3: scalar multiplication", "[mat3]") {
    Mat3 m = kIdentityMat3 * 3.0;
    CHECK(m.at(0,0) == Approx(3.0));
    CHECK(m.at(1,1) == Approx(3.0));
    CHECK(m.at(0,1) == Approx(0.0));
}

TEST_CASE("Mat3: outerProduct", "[mat3]") {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    Mat3 m = outerProduct(a, b);
    CHECK(m.at(0,1) == Approx(1.0));
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (!(i == 0 && j == 1))
                CHECK(m.at(i,j) == Approx(0.0));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Quaternion tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Quaternion: identity rotation", "[quaternion]") {
    Quaternion q;  // default: (1, 0, 0, 0)
    CHECK(q.real == Approx(1.0));
    CHECK(q.isUnit());
    Vec3 v{1.0, 2.0, 3.0};
    Vec3 r = q.rotate(v);
    CHECK(r.x == Approx(v.x).margin(LOOSE_EPS));
    CHECK(r.y == Approx(v.y).margin(LOOSE_EPS));
    CHECK(r.z == Approx(v.z).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: 90° rotation around z-axis", "[quaternion]") {
    Quaternion q(Vec3{0.0, 0.0, 1.0}, M_PI / 2.0);
    CHECK(q.isUnit());

    Vec3 r = q.rotate(Vec3{1.0, 0.0, 0.0});
    CHECK(r.x == Approx(0.0).margin(LOOSE_EPS));
    CHECK(r.y == Approx(1.0).margin(LOOSE_EPS));
    CHECK(r.z == Approx(0.0).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: 180° rotation around x-axis", "[quaternion]") {
    Quaternion q(Vec3{1.0, 0.0, 0.0}, M_PI);
    CHECK(q.isUnit());

    Vec3 r = q.rotate(Vec3{0.0, 1.0, 0.0});
    CHECK(r.x == Approx(0.0).margin(LOOSE_EPS));
    CHECK(r.y == Approx(-1.0).margin(LOOSE_EPS));
    CHECK(r.z == Approx(0.0).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: conjugate inverts rotation", "[quaternion]") {
    Quaternion q(Vec3{0.0, 1.0, 0.0}, 0.7);
    Vec3 v{1.0, 2.0, 3.0};
    Vec3 rotated = q.rotate(v);
    Vec3 back = q.conjugate().rotate(rotated);
    CHECK(back.x == Approx(v.x).margin(LOOSE_EPS));
    CHECK(back.y == Approx(v.y).margin(LOOSE_EPS));
    CHECK(back.z == Approx(v.z).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: composition = sequential rotation", "[quaternion]") {
    Quaternion q1(Vec3{0.0, 0.0, 1.0}, M_PI / 2.0);  // 90° around z
    Quaternion q2(Vec3{0.0, 1.0, 0.0}, M_PI / 2.0);  // 90° around y

    Vec3 v{1.0, 0.0, 0.0};
    Vec3 sequential = q2.rotate(q1.rotate(v));
    Vec3 composed   = (q2 * q1).rotate(v);

    CHECK(composed.x == Approx(sequential.x).margin(LOOSE_EPS));
    CHECK(composed.y == Approx(sequential.y).margin(LOOSE_EPS));
    CHECK(composed.z == Approx(sequential.z).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: angle() for 90° rotation", "[quaternion]") {
    Quaternion q(Vec3{1.0, 0.0, 0.0}, M_PI / 2.0);
    CHECK(q.angle() == Approx(M_PI / 2.0).margin(LOOSE_EPS));
}

TEST_CASE("Quaternion: smallestAngle for 270° = 90°", "[quaternion]") {
    Quaternion q(Vec3{1.0, 0.0, 0.0}, 3.0 * M_PI / 2.0);
    CHECK(q.smallestAngle() == Approx(M_PI / 2.0).margin(1e-6));
}

TEST_CASE("Quaternion: canonicalize ensures positive real part", "[quaternion]") {
    Quaternion q(-0.5, 0.5, 0.5, 0.5);
    Quaternion c = q.canonicalize();
    CHECK(c.real >= 0);
    CHECK(c.isUnit());
}

TEST_CASE("Quaternion: toSpherical round-trip", "[quaternion]") {
    stfloat ra = 1.2, dec = 0.3, roll = 0.5;
    Quaternion q = sphericalToQuaternion(ra, dec, roll);
    CHECK(q.isUnit());
    EulerAngles e = q.toSpherical();
    CHECK(e.ra == Approx(ra).margin(1e-6));
    CHECK(e.de == Approx(dec).margin(1e-6));
    CHECK(e.roll == Approx(roll).margin(1e-6));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Attitude wrapper
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Attitude: default is unknown", "[attitude]") {
    Attitude att;
    CHECK(!att.isKnown());
}

TEST_CASE("Attitude: from quaternion is known", "[attitude]") {
    Attitude att{Quaternion{}};
    CHECK(att.isKnown());
}

TEST_CASE("Attitude: rotate delegates to quaternion", "[attitude]") {
    Quaternion q(Vec3{0, 0, 1}, M_PI / 2.0);
    Attitude att(q);
    Vec3 r = att.rotate(Vec3{1.0, 0.0, 0.0});
    CHECK(r.x == Approx(0.0).margin(LOOSE_EPS));
    CHECK(r.y == Approx(1.0).margin(LOOSE_EPS));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Test helper: create a synthetic identified-star scene
// ═══════════════════════════════════════════════════════════════════════════

/// Build a catalog, apply a known rotation, project into camera, and produce
/// matched (Stars, StarIdentifiers) pairs.
struct SyntheticScene {
    Catalog        catalog;
    Camera         camera;
    Quaternion     trueAttitude;
    Stars          stars;
    StarIdentifiers ids;

    SyntheticScene(stfloat ra, stfloat dec, stfloat roll)
        : camera(fovToFocalLength(40.0 * M_PI / 180.0, 1024.0), 1024, 1024),
          trueAttitude(sphericalToQuaternion(ra, dec, roll))
    {
        // Catalog stars spread around ra≈0, dec≈0
        catalog.push_back(CatalogStar(0.00,  0.00,  100, 0));
        catalog.push_back(CatalogStar(0.05,  0.02,  100, 1));
        catalog.push_back(CatalogStar(0.10,  0.00,  100, 2));
        catalog.push_back(CatalogStar(0.03, -0.05,  100, 3));
        catalog.push_back(CatalogStar(0.08,  0.06,  100, 4));
        catalog.push_back(CatalogStar(0.12, -0.03,  100, 5));
        catalog.push_back(CatalogStar(0.06,  0.08,  100, 6));
        catalog.push_back(CatalogStar(0.15,  0.01,  100, 7));
        catalog.push_back(CatalogStar(0.02,  0.10,  100, 8));
        catalog.push_back(CatalogStar(0.11,  0.07,  100, 9));

        // Project each catalog star through the attitude into the camera
        for (int ci = 0; ci < (int)catalog.size(); ci++) {
            Vec3 bodyVec = trueAttitude.rotate(catalog[ci].spatial);
            if (bodyVec.x <= 0) continue;  // behind camera
            Vec2 px = camera.spatialToCamera(bodyVec);
            if (!camera.inSensor(px)) continue;
            stars.push_back(Star(px.x, px.y, 2.0));
            ids.push_back(StarIdentifier((int)stars.size() - 1, ci));
        }
    }
};

// ═══════════════════════════════════════════════════════════════════════════
//  QUEST algorithm tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("QUEST: recovers identity attitude", "[quest]") {
    // Camera pointing at ra=0, dec=0, roll=0 (identity rotation equivalent)
    // The trueAttitude quaternion should be close to identity.
    // Use catalog stars that lie along the +x axis so they're visible.
    Camera cam(fovToFocalLength(40.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    Catalog cat;
    cat.push_back(CatalogStar(0.00,  0.00,  100, 0));
    cat.push_back(CatalogStar(0.05,  0.02,  100, 1));
    cat.push_back(CatalogStar(0.10,  0.01,  100, 2));
    cat.push_back(CatalogStar(0.03, -0.04,  100, 3));

    // Project with identity rotation (catalog spatial IS body spatial)
    Stars stars;
    StarIdentifiers ids;
    for (int ci = 0; ci < (int)cat.size(); ci++) {
        const Vec3 &sp = cat[ci].spatial;
        if (sp.x <= 0) continue;
        Vec2 px = cam.spatialToCamera(sp);
        if (!cam.inSensor(px)) continue;
        stars.push_back(Star(px.x, px.y, 2.0));
        ids.push_back(StarIdentifier((int)stars.size() - 1, ci));
    }
    REQUIRE(ids.size() >= 2);

    QuestAlgorithm quest;
    Attitude att = quest.estimate(cam, stars, cat, ids);
    REQUIRE(att.isKnown());

    Quaternion q = att.getQuaternion().canonicalize();
    // Should be close to identity: smallestAngle ≈ 0
    CHECK(q.smallestAngle() < 0.001);
}

TEST_CASE("QUEST: recovers known rotation (small ra offset)", "[quest]") {
    SyntheticScene scene(0.05, 0.0, 0.0);
    REQUIRE(scene.ids.size() >= 4);

    QuestAlgorithm quest;
    Attitude att = quest.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    REQUIRE(att.isKnown());

    // Verify the recovered attitude matches the true one by comparing
    // how each catalog vector is rotated.
    for (const auto &id : scene.ids) {
        Vec3 expectedBody = scene.trueAttitude.rotate(scene.catalog[id.catalogIndex].spatial);
        Vec3 recoveredBody = att.getQuaternion().rotate(scene.catalog[id.catalogIndex].spatial);
        CHECK(angleUnit(expectedBody.normalize(), recoveredBody.normalize()) < 0.001);
    }
}

TEST_CASE("QUEST: recovers known rotation (ra + dec + roll)", "[quest]") {
    SyntheticScene scene(0.03, 0.02, 0.1);
    REQUIRE(scene.ids.size() >= 4);

    QuestAlgorithm quest;
    Attitude att = quest.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    REQUIRE(att.isKnown());

    Quaternion error = att.getQuaternion() * scene.trueAttitude.conjugate();
    CHECK(error.canonicalize().smallestAngle() < 0.001);
}

TEST_CASE("QUEST: more stars improves accuracy", "[quest]") {
    // Scene with all 10 stars
    SyntheticScene full(0.04, -0.02, 0.05);
    REQUIRE(full.ids.size() >= 6);

    // Scene with only the first 3
    Stars fewStars;
    StarIdentifiers fewIds;
    for (int i = 0; i < 3 && i < (int)full.ids.size(); i++) {
        fewStars.push_back(full.stars[full.ids[i].starIndex]);
        fewIds.push_back(StarIdentifier(i, full.ids[i].catalogIndex));
    }

    QuestAlgorithm quest;
    Attitude attFull = quest.estimate(full.camera, full.stars, full.catalog, full.ids);
    Attitude attFew  = quest.estimate(full.camera, fewStars, full.catalog, fewIds);

    REQUIRE(attFull.isKnown());
    REQUIRE(attFew.isKnown());

    // Both should be close, but full should be at least as good
    Quaternion errFull = attFull.getQuaternion() * full.trueAttitude.conjugate();
    Quaternion errFew  = attFew.getQuaternion()  * full.trueAttitude.conjugate();
    CHECK(errFull.canonicalize().smallestAngle() <= errFew.canonicalize().smallestAngle() + 0.01);
}

TEST_CASE("QUEST: too few stars returns unknown", "[quest]") {
    Camera cam(200.0, 512, 512);
    Catalog cat;
    cat.push_back(CatalogStar(0.0, 0.0, 100, 0));

    Stars stars;
    stars.push_back(Star(256.0, 256.0, 2.0));
    StarIdentifiers ids;
    ids.push_back(StarIdentifier(0, 0));

    QuestAlgorithm quest;
    Attitude att = quest.estimate(cam, stars, cat, ids);
    CHECK(!att.isKnown());
}

TEST_CASE("QUEST: result quaternion is unit", "[quest]") {
    SyntheticScene scene(0.06, 0.01, 0.2);
    REQUIRE(scene.ids.size() >= 2);

    QuestAlgorithm quest;
    Attitude att = quest.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    REQUIRE(att.isKnown());
    CHECK(att.getQuaternion().isUnit());
}

// ═══════════════════════════════════════════════════════════════════════════
//  TRIAD algorithm tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("TRIAD: recovers known rotation", "[triad]") {
    SyntheticScene scene(0.05, 0.02, 0.1);
    REQUIRE(scene.ids.size() >= 2);

    TriadAlgorithm triad;
    Attitude att = triad.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    REQUIRE(att.isKnown());

    Quaternion error = att.getQuaternion() * scene.trueAttitude.conjugate();
    // TRIAD is less accurate, but should still be within ~1° for clean data
    CHECK(error.canonicalize().smallestAngle() < 0.02);
}

TEST_CASE("TRIAD: too few stars returns unknown", "[triad]") {
    Camera cam(200.0, 512, 512);
    Catalog cat;
    cat.push_back(CatalogStar(0.0, 0.0, 100, 0));
    Stars stars;
    stars.push_back(Star(256.0, 256.0, 2.0));
    StarIdentifiers ids;
    ids.push_back(StarIdentifier(0, 0));

    TriadAlgorithm triad;
    Attitude att = triad.estimate(cam, stars, cat, ids);
    CHECK(!att.isKnown());
}

TEST_CASE("TRIAD: result quaternion is unit", "[triad]") {
    SyntheticScene scene(0.03, 0.01, 0.0);
    REQUIRE(scene.ids.size() >= 2);

    TriadAlgorithm triad;
    Attitude att = triad.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    REQUIRE(att.isKnown());
    CHECK(att.getQuaternion().isUnit(0.001));
}

// ═══════════════════════════════════════════════════════════════════════════
//  QUEST vs TRIAD comparison
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("QUEST and TRIAD agree on clean data", "[quest][triad]") {
    SyntheticScene scene(0.04, 0.03, 0.08);
    REQUIRE(scene.ids.size() >= 4);

    QuestAlgorithm quest;
    TriadAlgorithm triad;

    Attitude attQ = quest.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);
    Attitude attT = triad.estimate(scene.camera, scene.stars, scene.catalog, scene.ids);

    REQUIRE(attQ.isKnown());
    REQUIRE(attT.isKnown());

    // Both should produce similar results
    Quaternion diff = attQ.getQuaternion() * attT.getQuaternion().conjugate();
    CHECK(diff.canonicalize().smallestAngle() < 0.02);
}
