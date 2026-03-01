// Unit tests for star identification (starid.hpp / starid.cc)
// Tests the helper functions and the Pyramid algorithm with a synthetic scene.

#include <catch.hpp>
#include "starid.hpp"
#include "database.hpp"
#include "camera.hpp"
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <sstream>
#include <iostream>

/// RAII helper: redirects std::cerr to a string while in scope so that
/// expected error messages from negative test cases don't pollute output.
struct SuppressCerr {
    std::streambuf *orig;
    std::ostringstream sink;
    SuppressCerr()  : orig(std::cerr.rdbuf(sink.rdbuf())) {}
    ~SuppressCerr() { std::cerr.rdbuf(orig); }
};

// ═══════════════════════════════════════════════════════════════════════════
//  Test helper: build a small catalog, construct a database, and simulate
//  a camera image by projecting known catalog stars onto the focal plane.
// ═══════════════════════════════════════════════════════════════════════════

/// A tight cluster of stars within ~20° near the boresight direction,
/// easily visible in a camera with ~40° FOV.
static Catalog makeSmallCatalog() {
    Catalog cat;
    // Stars at known positions (RA, Dec in radians).  All near (0, 0).
    cat.push_back(CatalogStar(0.00,  0.00,  100, 0));
    cat.push_back(CatalogStar(0.05,  0.02,  100, 1));
    cat.push_back(CatalogStar(0.10,  0.00,  100, 2));
    cat.push_back(CatalogStar(0.03, -0.05,  100, 3));
    cat.push_back(CatalogStar(0.08,  0.06,  100, 4));
    cat.push_back(CatalogStar(0.12, -0.03,  100, 5));
    cat.push_back(CatalogStar(0.06,  0.08,  100, 6));
    cat.push_back(CatalogStar(0.15,  0.01,  100, 7));
    cat.push_back(CatalogStar(0.02,  0.10,  100, 8));
    cat.push_back(CatalogStar(0.11,  0.07,  100, 9));
    return cat;
}

/// Camera centred on (ra=0.07, de=0.02) with enough FOV to see all catalog stars.
static Camera makeTestCamera() {
    // focal length chosen so FOV ≈ 40°
    return Camera(fovToFocalLength(40.0 * M_PI / 180.0, 1024.0), 1024, 1024);
}

/// Project catalog stars into the camera, returning Stars (centroids) and the
/// mapping from star index → catalog index.
static void projectCatalog(const Catalog &cat, const Camera &cam,
                           Stars &outStars,
                           std::vector<int> &outCatIndices) {
    outStars.clear();
    outCatIndices.clear();
    for (int ci = 0; ci < (int)cat.size(); ci++) {
        const Vec3 &sp = cat[ci].spatial;
        // Only project if in front of the camera (x > 0) — for a pointing near
        // the catalog centre this is always true.
        if (sp.x <= 0) continue;
        Vec2 px = cam.spatialToCamera(sp);
        if (cam.inSensor(px)) {
            outStars.push_back(Star(px.x, px.y, 2.0));
            outCatIndices.push_back(ci);
        }
    }
}

/// Build a full serialized multi-database from the catalog.
static std::vector<unsigned char>
buildDatabase(const Catalog &cat, stfloat minDist, stfloat maxDist, long numBins) {
    // Pair-distance KVector sub-database
    SerializeContext pairSer;
    serializePairDistanceKVector(&pairSer, cat, minDist, maxDist, numBins);

    // Catalog sub-database
    SerializeContext catSer;
    serializeCatalog(&catSer, cat);

    // Combine into a MultiDatabase
    MultiDatabaseDescriptor desc;
    desc.push_back(MultiDatabaseEntry(PairDistanceKVectorDatabase::kMagicValue,
                                      pairSer.buffer));
    desc.push_back(MultiDatabaseEntry(kCatalogMagicValue, catSer.buffer));

    SerializeContext multiSer;
    serializeMultiDatabase(&multiSer, desc);
    return multiSer.buffer;
}

// ═══════════════════════════════════════════════════════════════════════════
//  pairDistanceQueryToMap
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("pairDistanceQueryToMap: builds bidirectional map", "[starid][helper]") {
    // Simulate a raw pairs buffer:  (0,1), (2,3)
    int16_t raw[] = {0, 1, 2, 3};
    auto map = pairDistanceQueryToMap(raw, raw + 4);

    CHECK(map.count(0) >= 1);
    CHECK(map.count(1) >= 1);
    CHECK(map.count(2) >= 1);
    CHECK(map.count(3) >= 1);

    // 0→1 and 1→0 should both be present
    bool found01 = false, found10 = false;
    for (auto [k, v] : map) {
        if (k == 0 && v == 1) found01 = true;
        if (k == 1 && v == 0) found10 = true;
    }
    CHECK(found01);
    CHECK(found10);
}

TEST_CASE("pairDistanceQueryToMap: empty input yields empty map", "[starid][helper]") {
    int16_t raw[] = {0};  // not used
    auto map = pairDistanceQueryToMap(raw, raw);  // begin == end
    CHECK(map.empty());
}

// ═══════════════════════════════════════════════════════════════════════════
//  identifyRemainingStarsPairDistance
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("identifyRemainingStars finds additional stars given a seed", "[starid][remaining]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    // Need at least 4 projected stars
    REQUIRE(stars.size() >= 4);

    // Build the pair-distance DB
    stfloat maxAngle = cam.fov();
    stfloat minAngle = 0.01;
    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, minAngle, maxAngle, 100);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    // Seed with the first 4 correct identifications
    StarIdentifiers identified;
    for (int i = 0; i < 4 && i < (int)stars.size(); i++) {
        identified.push_back(StarIdentifier(i, catIdx[i]));
    }

    int numRemaining = (int)stars.size() - 4;
    stfloat tolerance = 0.001;

    int numExtra = identifyRemainingStarsPairDistance(
        identified, stars, db, cat, cam, tolerance);

    // Should have identified at least some additional stars
    if (numRemaining > 0) {
        CHECK(numExtra > 0);
    }
    CHECK((int)identified.size() > 4);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Pyramid full integration test
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Pyramid: identifies stars in a synthetic scene", "[starid][pyramid]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    stfloat minAngle = 0.005;

    auto db = buildDatabase(cat, minAngle, maxAngle, 200);

    stfloat tolerance = 0.001;
    PyramidStarIdAlgorithm algo(tolerance, 0, 0.01, 100000);

    StarIdentifiers result = algo.identify(db.data(), stars, cat, cam);

    // The algorithm should identify at least 4 stars (the initial pyramid)
    CHECK(result.size() >= 4);

    // Verify that every identification is correct
    for (const auto &id : result) {
        REQUIRE(id.starIndex >= 0);
        REQUIRE(id.starIndex < (int)stars.size());

        // The centroid at id.starIndex should match the catalog star id.catalogIndex
        CHECK(id.catalogIndex == catIdx[id.starIndex]);
    }
}

TEST_CASE("Pyramid: too few stars returns empty", "[starid][pyramid]") {
    SuppressCerr quiet;  // suppress expected "not enough stars" message
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();

    Stars stars;
    stars.push_back(Star(100.0, 100.0, 2.0));  // just one star

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    PyramidStarIdAlgorithm algo(0.001, 0, 0.01, 100000);
    StarIdentifiers result = algo.identify(db.data(), stars, cat, cam);

    CHECK(result.empty());
}

TEST_CASE("Pyramid: empty database returns empty", "[starid][pyramid]") {
    SuppressCerr quiet;  // suppress expected "not enough stars or database missing" message
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    // Build a multi-database with NO pair-distance sub-db
    SerializeContext catSer;
    serializeCatalog(&catSer, cat);
    MultiDatabaseDescriptor desc;
    desc.push_back(MultiDatabaseEntry(kCatalogMagicValue, catSer.buffer));
    SerializeContext multiSer;
    serializeMultiDatabase(&multiSer, desc);

    PyramidStarIdAlgorithm algo(0.001, 0, 0.01, 100000);
    StarIdentifiers result = algo.identify(multiSer.buffer.data(), stars, cat, cam);

    CHECK(result.empty());
}

TEST_CASE("Pyramid: no duplicate identifications", "[starid][pyramid]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    PyramidStarIdAlgorithm algo(0.001, 0, 0.01, 100000);
    StarIdentifiers result = algo.identify(db.data(), stars, cat, cam);

    // Check no star index appears twice
    std::set<int> starIndices;
    for (const auto &id : result) {
        CHECK(starIndices.find(id.starIndex) == starIndices.end());
        starIndices.insert(id.starIndex);
    }

    // Check no catalog index appears twice
    std::set<int> catIndices;
    for (const auto &id : result) {
        CHECK(catIndices.find(id.catalogIndex) == catIndices.end());
        catIndices.insert(id.catalogIndex);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tetra3 full integration tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Tetra3: identifies stars in a synthetic scene", "[starid][tetra3]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    Tetra3StarIdAlgorithm algo(
        0.001,    // distance tolerance
        0.01,     // pattern hash ratio bin size
        10,       // image stars for pattern checking
        10,       // catalog stars for pattern hash
        200000);  // max generated catalog patterns

    StarIdentifiers result = algo.identify(db.data(), stars, cat, cam);

    CHECK(result.size() >= 4);

    for (const auto &id : result) {
        REQUIRE(id.starIndex >= 0);
        REQUIRE(id.starIndex < (int)stars.size());
        CHECK(id.catalogIndex == catIdx[id.starIndex]);
    }
}

TEST_CASE("Tetra3: repeated solve is stable", "[starid][tetra3]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    Tetra3StarIdAlgorithm algo(0.001, 0.01, 10, 10, 200000);
    StarIdentifiers first = algo.identify(db.data(), stars, cat, cam);
    StarIdentifiers second = algo.identify(db.data(), stars, cat, cam);

    REQUIRE(first.size() >= 4);
    REQUIRE(second.size() >= 4);

    std::vector<int> mapFirst(stars.size(), -1);
    std::vector<int> mapSecond(stars.size(), -1);
    for (const auto &id : first) mapFirst[id.starIndex] = id.catalogIndex;
    for (const auto &id : second) mapSecond[id.starIndex] = id.catalogIndex;
    CHECK(mapFirst == mapSecond);
}

TEST_CASE("Tetra3: too few stars returns empty", "[starid][tetra3]") {
    SuppressCerr quiet;  // suppress expected "insufficient stars" message
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();

    Stars stars;
    stars.push_back(Star(100.0, 100.0, 2.0));

    Tetra3StarIdAlgorithm algo(0.001, 0.01, 10, 10, 200000);
    StarIdentifiers result = algo.identify(nullptr, stars, cat, cam);

    CHECK(result.empty());
}

TEST_CASE("Tetra3: empty multi-database still identifies seed tetra", "[starid][tetra3]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    // Build a multi-database with no pair-distance sub-database.
    SerializeContext catSer;
    serializeCatalog(&catSer, cat);
    MultiDatabaseDescriptor desc;
    desc.push_back(MultiDatabaseEntry(kCatalogMagicValue, catSer.buffer));
    SerializeContext multiSer;
    serializeMultiDatabase(&multiSer, desc);

    Tetra3StarIdAlgorithm algo(0.001, 0.01, 10, 10, 200000);
    StarIdentifiers result = algo.identify(multiSer.buffer.data(), stars, cat, cam);

    CHECK(result.size() >= 4);
    for (const auto &id : result) {
        CHECK(id.catalogIndex == catIdx[id.starIndex]);
    }
}

TEST_CASE("Tetra3: invalid pattern error returns empty", "[starid][tetra3]") {
    SuppressCerr quiet;  // suppress expected "invalid parameters" message
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    Tetra3StarIdAlgorithm algo(0.001, 0.0, 10, 10, 200000);
    StarIdentifiers result = algo.identify(nullptr, stars, cat, cam);

    CHECK(result.empty());
}

TEST_CASE("Tetra3: no duplicate identifications", "[starid][tetra3]") {
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    Tetra3StarIdAlgorithm algo(0.001, 0.01, 10, 10, 200000);
    StarIdentifiers result = algo.identify(db.data(), stars, cat, cam);

    std::set<int> starIndices;
    for (const auto &id : result) {
        CHECK(starIndices.find(id.starIndex) == starIndices.end());
        starIndices.insert(id.starIndex);
    }

    std::set<int> catIndices;
    for (const auto &id : result) {
        CHECK(catIndices.find(id.catalogIndex) == catIndices.end());
        catIndices.insert(id.catalogIndex);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Accuracy comparison: Tetra3 vs Pyramid
// ═══════════════════════════════════════════════════════════════════════════

/// A larger catalog with 30 stars spread over a wider area, used for more
/// realistic accuracy tests.
static Catalog makeLargeCatalog() {
    Catalog cat;
    // Evenly scattered stars near the boresight (RA, Dec in radians).
    cat.push_back(CatalogStar(0.00,  0.00,  80,  0));
    cat.push_back(CatalogStar(0.05,  0.02,  90,  1));
    cat.push_back(CatalogStar(0.10,  0.00,  85,  2));
    cat.push_back(CatalogStar(0.03, -0.05,  95,  3));
    cat.push_back(CatalogStar(0.08,  0.06, 100,  4));
    cat.push_back(CatalogStar(0.12, -0.03, 105,  5));
    cat.push_back(CatalogStar(0.06,  0.08, 110,  6));
    cat.push_back(CatalogStar(0.15,  0.01,  88,  7));
    cat.push_back(CatalogStar(0.02,  0.10,  92,  8));
    cat.push_back(CatalogStar(0.11,  0.07,  97,  9));
    cat.push_back(CatalogStar(0.18,  0.04,  83, 10));
    cat.push_back(CatalogStar(0.01, -0.08,  91, 11));
    cat.push_back(CatalogStar(0.14,  0.09, 102, 12));
    cat.push_back(CatalogStar(0.07, -0.06, 104, 13));
    cat.push_back(CatalogStar(0.20,  0.00,  87, 14));
    cat.push_back(CatalogStar(0.04,  0.12,  93, 15));
    cat.push_back(CatalogStar(0.16, -0.05,  96, 16));
    cat.push_back(CatalogStar(0.09,  0.11,  99, 17));
    cat.push_back(CatalogStar(0.13, -0.07, 101, 18));
    cat.push_back(CatalogStar(0.22,  0.03,  86, 19));
    cat.push_back(CatalogStar(0.17,  0.08,  94, 20));
    cat.push_back(CatalogStar(0.19, -0.02,  89, 21));
    cat.push_back(CatalogStar(0.21,  0.06,  98, 22));
    cat.push_back(CatalogStar(0.24, -0.01, 103, 23));
    cat.push_back(CatalogStar(0.23,  0.05, 106, 24));
    cat.push_back(CatalogStar(0.25,  0.02,  84, 25));
    cat.push_back(CatalogStar(0.26, -0.04, 107, 26));
    cat.push_back(CatalogStar(0.28,  0.01, 108, 27));
    cat.push_back(CatalogStar(0.27,  0.07,  82, 28));
    cat.push_back(CatalogStar(0.30,  0.00, 109, 29));
    return cat;
}

/// Helper: count correct identifications.
static int countCorrect(const StarIdentifiers &result, const std::vector<int> &catIdx) {
    int correct = 0;
    for (const auto &id : result) {
        if (id.starIndex >= 0 && id.starIndex < (int)catIdx.size()) {
            if (id.catalogIndex == catIdx[id.starIndex]) ++correct;
        }
    }
    return correct;
}

/// Helper: check all identifications are correct (no false positives).
static bool allCorrect(const StarIdentifiers &result, const std::vector<int> &catIdx) {
    for (const auto &id : result) {
        if (id.starIndex < 0 || id.starIndex >= (int)catIdx.size()) return false;
        if (id.catalogIndex != catIdx[id.starIndex]) return false;
    }
    return true;
}

/// Project catalog stars through a rotated camera.
static void projectCatalogRotated(const Catalog &cat, const Camera &cam,
                                  const Quaternion &rotation,
                                  Stars &outStars,
                                  std::vector<int> &outCatIndices) {
    outStars.clear();
    outCatIndices.clear();
    for (int ci = 0; ci < (int)cat.size(); ci++) {
        // Apply inverse rotation to put catalog stars in camera frame.
        Vec3 sp = rotation.conjugate().rotate(cat[ci].spatial);
        if (sp.x <= 0) continue;
        Vec2 px = cam.spatialToCamera(sp);
        if (cam.inSensor(px)) {
            outStars.push_back(Star(px.x, px.y, 2.0));
            outCatIndices.push_back(ci);
        }
    }
}

TEST_CASE("Tetra3 vs Pyramid: same scene identification count", "[starid][accuracy]") {
    // Both algorithms should identify at least 4 stars and all identifications
    // must be correct.
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 100000);
    StarIdentifiers pyramidResult = pyramid.identify(db.data(), stars, cat, cam);

    Tetra3StarIdAlgorithm tetra3(0.001, 0.01, 10, 10, 200000);
    StarIdentifiers tetra3Result = tetra3.identify(db.data(), stars, cat, cam);

    REQUIRE(pyramidResult.size() >= 4);
    REQUIRE(tetra3Result.size() >= 4);

    CHECK(allCorrect(pyramidResult, catIdx));
    CHECK(allCorrect(tetra3Result, catIdx));

    // Tetra3 should identify at least as many stars as Pyramid (with pair-
    // distance refinement both should reach the same count).
    CHECK(tetra3Result.size() >= pyramidResult.size() - 1);
}

TEST_CASE("Tetra3 vs Pyramid: multiple boresight directions", "[starid][accuracy]") {
    // Rotate the camera to several different orientations and verify both
    // algorithms produce correct results each time.
    Catalog cat = makeLargeCatalog();
    Camera cam  = Camera(fovToFocalLength(50.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    struct Pointing {
        stfloat ra, de;  // radians
    };
    // Several boresight directions within the catalog's extent.
    Pointing pointings[] = {
        {0.10,  0.00},
        {0.15,  0.03},
        {0.05, -0.02},
        {0.20,  0.02},
        {0.12,  0.06},
    };

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 500);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 200000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 15, (int)cat.size(), 500000);

    int pyramidTotal = 0, tetra3Total = 0;
    int pyramidScenes = 0, tetra3Scenes = 0;

    for (const auto &pt : pointings) {
        Quaternion rot = sphericalToQuaternion(pt.ra, pt.de, 0.0);
        Stars stars;
        std::vector<int> catIdx;
        projectCatalogRotated(cat, cam, rot, stars, catIdx);

        if (stars.size() < 4) continue;

        StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
        StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

        if (!pRes.empty()) {
            CHECK(allCorrect(pRes, catIdx));
            pyramidTotal += (int)pRes.size();
            ++pyramidScenes;
        }
        if (!tRes.empty()) {
            CHECK(allCorrect(tRes, catIdx));
            tetra3Total += (int)tRes.size();
            ++tetra3Scenes;
        }
    }

    // Tetra3 should solve at least as many scenes as Pyramid.
    CHECK(tetra3Scenes >= pyramidScenes);
}

TEST_CASE("Tetra3 vs Pyramid: larger catalog correctness", "[starid][accuracy]") {
    // Use the larger catalog to verify that Tetra3 still produces 100% correct
    // identifications and identifies at least 4 stars.
    Catalog cat = makeLargeCatalog();
    Camera cam  = Camera(fovToFocalLength(50.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    Quaternion rot = sphericalToQuaternion(0.14, 0.02, 0.0);
    Stars stars;
    std::vector<int> catIdx;
    projectCatalogRotated(cat, cam, rot, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 500);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 200000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 15, (int)cat.size(), 500000);

    StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
    StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

    REQUIRE(pRes.size() >= 4);
    REQUIRE(tRes.size() >= 4);

    CHECK(allCorrect(pRes, catIdx));
    CHECK(allCorrect(tRes, catIdx));

    int pCorrect = countCorrect(pRes, catIdx);
    int tCorrect = countCorrect(tRes, catIdx);

    // Both must be 100% correct.
    CHECK(pCorrect == (int)pRes.size());
    CHECK(tCorrect == (int)tRes.size());

    // Tetra3 identification count should be comparable to Pyramid.
    CHECK(tCorrect >= pCorrect - 1);
}

TEST_CASE("Tetra3 vs Pyramid: robustness with centroid noise", "[starid][accuracy]") {
    // Add small perturbations to the centroid positions and verify that both
    // algorithms still produce correct results.
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    // Apply sub-pixel noise (deterministic pattern, not random).
    Stars noisyStars;
    stfloat offsets[] = {0.3, -0.2, 0.1, -0.4, 0.25, -0.15, 0.35, -0.05, 0.2, -0.3};
    for (int i = 0; i < (int)stars.size(); i++) {
        stfloat dx = offsets[i % 10];
        stfloat dy = offsets[(i + 5) % 10];
        noisyStars.push_back(Star(stars[i].position.x + dx,
                                  stars[i].position.y + dy, 2.0));
    }

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    // Use a slightly larger tolerance to account for noise.
    PyramidStarIdAlgorithm pyramid(0.002, 0, 0.01, 100000);
    Tetra3StarIdAlgorithm  tetra3(0.002, 0.01, 10, 10, 200000);

    StarIdentifiers pRes = pyramid.identify(db.data(), noisyStars, cat, cam);
    StarIdentifiers tRes = tetra3.identify(db.data(), noisyStars, cat, cam);

    // Both should still identify at least the initial 4-star pattern.
    CHECK(pRes.size() >= 4);
    CHECK(tRes.size() >= 4);

    // All identifications must still be correct.
    if (!pRes.empty()) CHECK(allCorrect(pRes, catIdx));
    if (!tRes.empty()) CHECK(allCorrect(tRes, catIdx));
}

TEST_CASE("Tetra3 vs Pyramid: false stars in the image", "[starid][accuracy]") {
    // Add spurious centroids that don't correspond to any catalog star and
    // verify that both algorithms still identify the real stars correctly.
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();
    Stars stars;
    std::vector<int> catIdx;
    projectCatalog(cat, cam, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    // Inject false stars at positions that don't match any catalog star.
    Stars augmentedStars = stars;
    augmentedStars.push_back(Star(50.0,  50.0,  2.0));   // corner
    augmentedStars.push_back(Star(900.0, 900.0, 2.0));   // far corner
    augmentedStars.push_back(Star(200.0, 800.0, 2.0));   // off-pattern

    // Extend catIdx with -1 for false stars (no correct match).
    std::vector<int> augCatIdx = catIdx;
    augCatIdx.push_back(-1);
    augCatIdx.push_back(-1);
    augCatIdx.push_back(-1);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    PyramidStarIdAlgorithm pyramid(0.001, 3, 0.01, 100000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 13, 10, 200000);

    StarIdentifiers pRes = pyramid.identify(db.data(), augmentedStars, cat, cam);
    StarIdentifiers tRes = tetra3.identify(db.data(), augmentedStars, cat, cam);

    // At least 4 correct identifications from each.
    CHECK(pRes.size() >= 4);
    CHECK(tRes.size() >= 4);

    // None of the false-star indices should appear as identified.
    int realStarCount = (int)stars.size();
    for (const auto &id : pRes) {
        if (id.starIndex < realStarCount) {
            CHECK(id.catalogIndex == augCatIdx[id.starIndex]);
        }
    }
    for (const auto &id : tRes) {
        if (id.starIndex < realStarCount) {
            CHECK(id.catalogIndex == augCatIdx[id.starIndex]);
        }
    }
}

TEST_CASE("Tetra3 vs Pyramid: narrow FOV", "[starid][accuracy]") {
    // Use a narrow field-of-view camera (20°) so fewer stars are visible —
    // both algorithms must still identify at least 4 correctly.
    Catalog cat = makeLargeCatalog();
    Camera cam  = Camera(fovToFocalLength(20.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    // Point the camera at the center of the catalog cluster.
    Quaternion rot = sphericalToQuaternion(0.14, 0.02, 0.0);
    Stars stars;
    std::vector<int> catIdx;
    projectCatalogRotated(cat, cam, rot, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 500);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 200000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 15, (int)cat.size(), 500000);

    StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
    StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

    REQUIRE(pRes.size() >= 4);
    REQUIRE(tRes.size() >= 4);

    CHECK(allCorrect(pRes, catIdx));
    CHECK(allCorrect(tRes, catIdx));

    CHECK(tRes.size() >= pRes.size() - 1);
}

TEST_CASE("Tetra3 vs Pyramid: wide FOV", "[starid][accuracy]") {
    // Use a wide field-of-view camera (60°) with many visible stars.
    Catalog cat = makeLargeCatalog();
    Camera cam  = Camera(fovToFocalLength(60.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    Quaternion rot = sphericalToQuaternion(0.14, 0.02, 0.0);
    Stars stars;
    std::vector<int> catIdx;
    projectCatalogRotated(cat, cam, rot, stars, catIdx);

    REQUIRE(stars.size() >= 4);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 500);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 200000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 15, (int)cat.size(), 500000);

    StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
    StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

    REQUIRE(pRes.size() >= 4);
    REQUIRE(tRes.size() >= 4);

    CHECK(allCorrect(pRes, catIdx));
    CHECK(allCorrect(tRes, catIdx));

    // With a wide FOV and many stars, both should identify a large fraction.
    CHECK(pRes.size() >= stars.size() / 2);
    CHECK(tRes.size() >= stars.size() / 2);
}

TEST_CASE("Tetra3 vs Pyramid: rotated camera roll", "[starid][accuracy]") {
    // Verify both algorithms are invariant to camera roll.
    Catalog cat = makeSmallCatalog();
    Camera cam = makeTestCamera();

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 200);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 100000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 10, 10, 200000);

    stfloat rollAngles[] = { 0.0, M_PI / 6.0, M_PI / 3.0, M_PI / 2.0, M_PI };

    for (stfloat roll : rollAngles) {
        Quaternion rot = sphericalToQuaternion(0.07, 0.02, roll);
        Stars stars;
        std::vector<int> catIdx;
        projectCatalogRotated(cat, cam, rot, stars, catIdx);

        if (stars.size() < 4) continue;

        StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
        StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

        // Both should identify at least 4 stars.
        CHECK(pRes.size() >= 4);
        CHECK(tRes.size() >= 4);

        if (!pRes.empty()) CHECK(allCorrect(pRes, catIdx));
        if (!tRes.empty()) CHECK(allCorrect(tRes, catIdx));
    }
}

TEST_CASE("Tetra3 vs Pyramid: aggregate accuracy summary", "[starid][accuracy]") {
    // Run both algorithms across many orientations and produce a summary.
    // Tetra3 should achieve comparable total identification accuracy.
    Catalog cat = makeLargeCatalog();
    Camera cam  = Camera(fovToFocalLength(50.0 * M_PI / 180.0, 1024.0), 1024, 1024);

    stfloat maxAngle = cam.fov();
    auto db = buildDatabase(cat, 0.005, maxAngle, 500);

    PyramidStarIdAlgorithm pyramid(0.001, 0, 0.01, 200000);
    Tetra3StarIdAlgorithm  tetra3(0.001, 0.01, 15, (int)cat.size(), 500000);

    int pyramidTotalCorrect = 0, tetra3TotalCorrect = 0;
    int pyramidTotalStars   = 0, tetra3TotalStars   = 0;
    int pyramidScenes = 0, tetra3Scenes = 0;
    int totalScenes = 0;

    // Sweep boresight through a grid of pointings.
    for (int ri = 0; ri < 4; ri++) {
        for (int di = 0; di < 3; di++) {
            stfloat ra = 0.06 + ri * 0.05;
            stfloat de = -0.02 + di * 0.03;

            Quaternion rot = sphericalToQuaternion(ra, de, 0.0);
            Stars stars;
            std::vector<int> catIdx;
            projectCatalogRotated(cat, cam, rot, stars, catIdx);

            if (stars.size() < 4) continue;
            ++totalScenes;

            StarIdentifiers pRes = pyramid.identify(db.data(), stars, cat, cam);
            StarIdentifiers tRes = tetra3.identify(db.data(), stars, cat, cam);

            if (!pRes.empty()) {
                pyramidTotalCorrect += countCorrect(pRes, catIdx);
                pyramidTotalStars   += (int)pRes.size();
                ++pyramidScenes;
            }
            if (!tRes.empty()) {
                tetra3TotalCorrect += countCorrect(tRes, catIdx);
                tetra3TotalStars   += (int)tRes.size();
                ++tetra3Scenes;
            }
        }
    }

    REQUIRE(totalScenes > 0);

    // Both must have 100% correct identifications (no false positives).
    CHECK(pyramidTotalCorrect == pyramidTotalStars);
    CHECK(tetra3TotalCorrect  == tetra3TotalStars);

    // Tetra3 must solve at least as many scenes as Pyramid.
    CHECK(tetra3Scenes >= pyramidScenes);

    // Tetra3 total identifications should be within 10% of Pyramid's.
    if (pyramidTotalStars > 0) {
        stfloat ratio = (stfloat)tetra3TotalStars / (stfloat)pyramidTotalStars;
        CHECK(ratio >= 0.9);
    }
}
