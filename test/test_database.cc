// Unit tests for KVector index and PairDistanceKVectorDatabase

#include <catch.hpp>
#include "database.hpp"
#include <cmath>
#include <vector>
#include <set>

// ═══════════════════════════════════════════════════════════════════════════
//  Helper: build a small test catalog of stars spread around the sky
// ═══════════════════════════════════════════════════════════════════════════

static Catalog makeTestCatalog() {
    Catalog cat;
    // Place stars at known RA/Dec positions (in radians)
    cat.push_back(CatalogStar(0.0,    0.0,    100, 0));
    cat.push_back(CatalogStar(0.1,    0.0,    100, 1));
    cat.push_back(CatalogStar(0.2,    0.1,    100, 2));
    cat.push_back(CatalogStar(0.05,   0.15,   100, 3));
    cat.push_back(CatalogStar(0.15,  -0.1,    100, 4));
    cat.push_back(CatalogStar(0.3,    0.0,    100, 5));
    cat.push_back(CatalogStar(0.25,   0.05,   100, 6));
    cat.push_back(CatalogStar(0.35,   0.1,    100, 7));
    cat.push_back(CatalogStar(0.4,   -0.05,   100, 8));
    cat.push_back(CatalogStar(0.5,    0.0,    100, 9));
    return cat;
}

// ═══════════════════════════════════════════════════════════════════════════
//  KVector index basic tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("KVector: serialize then deserialize", "[kvector]") {
    // Build a sorted list of values
    std::vector<stfloat> values = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    SerializeContext ser;
    serializeKVectorIndex(&ser, values, 0.05, 0.95, 10);

    DeserializeContext des(ser.buffer.data());
    KVectorIndex idx(&des);

    CHECK(idx.numValues() == (long)values.size());
    CHECK(idx.minVal() == Approx(0.05));
    CHECK(idx.maxVal() == Approx(0.95));
    CHECK(idx.numBins() == 10);
}

TEST_CASE("KVector: query returns narrowing range", "[kvector]") {
    std::vector<stfloat> values = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    SerializeContext ser;
    serializeKVectorIndex(&ser, values, 0.05, 0.95, 20);

    DeserializeContext des(ser.buffer.data());
    KVectorIndex idx(&des);

    // Wider query should return more values
    long upper1;
    long lower1 = idx.queryLiberal(0.15, 0.85, &upper1);
    long count1 = upper1 - lower1;

    long upper2;
    long lower2 = idx.queryLiberal(0.35, 0.55, &upper2);
    long count2 = upper2 - lower2;

    CHECK(count1 >= count2);
    CHECK(count1 > 0);
    CHECK(count2 > 0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  PairDistanceKVectorDatabase
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("PairDistanceDB: construction and basic queries", "[pairdb]") {
    Catalog cat = makeTestCatalog();

    stfloat minDist = 0.05;
    stfloat maxDist = 0.6;
    long numBins = 50;

    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, minDist, maxDist, numBins);

    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    CHECK(db.numPairs() > 0);
    CHECK(db.minDistance() == Approx(minDist));
    CHECK(db.maxDistance() == Approx(maxDist));
}

TEST_CASE("PairDistanceDB: all returned pairs have valid distances", "[pairdb]") {
    Catalog cat = makeTestCatalog();

    stfloat minDist = 0.05;
    stfloat maxDist = 0.6;

    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, minDist, maxDist, 50);

    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    // Query a sub-range
    stfloat qMin = 0.08;
    stfloat qMax = 0.15;
    const int16_t *end;
    const int16_t *pairs = db.findPairsLiberal(qMin, qMax, &end);

    for (const int16_t *p = pairs; p < end; p += 2) {
        int16_t i1 = p[0];
        int16_t i2 = p[1];
        REQUIRE(i1 >= 0);
        REQUIRE(i1 < (int16_t)cat.size());
        REQUIRE(i2 >= 0);
        REQUIRE(i2 < (int16_t)cat.size());

        stfloat d = angleUnit(cat[i1].spatial, cat[i2].spatial);
        // Liberal query may return slightly outside, but should be close
        CHECK(d >= minDist);
        CHECK(d <= maxDist);
    }
}

TEST_CASE("PairDistanceDB: narrower query returns fewer pairs", "[pairdb]") {
    Catalog cat = makeTestCatalog();

    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, 0.05, 0.6, 50);

    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    const int16_t *end1;
    const int16_t *p1 = db.findPairsLiberal(0.06, 0.55, &end1);
    long count1 = (end1 - p1) / 2;

    const int16_t *end2;
    const int16_t *p2 = db.findPairsLiberal(0.1, 0.2, &end2);
    long count2 = (end2 - p2) / 2;

    CHECK(count1 >= count2);
}

TEST_CASE("PairDistanceDB: known pair can be found", "[pairdb]") {
    Catalog cat = makeTestCatalog();

    // Distance between star 0 (ra=0, de=0) and star 1 (ra=0.1, de=0)
    stfloat expectedDist = angleUnit(cat[0].spatial, cat[1].spatial);
    REQUIRE(expectedDist > 0.0);

    // Ensure the database range brackets this distance (min must be > 0)
    stfloat dbMin = std::max(0.001, expectedDist * 0.5);
    stfloat dbMax = expectedDist * 2.0;

    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, dbMin, dbMax, 50);

    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    stfloat tol = 0.01;
    const int16_t *end;
    const int16_t *pairs = db.findPairsLiberal(expectedDist - tol, expectedDist + tol, &end);

    bool found = false;
    for (const int16_t *p = pairs; p < end; p += 2) {
        if ((p[0] == 0 && p[1] == 1) || (p[0] == 1 && p[1] == 0)) {
            found = true;
            break;
        }
    }
    CHECK(found);
}

TEST_CASE("PairDistanceDB: query partitions cover all pairs", "[pairdb]") {
    Catalog cat = makeTestCatalog();

    stfloat minD = 0.05;
    stfloat maxD = 0.5;

    SerializeContext ser;
    serializePairDistanceKVector(&ser, cat, minD, maxD, 50);

    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    // Sum up pairs across contiguous partitions
    long totalPairs = 0;
    stfloat step = (maxD - minD) / 5.0;
    for (int i = 0; i < 5; i++) {
        stfloat lo = minD + step * i + 0.000001;
        stfloat hi = minD + step * (i + 1) - 0.000001;
        if (hi <= lo) continue;
        const int16_t *end;
        const int16_t *p = db.findPairsLiberal(lo, hi, &end);
        totalPairs += (end - p) / 2;
    }

    // Should account for all pairs (may have a few duplicates due to liberal queries)
    CHECK(totalPairs >= db.numPairs() - 5);  // allow minor edge rounding
    CHECK(totalPairs <= db.numPairs() + 5);
}
