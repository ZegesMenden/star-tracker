// Unit tests for serialization/deserialization helpers and catalog round-trip

#include <catch.hpp>
#include "database.hpp"
#include <cmath>
#include <cstdint>

// ═══════════════════════════════════════════════════════════════════════════
//  Primitive serialization
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Serialize/deserialize int32_t", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int32_t>(&ser, 42);
    DeserializeContext des(ser.buffer.data());
    int32_t val = deserializePrimitive<int32_t>(&des);
    CHECK(val == 42);
}

TEST_CASE("Serialize/deserialize stfloat", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<stfloat>(&ser, 3.14159);
    DeserializeContext des(ser.buffer.data());
    stfloat val = deserializePrimitive<stfloat>(&des);
    CHECK(val == Approx(3.14159));
}

TEST_CASE("Serialize/deserialize int16_t", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int16_t>(&ser, -1234);
    DeserializeContext des(ser.buffer.data());
    int16_t val = deserializePrimitive<int16_t>(&des);
    CHECK(val == -1234);
}

TEST_CASE("Serialize/deserialize int64_t", "[serialize]") {
    SerializeContext ser;
    int64_t bigVal = 27837492938LL;
    serializePrimitive<int64_t>(&ser, bigVal);
    DeserializeContext des(ser.buffer.data());
    int64_t val = deserializePrimitive<int64_t>(&des);
    CHECK(val == bigVal);
}

TEST_CASE("Multiple primitives round-trip", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int64_t>(&ser, 99999LL);
    serializePrimitive<stfloat>(&ser, 2.71828);

    DeserializeContext des(ser.buffer.data());
    CHECK(deserializePrimitive<int64_t>(&des) == 99999LL);
    CHECK(deserializePrimitive<stfloat>(&des) == Approx(2.71828));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Padding
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Padding: int8_t then int32_t aligns to 4 bytes", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int8_t>(&ser, 23);
    serializePrimitive<int32_t>(&ser, 1234567);
    // int8_t is 1 byte + 3 padding + 4 bytes = 8
    CHECK(ser.buffer.size() == 8);

    DeserializeContext des(ser.buffer.data());
    CHECK(deserializePrimitive<int8_t>(&des) == 23);
    CHECK(deserializePrimitive<int32_t>(&des) == 1234567);
}

TEST_CASE("Padding: no padding needed when already aligned", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int32_t>(&ser, 111);
    serializePrimitive<int32_t>(&ser, 222);
    CHECK(ser.buffer.size() == 8);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Array deserialization (zero-copy)
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("deserializeArray reads back correct values", "[serialize]") {
    SerializeContext ser;
    serializePrimitive<int32_t>(&ser, 10);
    serializePrimitive<int32_t>(&ser, 20);
    serializePrimitive<int32_t>(&ser, 30);

    DeserializeContext des(ser.buffer.data());
    const int32_t *arr = deserializeArray<int32_t>(&des, 3);
    CHECK(arr[0] == 10);
    CHECK(arr[1] == 20);
    CHECK(arr[2] == 30);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Catalog serialization round-trip
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("Catalog serialize -> deserialize round-trip", "[serialize][catalog]") {
    Catalog cat;
    cat.push_back(CatalogStar(0.0, 0.0, 100, 1));
    cat.push_back(CatalogStar(1.0, 0.5, 200, 2));
    cat.push_back(CatalogStar(M_PI, -0.3, 350, 3));

    SerializeContext ser;
    serializeCatalog(&ser, cat);

    DeserializeContext des(ser.buffer.data());
    Catalog result = deserializeCatalog(&des);

    REQUIRE(result.size() == cat.size());
    for (size_t i = 0; i < cat.size(); i++) {
        CHECK(result[i].spatial.x == Approx(cat[i].spatial.x).margin(1e-6));
        CHECK(result[i].spatial.y == Approx(cat[i].spatial.y).margin(1e-6));
        CHECK(result[i].spatial.z == Approx(cat[i].spatial.z).margin(1e-6));
        CHECK(result[i].magnitude == cat[i].magnitude);
        CHECK(result[i].name == cat[i].name);
    }
}

TEST_CASE("Empty catalog round-trip", "[serialize][catalog]") {
    Catalog cat;
    SerializeContext ser;
    serializeCatalog(&ser, cat);

    DeserializeContext des(ser.buffer.data());
    Catalog result = deserializeCatalog(&des);
    CHECK(result.empty());
}

// ═══════════════════════════════════════════════════════════════════════════
//  MultiDatabase
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("MultiDatabase: store and retrieve sub-databases", "[serialize][multidb]") {
    // Build two sub-databases
    SerializeContext sub1;
    serializePrimitive<int32_t>(&sub1, 0xDEAD);

    SerializeContext sub2;
    serializePrimitive<int32_t>(&sub2, 0xBEEF);

    MultiDatabaseDescriptor desc;
    desc.push_back(MultiDatabaseEntry(0x1111, sub1.buffer));
    desc.push_back(MultiDatabaseEntry(0x2222, sub2.buffer));

    SerializeContext multiSer;
    serializeMultiDatabase(&multiSer, desc);

    MultiDatabase mdb(multiSer.buffer.data());

    // Retrieve first sub-database
    const unsigned char *p1 = mdb.subDatabasePointer(0x1111);
    REQUIRE(p1 != nullptr);
    DeserializeContext d1(p1);
    CHECK(deserializePrimitive<int32_t>(&d1) == 0xDEAD);

    // Retrieve second sub-database
    const unsigned char *p2 = mdb.subDatabasePointer(0x2222);
    REQUIRE(p2 != nullptr);
    DeserializeContext d2(p2);
    CHECK(deserializePrimitive<int32_t>(&d2) == 0xBEEF);
}

TEST_CASE("MultiDatabase: missing magic value returns nullptr", "[serialize][multidb]") {
    MultiDatabaseDescriptor desc;
    SerializeContext sub;
    serializePrimitive<int32_t>(&sub, 1);
    desc.push_back(MultiDatabaseEntry(0xAAAA, sub.buffer));

    SerializeContext multiSer;
    serializeMultiDatabase(&multiSer, desc);

    MultiDatabase mdb(multiSer.buffer.data());
    CHECK(mdb.subDatabasePointer(0xBBBB) == nullptr);
}
