#include "database.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
//  KVectorIndex
// ═══════════════════════════════════════════════════════════════════════════

const int32_t PairDistanceKVectorDatabase::kMagicValue = 0x2536f009;

KVectorIndex::KVectorIndex(DeserializeContext *des) {
    numValues_ = deserializePrimitive<int32_t>(des);
    min_       = deserializePrimitive<stfloat>(des);
    max_       = deserializePrimitive<stfloat>(des);
    numBins_   = deserializePrimitive<int32_t>(des);

    assert(min_ >= (stfloat)0.0);
    assert(max_ > min_);
    binWidth_ = (max_ - min_) / (stfloat)numBins_;

    bins_ = deserializeArray<int32_t>(des, numBins_ + 1);
}

long KVectorIndex::binFor(stfloat query) const {
    long result = (long)stceil((query - min_) / binWidth_);
    if (result < 0)          result = 0;
    if (result > numBins_)   result = numBins_;
    return result;
}

long KVectorIndex::queryLiberal(stfloat minDist, stfloat maxDist, long *upperIndex) const {
    assert(maxDist > minDist);

    if (maxDist >= max_) maxDist = max_ - (stfloat)0.00001;
    if (minDist <= min_) minDist = min_ + (stfloat)0.00001;
    if (minDist > max_ || maxDist < min_) {
        *upperIndex = 0;
        return 0;
    }

    long lowerBin = binFor(minDist);
    long upperBin = binFor(maxDist);
    assert(upperBin >= lowerBin);
    assert(upperBin <= numBins_);

    long lowerIdx = bins_[lowerBin - 1];
    if (lowerIdx >= numValues_) {
        *upperIndex = 0;
        return 0;
    }
    *upperIndex = bins_[upperBin];
    return lowerIdx;
}

void serializeKVectorIndex(SerializeContext *ser,
                           const std::vector<stfloat> &values,
                           stfloat min, stfloat max, long numBins) {
    stfloat binWidth = (max - min) / (stfloat)numBins;
    std::vector<int32_t> kVector(numBins + 1);

    long lastBin = 0;
    for (int32_t i = 0; i < (int32_t)values.size(); i++) {
        long thisBin = (long)stceil((values[i] - min) / binWidth);
        if (thisBin < 0)        thisBin = 0;
        if (thisBin > numBins)  thisBin = numBins;
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin] = i;
        }
        lastBin = thisBin;
    }
    for (long bin = lastBin; bin <= numBins; bin++) {
        kVector[bin] = (int32_t)values.size();
    }

    serializePrimitive<int32_t>(ser, (int32_t)values.size());
    serializePrimitive<stfloat>(ser, min);
    serializePrimitive<stfloat>(ser, max);
    serializePrimitive<int32_t>(ser, numBins);

    for (const int32_t &bin : kVector) {
        serializePrimitive<int32_t>(ser, bin);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  PairDistanceKVectorDatabase
// ═══════════════════════════════════════════════════════════════════════════

struct KVectorPair {
    int16_t index1;
    int16_t index2;
    stfloat distance;
};

static bool compareKVectorPairs(const KVectorPair &a, const KVectorPair &b) {
    return a.distance < b.distance;
}

static std::vector<KVectorPair> catalogToPairDistances(const Catalog &catalog,
                                                       stfloat minDist,
                                                       stfloat maxDist) {
    std::vector<KVectorPair> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i + 1; k < (int16_t)catalog.size(); k++) {
            stfloat d = angleUnit(catalog[i].spatial, catalog[k].spatial);
            if (d >= minDist && d <= maxDist) {
                result.push_back({i, k, d});
            }
        }
    }
    return result;
}

PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(DeserializeContext *des)
    : index_(des) {
    pairs_ = deserializeArray<int16_t>(des, 2 * index_.numValues());
}

const int16_t *PairDistanceKVectorDatabase::findPairsLiberal(
    stfloat minDist, stfloat maxDist, const int16_t **end) const {

    long upperIdx = -1;
    long lowerIdx = index_.queryLiberal(minDist, maxDist, &upperIdx);
    *end = &pairs_[upperIdx * 2];
    return &pairs_[lowerIdx * 2];
}

void serializePairDistanceKVector(SerializeContext *ser,
                                  const Catalog &catalog,
                                  stfloat minDistance, stfloat maxDistance,
                                  long numBins) {
    auto pairs = catalogToPairDistances(catalog, minDistance, maxDistance);
    std::sort(pairs.begin(), pairs.end(), compareKVectorPairs);

    std::vector<stfloat> distances;
    distances.reserve(pairs.size());
    for (const auto &p : pairs) distances.push_back(p.distance);

    serializeKVectorIndex(ser, distances, minDistance, maxDistance, numBins);

    for (const auto &p : pairs) {
        serializePrimitive<int16_t>(ser, p.index1);
        serializePrimitive<int16_t>(ser, p.index2);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  MultiDatabase
// ═══════════════════════════════════════════════════════════════════════════
//
//  Layout (repeating):
//    int32_t  magicValue
//    uint32_t length          (byte count of the sub-database)
//    <padding to 8-byte boundary>
//    unsigned char[length]    sub-database bytes
//    ...
//    int32_t  0               (caboose)
//

const unsigned char *MultiDatabase::subDatabasePointer(int32_t magicValue) const {
    DeserializeContext des(buffer_);
    while (true) {
        int32_t curMagic = deserializePrimitive<int32_t>(&des);
        if (curMagic == 0) return nullptr;                   // caboose

        uint32_t dbLength = deserializePrimitive<uint32_t>(&des);
        deserializePadding<uint64_t>(&des);                  // align to 8 bytes
        const unsigned char *subPtr = deserializeArray<unsigned char>(&des, dbLength);
        if (curMagic == magicValue) return subPtr;
    }
}

void serializeMultiDatabase(SerializeContext *ser,
                            const MultiDatabaseDescriptor &dbs) {
    for (const auto &entry : dbs) {
        serializePrimitive<int32_t>(ser, entry.magicValue);
        serializePrimitive<uint32_t>(ser, (uint32_t)entry.bytes.size());
        serializePadding<uint64_t>(ser);
        ser->buffer.insert(ser->buffer.end(),
                           entry.bytes.begin(), entry.bytes.end());
    }
    serializePrimitive<int32_t>(ser, 0);   // caboose
}

// ═══════════════════════════════════════════════════════════════════════════
//  Catalog serialization
// ═══════════════════════════════════════════════════════════════════════════

void serializeCatalog(SerializeContext *ser, const Catalog &catalog) {
    serializePrimitive<int16_t>(ser, (int16_t)catalog.size());
    for (const CatalogStar &s : catalog) {
        serializePrimitive<stfloat>(ser, s.spatial.x);
        serializePrimitive<stfloat>(ser, s.spatial.y);
        serializePrimitive<stfloat>(ser, s.spatial.z);
        serializePrimitive<int32_t>(ser, s.magnitude);
        serializePrimitive<int16_t>(ser, s.name);
    }
}

Catalog deserializeCatalog(DeserializeContext *des) {
    int16_t n = deserializePrimitive<int16_t>(des);
    Catalog result;
    result.reserve(n);
    for (int16_t i = 0; i < n; i++) {
        CatalogStar s;
        s.spatial.x = deserializePrimitive<stfloat>(des);
        s.spatial.y = deserializePrimitive<stfloat>(des);
        s.spatial.z = deserializePrimitive<stfloat>(des);
        s.magnitude = deserializePrimitive<int32_t>(des);
        s.name      = deserializePrimitive<int16_t>(des);
        result.push_back(s);
    }
    return result;
}
