#pragma once

#include "star.hpp"

#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

// ═══════════════════════════════════════════════════════════════════════════
//  Serialization / Deserialization helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Accumulates bytes during database construction.
struct SerializeContext {
    std::vector<unsigned char> buffer;
};

/// Reads from a flat byte buffer produced by SerializeContext.
struct DeserializeContext {
    explicit DeserializeContext(const unsigned char *buf)
        : buffer(buf), cursor(buf) {}

    size_t offset() const { return (size_t)(cursor - buffer); }

    void advance(size_t n) { cursor += n; }

    const unsigned char *getCursor() const { return cursor; }

private:
    const unsigned char *buffer;
    const unsigned char *cursor;
};

/// Align the cursor forward so it sits on a sizeof(T) boundary.
template <typename T>
void deserializePadding(DeserializeContext *des) {
    size_t rem = des->offset() % sizeof(T);
    if (rem != 0) des->advance(sizeof(T) - rem);
}

/// Read one primitive from the buffer.
template <typename T>
T deserializePrimitive(DeserializeContext *des) {
    deserializePadding<T>(des);
    T val;
    std::memcpy(&val, des->getCursor(), sizeof(T));
    des->advance(sizeof(T));
    return val;
}

/// Return a pointer into the existing buffer (zero-copy).
template <typename T>
const T *deserializeArray(DeserializeContext *des, long count) {
    deserializePadding<T>(des);
    const T *p = reinterpret_cast<const T *>(des->getCursor());
    des->advance(sizeof(T) * count);
    return p;
}

/// Pad the output buffer to a sizeof(T) boundary.
template <typename T>
void serializePadding(SerializeContext *ser) {
    while (ser->buffer.size() % sizeof(T) != 0)
        ser->buffer.push_back(0);
}

/// Append one primitive value to the output buffer.
template <typename T>
void serializePrimitive(SerializeContext *ser, const T &val) {
    serializePadding<T>(ser);
    const unsigned char *p = reinterpret_cast<const unsigned char *>(&val);
    ser->buffer.insert(ser->buffer.end(), p, p + sizeof(T));
}

// ═══════════════════════════════════════════════════════════════════════════
//  K-Vector index  (enables O(1) range queries on sorted numerical data)
// ═══════════════════════════════════════════════════════════════════════════

class KVectorIndex {
public:
    /// Deserialize from an existing buffer.
    explicit KVectorIndex(DeserializeContext *des);

    /// Return the lower index into the sorted values array such that all
    /// entries in [lower, *upperIndex) *may* fall in [minDist, maxDist].
    long queryLiberal(stfloat minDist, stfloat maxDist, long *upperIndex) const;

    long    numValues() const { return numValues_; }
    long    numBins()   const { return numBins_; }
    stfloat maxVal()    const { return max_; }
    stfloat minVal()    const { return min_; }

private:
    long binFor(stfloat dist) const;

    long     numValues_;
    stfloat  min_;
    stfloat  max_;
    stfloat  binWidth_;
    long     numBins_;
    const int32_t *bins_;          // points into the deserialized buffer
};

/// Serialize a KVector index for a pre-sorted array of values.
void serializeKVectorIndex(SerializeContext *ser,
                           const std::vector<stfloat> &sortedValues,
                           stfloat min, stfloat max, long numBins);

// ═══════════════════════════════════════════════════════════════════════════
//  Pair-distance K-Vector database
// ═══════════════════════════════════════════════════════════════════════════

class PairDistanceKVectorDatabase {
public:
    /// Deserialize from buffer.
    explicit PairDistanceKVectorDatabase(DeserializeContext *des);

    /// Return pointer-range of (index1, index2) pairs whose angular distance
    /// is approximately within [min, max].  Each pair is two consecutive int16_t.
    const int16_t *findPairsLiberal(stfloat min, stfloat max,
                                    const int16_t **end) const;

    stfloat maxDistance() const { return index_.maxVal(); }
    stfloat minDistance() const { return index_.minVal(); }
    long    numPairs()    const { return index_.numValues(); }

    static const int32_t kMagicValue;

private:
    KVectorIndex    index_;
    const int16_t  *pairs_;   // points into the deserialized buffer
};

/// Build and serialize a pair-distance KVector from a catalog.
void serializePairDistanceKVector(SerializeContext *ser,
                                  const Catalog &catalog,
                                  stfloat minDistance, stfloat maxDistance,
                                  long numBins);

// ═══════════════════════════════════════════════════════════════════════════
//  Multi-database container  (maps magic values → sub-database buffers)
// ═══════════════════════════════════════════════════════════════════════════

/// One entry destined for a MultiDatabase.
struct MultiDatabaseEntry {
    int32_t                    magicValue;
    std::vector<unsigned char> bytes;

    MultiDatabaseEntry(int32_t mv, std::vector<unsigned char> b)
        : magicValue(mv), bytes(std::move(b)) {}
};

using MultiDatabaseDescriptor = std::vector<MultiDatabaseEntry>;

/// Container that holds one or more sub-databases keyed by magic value.
class MultiDatabase {
public:
    explicit MultiDatabase(const unsigned char *buf) : buffer_(buf) {}

    /// Returns a pointer to the start of the sub-database identified by
    /// magicValue, or nullptr if not found.
    const unsigned char *subDatabasePointer(int32_t magicValue) const;

private:
    const unsigned char *buffer_;
};

/// Serialize a set of sub-databases into one flat buffer.
void serializeMultiDatabase(SerializeContext *ser,
                            const MultiDatabaseDescriptor &dbs);

// ═══════════════════════════════════════════════════════════════════════════
//  Catalog serialization helpers
// ═══════════════════════════════════════════════════════════════════════════

static const int32_t kCatalogMagicValue = 0xF9A283BC;

void serializeCatalog(SerializeContext *ser, const Catalog &catalog);
Catalog deserializeCatalog(DeserializeContext *des);
