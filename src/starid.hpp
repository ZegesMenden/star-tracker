#pragma once

#include "star.hpp"
#include "camera.hpp"
#include "database.hpp"

#include <vector>
#include <unordered_map>
#include <cstdint>

// ═══════════════════════════════════════════════════════════════════════════
//  Star-identification algorithm base class
// ═══════════════════════════════════════════════════════════════════════════

/// Abstract base – different algorithms inherit and override `identify()`.
class StarIdAlgorithm {
public:
    virtual ~StarIdAlgorithm() = default;

    /// Given a serialized database, detected centroids, catalog, and camera
    /// model, return a list of (centroid → catalog) identifications.
    virtual StarIdentifiers identify(const unsigned char *database,
                                     const Stars          &stars,
                                     const Catalog        &catalog,
                                     const Camera         &camera) const = 0;
};

// ═══════════════════════════════════════════════════════════════════════════
//  Pyramid star-id
// ═══════════════════════════════════════════════════════════════════════════

/// The Pyramid algorithm searches 4-star patterns (pyramids) in the image,
/// matches them against a pair-distance KVector database, and then identifies
/// all remaining stars by triangulation.
///
/// Reference: Mortari, Samaan, Bruccoleri, Junkins.
/// "The Pyramid Star Identification Technique", Navigation 51(3), 2004.
class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    /// @param tolerance   Angular tolerance in radians for distance matching.
    /// @param numFalseStars  Estimated false-star count over the full sphere.
    /// @param maxMismatchProb  Maximum acceptable probability of a false match.
    /// @param cutoff      Maximum number of pyramid iterations before giving up.
    PyramidStarIdAlgorithm(stfloat tolerance, int numFalseStars,
                           stfloat maxMismatchProb, long cutoff)
        : tolerance_(tolerance),
          numFalseStars_(numFalseStars),
          maxMismatchProb_(maxMismatchProb),
          cutoff_(cutoff) {}

    StarIdentifiers identify(const unsigned char *database,
                             const Stars          &stars,
                             const Catalog        &catalog,
                             const Camera         &camera) const override;

private:
    stfloat tolerance_;
    int     numFalseStars_;
    stfloat maxMismatchProb_;
    long    cutoff_;
};

// ═══════════════════════════════════════════════════════════════════════════
//  Tetra3 star-id
// ═══════════════════════════════════════════════════════════════════════════

/// Tetra3-style geometric hashing over 4-star patterns.
///
/// For each 4-star pattern, we compute six angular distances, divide the five
/// smaller edges by the largest edge, and quantize those ratios into a hash
/// key. Image patterns query neighboring hash buckets and validate candidates
/// by checking all six angular distances.
///
/// Reference implementation: https://github.com/esa/tetra3
class Tetra3StarIdAlgorithm : public StarIdAlgorithm {
public:
    /// @param tolerance Angular tolerance in radians for final distance checks.
    /// @param patternMaxError Quantization step for normalized edge-ratio hash.
    /// @param patternCheckingStars Number of image stars to use for tetra search.
    /// @param catalogPatternStars Number of brightest catalog stars to hash.
    /// @param patternCutoff Maximum catalog tetra patterns to generate.
    Tetra3StarIdAlgorithm(stfloat tolerance,
                          stfloat patternMaxError = (stfloat)0.01,
                          int patternCheckingStars = 10,
                          int catalogPatternStars = 48,
                          long patternCutoff = 200000)
        : tolerance_(tolerance),
          patternMaxError_(patternMaxError),
          patternCheckingStars_(patternCheckingStars),
          catalogPatternStars_(catalogPatternStars),
          patternCutoff_(patternCutoff) {}

    StarIdentifiers identify(const unsigned char *database,
                             const Stars          &stars,
                             const Catalog        &catalog,
                             const Camera         &camera) const override;

private:
    stfloat tolerance_;
    stfloat patternMaxError_;
    int     patternCheckingStars_;
    int     catalogPatternStars_;
    long    patternCutoff_;
};

// ═══════════════════════════════════════════════════════════════════════════
//  Internal helpers  (exposed for unit testing)
// ═══════════════════════════════════════════════════════════════════════════

/// Build a multi-map  catalogStar → otherCatalogStar  from a raw pairs range.
std::unordered_multimap<int16_t, int16_t>
pairDistanceQueryToMap(const int16_t *pairs, const int16_t *end);

/// After an initial 4-star match, attempt to identify remaining centroids
/// using triangulated pair-distance queries.
/// Returns the number of *additional* stars identified.
int identifyRemainingStarsPairDistance(StarIdentifiers                   &identified,
                                      const Stars                       &stars,
                                      const PairDistanceKVectorDatabase &db,
                                      const Catalog                     &catalog,
                                      const Camera                      &camera,
                                      stfloat                            tolerance);
