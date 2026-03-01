#include "starid.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
//  Utility helpers
// ═══════════════════════════════════════════════════════════════════════════

std::unordered_multimap<int16_t, int16_t>
pairDistanceQueryToMap(const int16_t *pairs, const int16_t *end) {
    std::unordered_multimap<int16_t, int16_t> map;
    for (const int16_t *p = pairs; p < end; p += 2) {
        map.emplace(p[0], p[1]);
        map.emplace(p[1], p[0]);
    }
    return map;
}

// ═══════════════════════════════════════════════════════════════════════════
//  IdentifyRemainingStarsPairDistance
// ═══════════════════════════════════════════════════════════════════════════
//
//  After the pyramid gives us 4 identified stars, this routine walks through
//  every unidentified centroid.  For each one it picks the two nearest
//  already-identified centroids, measures the angular distances, queries the
//  pair-distance database twice (once per identified star), and intersects the
//  results — ideally yielding exactly one catalog candidate.
//
//  This is a simplified but functional version of the LOST
//  IdentifyRemainingStarsPairDistance routine.
// ───────────────────────────────────────────────────────────────────────────

/// Given two identified catalog indices and the measured angular distances
/// from each to an unknown centroid, return all catalog stars that are
/// consistent with both distances (within tolerance) AND whose spectrality
/// (chirality) matches the observed configuration.
static std::vector<int16_t>
identifyThirdStar(const PairDistanceKVectorDatabase &db,
                  const Catalog &catalog,
                  int16_t catIdx1, int16_t catIdx2,
                  stfloat dist1, stfloat dist2,
                  stfloat tolerance) {

    const int16_t *end1;
    const int16_t *q1 = db.findPairsLiberal(dist1 - tolerance,
                                             dist1 + tolerance, &end1);
    const int16_t *end2;
    const int16_t *q2 = db.findPairsLiberal(dist2 - tolerance,
                                             dist2 + tolerance, &end2);

    // Collect all catalog stars paired with catIdx1 at distance ~dist1.
    std::unordered_set<int16_t> map1;
    map1.reserve((size_t)((end1 - q1) / 2));
    for (const int16_t *p = q1; p < end1; p += 2) {
        if (p[0] == catIdx1)      map1.insert(p[1]);
        else if (p[1] == catIdx1) map1.insert(p[0]);
    }

    // From those, keep only the ones also appearing in catIdx2's distance query.
    std::vector<int16_t> candidates;
    for (const int16_t *p = q2; p < end2; p += 2) {
        int16_t other = -1;
        if (p[0] == catIdx2)      other = p[1];
        else if (p[1] == catIdx2) other = p[0];
        if (other < 0) continue;

        // Check that `other` also appeared in the first query.
        if (map1.find(other) != map1.end()) {
            // Spectrality (chirality) check:
            // (cat1 × cat2) · candidate  must be positive.
            const Vec3 &s1 = catalog[catIdx1].spatial;
            const Vec3 &s2 = catalog[catIdx2].spatial;
            const Vec3 &sc = catalog[other].spatial;
            if (s1.cross(s2) * sc > 0) {
                candidates.push_back(other);
            }
        }
    }
    return candidates;
}

int identifyRemainingStarsPairDistance(StarIdentifiers                   &identified,
                                      const Stars                       &stars,
                                      const PairDistanceKVectorDatabase &db,
                                      const Catalog                     &catalog,
                                      const Camera                      &camera,
                                      stfloat                            tolerance) {
    int numExtra = 0;

    // We may loop several passes: each newly identified star can help
    // identify more stars in the next pass.
    bool progress = true;
    while (progress) {
        progress = false;

        for (int si = 0; si < (int)stars.size(); si++) {
            // Skip if already identified.
            bool alreadyDone = false;
            for (const auto &id : identified) {
                if (id.starIndex == si) { alreadyDone = true; break; }
            }
            if (alreadyDone) continue;

            Vec3 unidSpatial = camera.cameraToSpatial(stars[si].position).normalize();

            // Find the two closest already-identified centroids (by angular
            // distance in the image) that are within the database bounds.
            struct Neighbor { int starIdx; int catIdx; stfloat dist; };
            std::vector<Neighbor> neighbors;
            for (const auto &id : identified) {
                Vec3 idSpatial = camera.cameraToSpatial(stars[id.starIndex].position).normalize();
                stfloat d = angleUnit(unidSpatial, idSpatial);
                if (d >= db.minDistance() && d <= db.maxDistance()) {
                    neighbors.push_back({id.starIndex, id.catalogIndex, d});
                }
            }
            if (neighbors.size() < 2) continue;

            // Sort by distance and pick the two closest
            std::sort(neighbors.begin(), neighbors.end(),
                      [](const Neighbor &a, const Neighbor &b){ return a.dist < b.dist; });

            const Neighbor &n1 = neighbors[0];
            const Neighbor &n2 = neighbors[1];

            // Determine spectrality of the observed triple.
            Vec3 s1 = camera.cameraToSpatial(stars[n1.starIdx].position).normalize();
            Vec3 s2 = camera.cameraToSpatial(stars[n2.starIdx].position).normalize();
            stfloat spectral = s1.cross(s2) * unidSpatial;

            std::vector<int16_t> candidates =
                spectral > 0
                ? identifyThirdStar(db, catalog,
                                    (int16_t)n1.catIdx, (int16_t)n2.catIdx,
                                    n1.dist, n2.dist, tolerance)
                : identifyThirdStar(db, catalog,
                                    (int16_t)n2.catIdx, (int16_t)n1.catIdx,
                                    n2.dist, n1.dist, tolerance);

            if (candidates.size() == 1) {
                identified.push_back(StarIdentifier(si, candidates[0]));
                ++numExtra;
                progress = true;  // new identification may unlock others
            }
        }
    }

    return numExtra;
}

// ═══════════════════════════════════════════════════════════════════════════
//  PyramidStarIdAlgorithm
// ═══════════════════════════════════════════════════════════════════════════

StarIdentifiers PyramidStarIdAlgorithm::identify(
    const unsigned char *database,
    const Stars          &stars,
    const Catalog        &catalog,
    const Camera         &camera) const {

    StarIdentifiers identified;

    MultiDatabase multiDb(database);
    const unsigned char *dbBuf =
        multiDb.subDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (!dbBuf || stars.size() < 4) {
        std::cerr << "Pyramid: not enough stars or database missing." << std::endl;
        return identified;
    }

    DeserializeContext des(dbBuf);
    PairDistanceKVectorDatabase vectorDb(&des);

    // Constant used to estimate the probability that a random false-positive
    // pyramid match would occur.  See the Pyramid paper / HSL wiki
    // "Analytic_Star_Pattern_Probability."
    stfloat expectedMismatchConst =
        stpow((stfloat)numFalseStars_, 4)
        * stpow(tolerance_, 5) / (stfloat)2.0
        / (stfloat)(M_PI * M_PI);

    int numStars = (int)stars.size();
    int across       = (int)(stsqrt((stfloat)numStars)) * 2;
    int halfAcross   = (int)(stsqrt((stfloat)numStars) / (stfloat)2.0);
    long totalIter   = 0;

    int jMax = numStars - 3;
    for (int jIter = 0; jIter < jMax; jIter++) {
        int dj = 1 + (jIter + halfAcross) % jMax;

        int kMax = numStars - dj - 2;
        for (int kIter = 0; kIter < kMax; kIter++) {
            int dk = 1 + (kIter + across) % kMax;

            int rMax = numStars - dj - dk - 1;
            for (int rIter = 0; rIter < rMax; rIter++) {
                int dr = 1 + (rIter + halfAcross) % rMax;

                int iMax = numStars - dj - dk - dr - 1;
                for (int iIter = 0; iIter <= iMax; iIter++) {
                    int i = (iIter + iMax / 2) % (iMax + 1);

                    if (++totalIter > cutoff_) {
                        std::cerr << "Pyramid: cutoff reached." << std::endl;
                        return identified;
                    }

                    int j = i + dj;
                    int k = j + dk;
                    int r = k + dr;

                    // Project the four centroid positions to 3-D unit vectors.
                    Vec3 iSp = camera.cameraToSpatial(stars[i].position).normalize();
                    Vec3 jSp = camera.cameraToSpatial(stars[j].position).normalize();
                    Vec3 kSp = camera.cameraToSpatial(stars[k].position).normalize();

                    stfloat ijDist = angleUnit(iSp, jSp);

                    // Triangle inner angles (used for mismatch probability).
                    stfloat iSinInner = stsin(angle(jSp - iSp, kSp - iSp));
                    stfloat jSinInner = stsin(angle(iSp - jSp, kSp - jSp));
                    stfloat kSinInner = stsin(angle(iSp - kSp, jSp - kSp));

                    stfloat expectedMismatches =
                        expectedMismatchConst
                        * stsin(ijDist)
                        / kSinInner
                        / std::max({iSinInner, jSinInner, kSinInner});

                    if (expectedMismatches > maxMismatchProb_) continue;

                    Vec3 rSp = camera.cameraToSpatial(stars[r].position).normalize();

                    // Spectrality: sign of the determinant  [i, j, k].
                    bool spectral = (iSp.cross(jSp) * kSp) > 0;

                    stfloat ikDist = angleUnit(iSp, kSp);
                    stfloat irDist = angleUnit(iSp, rSp);
                    stfloat jkDist = angleUnit(jSp, kSp);
                    stfloat jrDist = angleUnit(jSp, rSp);
                    stfloat krDist = angleUnit(kSp, rSp);

                    // All six distances must lie within the database range
                    // (with tolerance margin) to avoid edge-case mismatches.
                    stfloat lo = vectorDb.minDistance() + tolerance_;
                    stfloat hi = vectorDb.maxDistance() - tolerance_;
                    auto inRange = [lo, hi](stfloat d) {
                        return d >= lo && d <= hi;
                    };
                    if (!inRange(ikDist) || !inRange(irDist) ||
                        !inRange(jkDist) || !inRange(jrDist) ||
                        !inRange(krDist))
                        continue;

                    // ── Query the database ───────────────────────────────
                    const int16_t *ijEnd, *ikEnd, *irEnd;
                    const int16_t *ijQ = vectorDb.findPairsLiberal(
                        ijDist - tolerance_, ijDist + tolerance_, &ijEnd);
                    const int16_t *ikQ = vectorDb.findPairsLiberal(
                        ikDist - tolerance_, ikDist + tolerance_, &ikEnd);
                    const int16_t *irQ = vectorDb.findPairsLiberal(
                        irDist - tolerance_, irDist + tolerance_, &irEnd);

                    auto ikMap = pairDistanceQueryToMap(ikQ, ikEnd);
                    auto irMap = pairDistanceQueryToMap(irQ, irEnd);

                    // ── Search for a unique 4-star match ─────────────────
                    int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
                    bool duplicate = false;

                    for (const int16_t *p = ijQ; p < ijEnd; p += 1) {
                        int iCand = *p;
                        // Pairs are stored as consecutive (a, b).
                        // Depending on parity, one is the "i" candidate and
                        // the other is the "j" candidate.
                        int jCand = ((p - ijQ) % 2 == 0) ? p[1] : p[-1];

                        const Vec3 &iCandSp = catalog[iCand].spatial;
                        const Vec3 &jCandSp = catalog[jCand].spatial;
                        Vec3 ijCross = iCandSp.cross(jCandSp);

                        // ── k candidates (via ik map) ────────────────────
                        auto kRange = ikMap.equal_range((int16_t)iCand);
                        for (auto kIt = kRange.first; kIt != kRange.second; ++kIt) {
                            int kCand = kIt->second;
                            const Vec3 &kCandSp = catalog[kCand].spatial;

                            // Spectrality must match.
                            if ((ijCross * kCandSp > 0) != spectral) continue;

                            // Verify jk distance.
                            stfloat jkCandDist = angleUnit(jCandSp, kCandSp);
                            if (stfabs(jkCandDist - jkDist) > tolerance_) continue;

                            // ── r candidates (via ir map) ────────────────
                            auto rRange = irMap.equal_range((int16_t)iCand);
                            for (auto rIt = rRange.first; rIt != rRange.second; ++rIt) {
                                int rCand = rIt->second;
                                const Vec3 &rCandSp = catalog[rCand].spatial;

                                stfloat jrCandDist = angleUnit(jCandSp, rCandSp);
                                if (stfabs(jrCandDist - jrDist) > tolerance_) continue;

                                stfloat krCandDist = angleUnit(kCandSp, rCandSp);
                                if (stfabs(krCandDist - krDist) > tolerance_) continue;

                                // All 6 distances verified!
                                if (iMatch == -1) {
                                    iMatch = iCand;
                                    jMatch = jCand;
                                    kMatch = kCand;
                                    rMatch = rCand;
                                } else {
                                    // Non-unique — skip this pyramid.
                                    duplicate = true;
                                    goto nextPyramid;
                                }
                            }
                        }
                    }

                    if (iMatch != -1 && !duplicate) {
                        identified.push_back(StarIdentifier(i, iMatch));
                        identified.push_back(StarIdentifier(j, jMatch));
                        identified.push_back(StarIdentifier(k, kMatch));
                        identified.push_back(StarIdentifier(r, rMatch));

                        identifyRemainingStarsPairDistance(
                            identified, stars, vectorDb, catalog,
                            camera, tolerance_);

                        return identified;
                    }

                nextPyramid:;
                }
            }
        }
    }

    std::cerr << "Pyramid: no unique match found." << std::endl;
    return identified;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tetra3StarIdAlgorithm internals
// ═══════════════════════════════════════════════════════════════════════════

namespace {

struct RatioHashKey {
    std::array<int16_t, 5> bins;

    bool operator==(const RatioHashKey &o) const {
        return bins == o.bins;
    }
};

struct RatioHashKeyHasher {
    size_t operator()(const RatioHashKey &k) const {
        // Simple FNV-1a style mix over 5 int16 bins.
        size_t h = (size_t)1469598103934665603ULL;
        for (int16_t b : k.bins) {
            h ^= (size_t)(uint16_t)b;
            h *= (size_t)1099511628211ULL;
        }
        return h;
    }
};

using CatalogPattern = std::array<int16_t, 4>;

struct CatalogPatternEntry {
    CatalogPattern         stars;
    std::array<stfloat, 6> edges;
};

using PatternHash =
    std::unordered_multimap<RatioHashKey, CatalogPatternEntry, RatioHashKeyHasher>;

struct PatternEdge {
    int a;
    int b;
};

static constexpr PatternEdge kPatternEdges[6] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
};

static constexpr int kPatternEdgeIndexByPair[4][4] = {
    {-1, 0, 1, 2},
    {0, -1, 3, 4},
    {1, 3, -1, 5},
    {2, 4, 5, -1},
};

static std::array<stfloat, 6> tetraEdges(const std::array<Vec3, 4> &vectors) {
    std::array<stfloat, 6> edges{};
    for (int e = 0; e < 6; ++e) {
        edges[e] = angleUnit(vectors[kPatternEdges[e].a], vectors[kPatternEdges[e].b]);
    }
    return edges;
}

static bool tetraHashKeyFromEdges(const std::array<stfloat, 6> &edges,
                                  stfloat patternMaxError,
                                  RatioHashKey *outKey,
                                  stfloat *largestEdgeOut = nullptr) {
    if (patternMaxError <= (stfloat)0.0) return false;
    int largestEdgeIndex = 0;
    for (int i = 1; i < 6; ++i) {
        if (edges[i] > edges[largestEdgeIndex]) largestEdgeIndex = i;
    }

    stfloat largestEdge = edges[largestEdgeIndex];
    if (largestEdge <= (stfloat)0.0) return false;

    std::array<stfloat, 5> ratios{};
    int r = 0;
    for (int i = 0; i < 6; ++i) {
        if (i == largestEdgeIndex) continue;
        ratios[r++] = edges[i] / largestEdge;
    }
    std::sort(ratios.begin(), ratios.end());

    RatioHashKey key{};
    for (int i = 0; i < 5; ++i) {
        int bin = (int)std::floor((double)(ratios[i] / patternMaxError));
        if (bin < 0 || bin > (int)std::numeric_limits<int16_t>::max()) return false;
        key.bins[i] = (int16_t)bin;
    }

    if (outKey) *outKey = key;
    if (largestEdgeOut) *largestEdgeOut = largestEdge;
    return true;
}

static bool tetraHashKey(const std::array<Vec3, 4> &vectors,
                         stfloat patternMaxError,
                         RatioHashKey *outKey,
                         stfloat *largestEdgeOut = nullptr,
                         std::array<stfloat, 6> *edgesOut = nullptr) {
    const std::array<stfloat, 6> edges = tetraEdges(vectors);
    if (edgesOut) *edgesOut = edges;
    return tetraHashKeyFromEdges(edges, patternMaxError, outKey, largestEdgeOut);
}

static bool mappingHasUniqueIndices(const std::array<int, 4> &map) {
    return map[0] != map[1] && map[0] != map[2] && map[0] != map[3]
        && map[1] != map[2] && map[1] != map[3]
        && map[2] != map[3];
}

static bool chiralityCompatible(const std::array<Vec3, 4> &imagePattern,
                                const std::array<Vec3, 4> &catalogPattern,
                                const std::array<int, 4> &perm) {
    const stfloat imageDet = imagePattern[0].cross(imagePattern[1]) * imagePattern[2];
    const stfloat catalogDet =
        catalogPattern[perm[0]].cross(catalogPattern[perm[1]]) * catalogPattern[perm[2]];

    // Degenerate triples carry little chirality information; accept them.
    if (stfabs(imageDet) < (stfloat)1e-9 || stfabs(catalogDet) < (stfloat)1e-9) {
        return true;
    }

    return (imageDet > (stfloat)0.0) == (catalogDet > (stfloat)0.0);
}

static bool edgesMatchPermutation(const std::array<stfloat, 6> &imageEdges,
                                  const std::array<stfloat, 6> &catalogEdges,
                                  const std::array<int, 4> &perm,
                                  stfloat tolerance) {
    for (int e = 0; e < 6; ++e) {
        const int ia = kPatternEdges[e].a;
        const int ib = kPatternEdges[e].b;
        const int ca = perm[ia];
        const int cb = perm[ib];
        const int edgeIdx = kPatternEdgeIndexByPair[ca][cb];
        if (edgeIdx < 0) return false;
        const stfloat catalogEdge = catalogEdges[edgeIdx];
        if (stfabs(catalogEdge - imageEdges[e]) > tolerance) return false;
    }
    return true;
}

struct CatalogHashCache {
    const Catalog *catalogPtr = nullptr;
    size_t         catalogSize = 0;
    stfloat        patternMaxError = (stfloat)0.0;
    int            catalogPatternStars = 0;
    long           patternCutoff = 0;
    stfloat        maxCatalogPatternEdge = (stfloat)0.0;
    PatternHash    hash;
};

static size_t estimatePatternCount(int n, long cutoff) {
    if (n < 4) return 0;
    long double total =
        (long double)n * (long double)(n - 1) * (long double)(n - 2) * (long double)(n - 3)
        / (long double)24.0;
    if (cutoff > 0 && total > (long double)cutoff) {
        total = (long double)cutoff;
    }
    if (total <= (long double)0.0) return 0;
    if (total > (long double)std::numeric_limits<size_t>::max()) {
        return std::numeric_limits<size_t>::max();
    }
    return (size_t)total;
}

static const PatternHash &buildOrGetCatalogHash(const Catalog &catalog,
                                                stfloat patternMaxError,
                                                int catalogPatternStars,
                                                long patternCutoff,
                                                stfloat maxCatalogPatternEdge) {
    static thread_local CatalogHashCache cache;

    const bool cacheHit =
        cache.catalogPtr == &catalog
        && cache.catalogSize == catalog.size()
        && cache.patternMaxError == patternMaxError
        && cache.catalogPatternStars == catalogPatternStars
        && cache.patternCutoff == patternCutoff
        && cache.maxCatalogPatternEdge == maxCatalogPatternEdge;

    if (cacheHit) return cache.hash;

    cache.catalogPtr = &catalog;
    cache.catalogSize = catalog.size();
    cache.patternMaxError = patternMaxError;
    cache.catalogPatternStars = catalogPatternStars;
    cache.patternCutoff = patternCutoff;
    cache.maxCatalogPatternEdge = maxCatalogPatternEdge;
    cache.hash.clear();

    std::vector<int> catalogOrder(catalog.size());
    std::iota(catalogOrder.begin(), catalogOrder.end(), 0);
    std::stable_sort(catalogOrder.begin(), catalogOrder.end(),
                     [&catalog](int a, int b) {
                         if (catalog[a].magnitude != catalog[b].magnitude)
                             return catalog[a].magnitude < catalog[b].magnitude;
                         return a < b;
                     });
    catalogOrder.resize(catalogPatternStars);

    const size_t reserveCount = estimatePatternCount(catalogPatternStars, patternCutoff);
    if (reserveCount > 0) {
        cache.hash.reserve(reserveCount);
    }

    long generatedPatterns = 0;
    for (int a = 0; a < catalogPatternStars - 3; ++a) {
        for (int b = a + 1; b < catalogPatternStars - 2; ++b) {
            for (int c = b + 1; c < catalogPatternStars - 1; ++c) {
                for (int d = c + 1; d < catalogPatternStars; ++d) {
                    if (patternCutoff > 0 && generatedPatterns >= patternCutoff) break;

                    const int ia = catalogOrder[a];
                    const int ib = catalogOrder[b];
                    const int ic = catalogOrder[c];
                    const int id = catalogOrder[d];

                    std::array<Vec3, 4> catVectors = {
                        catalog[ia].spatial, catalog[ib].spatial,
                        catalog[ic].spatial, catalog[id].spatial
                    };

                    RatioHashKey key{};
                    stfloat largestEdge = (stfloat)0.0;
                    std::array<stfloat, 6> catEdges{};
                    if (!tetraHashKey(catVectors, patternMaxError, &key,
                                      &largestEdge, &catEdges))
                        continue;

                    // Patterns larger than the image FOV are unlikely to appear.
                    if (largestEdge > maxCatalogPatternEdge) continue;

                    cache.hash.emplace(key, CatalogPatternEntry{
                        {(int16_t)ia, (int16_t)ib, (int16_t)ic, (int16_t)id},
                        catEdges
                    });
                    ++generatedPatterns;
                }
                if (patternCutoff > 0 && generatedPatterns >= patternCutoff) break;
            }
            if (patternCutoff > 0 && generatedPatterns >= patternCutoff) break;
        }
        if (patternCutoff > 0 && generatedPatterns >= patternCutoff) break;
    }

    return cache.hash;
}

}  // namespace

// ═══════════════════════════════════════════════════════════════════════════
//  Tetra3StarIdAlgorithm
// ═══════════════════════════════════════════════════════════════════════════

StarIdentifiers Tetra3StarIdAlgorithm::identify(
    const unsigned char *database,
    const Stars          &stars,
    const Catalog        &catalog,
    const Camera         &camera) const {

    StarIdentifiers identified;

    if (stars.size() < 4 || catalog.size() < 4 || patternMaxError_ <= (stfloat)0.0) {
        std::cerr << "Tetra3: insufficient stars/catalog or invalid parameters." << std::endl;
        return identified;
    }

    // Precompute normalized spatial vectors for all image stars.
    std::vector<Vec3> imageVectors(stars.size());
    for (size_t i = 0; i < stars.size(); ++i) {
        imageVectors[i] = camera.cameraToSpatial(stars[i].position).normalize();
    }

    // Tetra3 searches the brightest stars first.
    std::vector<int> imageOrder(stars.size());
    std::iota(imageOrder.begin(), imageOrder.end(), 0);
    std::stable_sort(imageOrder.begin(), imageOrder.end(),
                     [&stars](int a, int b) {
                         if (stars[a].magnitude != stars[b].magnitude)
                             return stars[a].magnitude < stars[b].magnitude;
                         return a < b;
                     });

    int imagePatternStars = (int)imageOrder.size();
    if (patternCheckingStars_ > 0) {
        imagePatternStars = std::min(imagePatternStars, patternCheckingStars_);
    }
    imageOrder.resize(imagePatternStars);

    int catalogPatternStars = (int)catalog.size();
    if (catalogPatternStars_ > 0) {
        catalogPatternStars = std::min(catalogPatternStars, catalogPatternStars_);
    }
    catalogPatternStars = std::min(catalogPatternStars,
                                   (int)std::numeric_limits<int16_t>::max());
    if (catalogPatternStars < 4) {
        std::cerr << "Tetra3: insufficient catalog stars for pattern hash." << std::endl;
        return identified;
    }

    const stfloat maxCatalogPatternEdge = camera.fov() * (stfloat)1.1;
    const PatternHash &catalogHash =
        buildOrGetCatalogHash(catalog, patternMaxError_, catalogPatternStars,
                              patternCutoff_, maxCatalogPatternEdge);

    if (catalogHash.empty()) {
        std::cerr << "Tetra3: catalog hash is empty." << std::endl;
        return identified;
    }

    // Try each image tetra pattern (bright stars first), query nearby hash bins,
    // and accept the first unique catalog match.
    for (int a = 0; a < imagePatternStars - 3; ++a) {
        for (int b = a + 1; b < imagePatternStars - 2; ++b) {
            for (int c = b + 1; c < imagePatternStars - 1; ++c) {
                for (int d = c + 1; d < imagePatternStars; ++d) {
                    std::array<int, 4> imageIdx = {
                        imageOrder[a], imageOrder[b], imageOrder[c], imageOrder[d]
                    };
                    std::array<Vec3, 4> imagePattern = {
                        imageVectors[imageIdx[0]], imageVectors[imageIdx[1]],
                        imageVectors[imageIdx[2]], imageVectors[imageIdx[3]]
                    };

                    RatioHashKey imageKey{};
                    stfloat largestImageEdge = (stfloat)0.0;
                    if (!tetraHashKey(imagePattern, patternMaxError_,
                                      &imageKey, &largestImageEdge))
                        continue;

                    if (largestImageEdge > maxCatalogPatternEdge) continue;

                    std::array<stfloat, 6> imageEdges = tetraEdges(imagePattern);

                    bool haveMatch = false;
                    bool duplicate = false;
                    std::array<int, 4> matchedCatalog = {-1, -1, -1, -1};

                    for (int o0 = -1; o0 <= 1 && !duplicate; ++o0) {
                        int b0 = imageKey.bins[0] + o0;
                        if (b0 < 0 || b0 > (int)std::numeric_limits<int16_t>::max()) continue;
                        for (int o1 = -1; o1 <= 1 && !duplicate; ++o1) {
                            int b1 = imageKey.bins[1] + o1;
                            if (b1 < 0 || b1 > (int)std::numeric_limits<int16_t>::max()) continue;
                            for (int o2 = -1; o2 <= 1 && !duplicate; ++o2) {
                                int b2 = imageKey.bins[2] + o2;
                                if (b2 < 0 || b2 > (int)std::numeric_limits<int16_t>::max()) continue;
                                for (int o3 = -1; o3 <= 1 && !duplicate; ++o3) {
                                    int b3 = imageKey.bins[3] + o3;
                                    if (b3 < 0 || b3 > (int)std::numeric_limits<int16_t>::max()) continue;
                                    for (int o4 = -1; o4 <= 1 && !duplicate; ++o4) {
                                        int b4 = imageKey.bins[4] + o4;
                                        if (b4 < 0 || b4 > (int)std::numeric_limits<int16_t>::max()) continue;

                                        RatioHashKey queryKey = {{
                                            (int16_t)b0, (int16_t)b1, (int16_t)b2,
                                            (int16_t)b3, (int16_t)b4
                                        }};

                                        auto range = catalogHash.equal_range(queryKey);
                                        for (auto it = range.first; it != range.second; ++it) {
                                            const CatalogPatternEntry &candidate = it->second;

                                            std::array<Vec3, 4> catalogPattern = {
                                                catalog[candidate.stars[0]].spatial,
                                                catalog[candidate.stars[1]].spatial,
                                                catalog[candidate.stars[2]].spatial,
                                                catalog[candidate.stars[3]].spatial
                                            };

                                            std::array<int, 4> perm = {0, 1, 2, 3};
                                            do {
                                                if (!chiralityCompatible(imagePattern,
                                                                         catalogPattern,
                                                                         perm))
                                                    continue;
                                                if (!edgesMatchPermutation(imageEdges,
                                                                           candidate.edges,
                                                                           perm,
                                                                           tolerance_))
                                                    continue;

                                                std::array<int, 4> candidateMap = {
                                                    candidate.stars[perm[0]],
                                                    candidate.stars[perm[1]],
                                                    candidate.stars[perm[2]],
                                                    candidate.stars[perm[3]]
                                                };
                                                if (!mappingHasUniqueIndices(candidateMap)) {
                                                    continue;
                                                }

                                                if (!haveMatch) {
                                                    matchedCatalog = candidateMap;
                                                    haveMatch = true;
                                                } else if (matchedCatalog != candidateMap) {
                                                    duplicate = true;
                                                    break;
                                                }
                                            } while (std::next_permutation(perm.begin(), perm.end()));

                                            if (duplicate) break;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (haveMatch && !duplicate) {
                        identified.push_back(StarIdentifier(imageIdx[0], matchedCatalog[0]));
                        identified.push_back(StarIdentifier(imageIdx[1], matchedCatalog[1]));
                        identified.push_back(StarIdentifier(imageIdx[2], matchedCatalog[2]));
                        identified.push_back(StarIdentifier(imageIdx[3], matchedCatalog[3]));

                        if (database) {
                            MultiDatabase multiDb(database);
                            const unsigned char *dbBuf =
                                multiDb.subDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
                            if (dbBuf) {
                                DeserializeContext des(dbBuf);
                                PairDistanceKVectorDatabase vectorDb(&des);
                                identifyRemainingStarsPairDistance(
                                    identified, stars, vectorDb, catalog, camera, tolerance_);
                            }
                        }

                        return identified;
                    }
                }
            }
        }
    }

    std::cerr << "Tetra3: no unique match found." << std::endl;
    return identified;
}
