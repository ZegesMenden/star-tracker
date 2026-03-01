#pragma once

#include "star.hpp"
#include "camera.hpp"

#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
//  Attitude estimation algorithm base class
// ═══════════════════════════════════════════════════════════════════════════

/// Abstract base – different algorithms inherit and override `estimate()`.
class AttitudeEstimationAlgorithm {
public:
    virtual ~AttitudeEstimationAlgorithm() = default;

    /// Given identified stars, their centroids, the catalog, and a camera model,
    /// compute the optimal spacecraft attitude (orientation).
    virtual Attitude estimate(const Camera         &camera,
                              const Stars          &stars,
                              const Catalog        &catalog,
                              const StarIdentifiers &ids) const = 0;
};

// ═══════════════════════════════════════════════════════════════════════════
//  QUEST algorithm
// ═══════════════════════════════════════════════════════════════════════════

/// The QUEST (QUaternion ESTimator) algorithm.
///
/// Computes the optimal rotation quaternion that maps catalog (inertial)
/// vectors to camera (body) vectors by finding the largest eigenvalue of
/// the Davenport K-matrix via Newton-Raphson iteration on its
/// characteristic polynomial.
///
/// Reference: Shuster & Oh, "Three-Axis Attitude Determination from
/// Vector Observations", J. Guidance & Control 4(1), 1981.
class QuestAlgorithm : public AttitudeEstimationAlgorithm {
public:
    QuestAlgorithm() = default;

    Attitude estimate(const Camera         &camera,
                      const Stars          &stars,
                      const Catalog        &catalog,
                      const StarIdentifiers &ids) const override;
};

// ═══════════════════════════════════════════════════════════════════════════
//  TRIAD algorithm
// ═══════════════════════════════════════════════════════════════════════════

/// The TRIAD algorithm: a fast attitude estimator using only two star pairs.
/// Less accurate than QUEST but very cheap.
class TriadAlgorithm : public AttitudeEstimationAlgorithm {
public:
    TriadAlgorithm() = default;

    Attitude estimate(const Camera         &camera,
                      const Stars          &stars,
                      const Catalog        &catalog,
                      const StarIdentifiers &ids) const override;
};
