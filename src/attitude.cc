#include "attitude.hpp"

#include <cmath>
#include <cassert>

// ═══════════════════════════════════════════════════════════════════════════
//  QUEST internals
// ═══════════════════════════════════════════════════════════════════════════

static const stfloat kQuestEpsilon = (stfloat)0.0001;

/// Characteristic polynomial of the Davenport K-matrix.
/// f(λ) = (λ² − a)(λ² − b) − cλ + cσ − d
/// @see Eq. 19b in https://arc.aiaa.org/doi/pdf/10.2514/1.62549
static stfloat questCharPoly(stfloat lam,
                             stfloat a, stfloat b, stfloat c,
                             stfloat d, stfloat sigma) {
    return (lam*lam - a) * (lam*lam - b)
         - c * lam + c * sigma - d;
}

/// Derivative of the characteristic polynomial.
/// f'(λ) = 4λ³ − 2(a+b)λ − c
static stfloat questCharPolyPrime(stfloat lam,
                                  stfloat a, stfloat b, stfloat c) {
    return (stfloat)4.0*lam*lam*lam - (stfloat)2.0*(a + b)*lam - c;
}

/// Newton-Raphson iteration to find the largest eigenvalue.
static stfloat questEigenvalue(stfloat guess,
                               stfloat a, stfloat b, stfloat c,
                               stfloat d, stfloat sigma) {
    stfloat h;
    int maxIter = 1000;
    do {
        stfloat fp = questCharPolyPrime(guess, a, b, c);
        if (stfabs(fp) < (stfloat)1e-15) break;   // avoid divide-by-zero
        h = questCharPoly(guess, a, b, c, d, sigma) / fp;
        guess -= h;
    } while (stfabs(h) >= kQuestEpsilon && --maxIter > 0);
    return guess;
}

// ═══════════════════════════════════════════════════════════════════════════
//  QUEST algorithm
// ═══════════════════════════════════════════════════════════════════════════

Attitude QuestAlgorithm::estimate(const Camera         &camera,
                                  const Stars          &stars,
                                  const Catalog        &catalog,
                                  const StarIdentifiers &ids) const {
    if (ids.size() < 2) {
        return Attitude();   // unknown
    }

    // -- Build the attitude profile matrix  B = Σ wᵢ rᵢ bᵢᵀ
    //    where rᵢ = catalog (inertial) unit vector
    //          bᵢ = camera  (body)     unit vector
    Mat3 B = {0,0,0, 0,0,0, 0,0,0};
    stfloat weightSum = 0;

    for (const StarIdentifier &s : ids) {
        Vec3 body = camera.cameraToSpatial(stars[s.starIndex].position).normalize();
        const Vec3 &ref = catalog[s.catalogIndex].spatial;
        B = B + outerProduct(ref, body) * s.weight;
        weightSum += s.weight;
    }

    // S = B + Bᵀ
    Mat3 S = B + B.transpose();

    // σ = trace(B)
    stfloat sigma = B.trace();

    // Z = skew-symmetric part of B
    Vec3 Z = {
        B.at(1,2) - B.at(2,1),
        B.at(2,0) - B.at(0,2),
        B.at(0,1) - B.at(1,0),
    };

    // -- Coefficients of the characteristic polynomial
    stfloat delta = S.det();
    stfloat kappa = (S.inverse() * delta).trace();
    stfloat a     = sigma * sigma - kappa;
    stfloat b     = sigma * sigma + (Z * Z);      // Z·Z = |Z|²
    stfloat c     = delta + (Z * (S * Z));         // Z · S Z
    stfloat d     = Z * ((S * S) * Z);             // Z · S² Z

    // -- Newton-Raphson for the largest eigenvalue of the K-matrix
    stfloat eig = questEigenvalue(weightSum, a, b, c, d, sigma);

    // -- Compute the optimal quaternion from the eigenvalue
    stfloat alpha = eig*eig - sigma*sigma + kappa;
    stfloat beta  = eig - sigma;
    stfloat gamma = (eig + sigma) * alpha - delta;

    // X = (αI + βS + S²) Z
    Mat3 term = (kIdentityMat3 * alpha) + (S * beta) + (S * S);
    Vec3 X = term * Z;

    // Normalize (γ, X) to get unit quaternion
    stfloat norm = (stfloat)1.0 / stsqrt(gamma * gamma + X.magnitudeSq());
    X = X * norm;
    gamma *= norm;

    return Attitude(Quaternion(gamma, X.x, X.y, X.z));
}

// ═══════════════════════════════════════════════════════════════════════════
//  TRIAD algorithm
// ═══════════════════════════════════════════════════════════════════════════

/// Build a coordinate frame (3×3 matrix whose columns are orthonormal axes)
/// from two non-parallel direction vectors.
static Mat3 triadFrame(Vec3 v1, Vec3 v2) {
    Vec3 d1 = v1.normalize();
    Vec3 d2 = v1.cross(v2).normalize();
    Vec3 d3 = d1.cross(d2).normalize();
    return {
        d1.x, d2.x, d3.x,
        d1.y, d2.y, d3.y,
        d1.z, d2.z, d3.z,
    };
}

Attitude TriadAlgorithm::estimate(const Camera         &camera,
                                  const Stars          &stars,
                                  const Catalog        &catalog,
                                  const StarIdentifiers &ids) const {
    if (ids.size() < 2) {
        return Attitude();
    }

    const StarIdentifier &a = ids[0];
    const StarIdentifier &b = ids[ids.size() / 2];

    Mat3 bodyFrame = triadFrame(
        camera.cameraToSpatial(stars[a.starIndex].position),
        camera.cameraToSpatial(stars[b.starIndex].position));

    Mat3 refFrame = triadFrame(
        catalog[a.catalogIndex].spatial,
        catalog[b.catalogIndex].spatial);

    // attitude = bodyFrame × refFrameᵀ
    Mat3 dcm = bodyFrame * refFrame.transpose();

    // Convert DCM → quaternion via Shepperd's method
    stfloat tr = dcm.trace();
    stfloat qr, qi, qj, qk;
    if (tr > 0) {
        stfloat s = stsqrt(tr + (stfloat)1.0) * (stfloat)2.0;
        qr = (stfloat)0.25 * s;
        qi = (dcm.at(2,1) - dcm.at(1,2)) / s;
        qj = (dcm.at(0,2) - dcm.at(2,0)) / s;
        qk = (dcm.at(1,0) - dcm.at(0,1)) / s;
    } else if (dcm.at(0,0) > dcm.at(1,1) && dcm.at(0,0) > dcm.at(2,2)) {
        stfloat s = stsqrt((stfloat)1.0 + dcm.at(0,0) - dcm.at(1,1) - dcm.at(2,2)) * (stfloat)2.0;
        qr = (dcm.at(2,1) - dcm.at(1,2)) / s;
        qi = (stfloat)0.25 * s;
        qj = (dcm.at(0,1) + dcm.at(1,0)) / s;
        qk = (dcm.at(0,2) + dcm.at(2,0)) / s;
    } else if (dcm.at(1,1) > dcm.at(2,2)) {
        stfloat s = stsqrt((stfloat)1.0 + dcm.at(1,1) - dcm.at(0,0) - dcm.at(2,2)) * (stfloat)2.0;
        qr = (dcm.at(0,2) - dcm.at(2,0)) / s;
        qi = (dcm.at(0,1) + dcm.at(1,0)) / s;
        qj = (stfloat)0.25 * s;
        qk = (dcm.at(1,2) + dcm.at(2,1)) / s;
    } else {
        stfloat s = stsqrt((stfloat)1.0 + dcm.at(2,2) - dcm.at(0,0) - dcm.at(1,1)) * (stfloat)2.0;
        qr = (dcm.at(1,0) - dcm.at(0,1)) / s;
        qi = (dcm.at(0,2) + dcm.at(2,0)) / s;
        qj = (dcm.at(1,2) + dcm.at(2,1)) / s;
        qk = (stfloat)0.25 * s;
    }

    return Attitude(Quaternion(qr, qi, qj, qk));
}
