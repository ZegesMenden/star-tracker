#pragma once

#include "stfloat.h"

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>

// ── Vector types ────────────────────────────────────────────────────────────

struct Vec2 {
    stfloat x;
    stfloat y;

    stfloat magnitude()    const;
    stfloat magnitudeSq()  const;
    Vec2    normalize()    const;

    stfloat operator*(const Vec2 &) const;   // dot
    Vec2    operator*(stfloat)      const;
    Vec2    operator-(const Vec2 &) const;
    Vec2    operator+(const Vec2 &) const;
};

struct Vec3 {
    stfloat x;
    stfloat y;
    stfloat z;

    stfloat magnitude()    const;
    stfloat magnitudeSq()  const;
    Vec3    normalize()    const;

    stfloat operator*(const Vec3 &) const;   // dot
    Vec3    operator*(stfloat)      const;
    Vec3    operator-(const Vec3 &) const;
    Vec3    operator+(const Vec3 &) const;
    Vec3    cross(const Vec3 &)     const;
};

// ── 3×3 Matrix ──────────────────────────────────────────────────────────────

class Mat3 {
public:
    stfloat x[9];   // row-major: x[3*row + col]

    stfloat at(int i, int j) const;

    Mat3 operator+(const Mat3 &) const;
    Mat3 operator*(const Mat3 &) const;
    Vec3  operator*(const Vec3 &) const;
    Mat3 operator*(stfloat)      const;

    Mat3 transpose() const;
    stfloat trace()  const;
    stfloat det()    const;
    Mat3 inverse()   const;

    Vec3 column(int j) const;
    Vec3 row(int i)    const;
};

extern const Mat3 kIdentityMat3;

/// Outer product of two Vec3 → 3×3 matrix.
Mat3 outerProduct(const Vec3 &a, const Vec3 &b);

// ── Quaternion ──────────────────────────────────────────────────────────────

struct EulerAngles {
    stfloat ra;     // right ascension (yaw)
    stfloat de;     // declination (pitch)
    stfloat roll;   // roll

    EulerAngles(stfloat ra, stfloat de, stfloat roll)
        : ra(ra), de(de), roll(roll) {}
};

/// Unit quaternion representing a 3-D rotation.
class Quaternion {
public:
    stfloat real;
    stfloat i, j, k;

    Quaternion() : real(1), i(0), j(0), k(0) {}
    Quaternion(stfloat real, stfloat i, stfloat j, stfloat k)
        : real(real), i(i), j(j), k(k) {}

    /// Construct from an axis and angle (radians).
    Quaternion(const Vec3 &axis, stfloat theta);

    Quaternion operator*(const Quaternion &) const;
    Quaternion conjugate()   const;

    Vec3       vector()      const;
    Vec3       rotate(const Vec3 &) const;

    stfloat    angle()       const;
    stfloat    smallestAngle() const;
    bool       isUnit(stfloat tolerance = (stfloat)1e-6) const;
    Quaternion canonicalize() const;

    EulerAngles toSpherical() const;
};

/// Wrapper that stores an orientation estimate (or "unknown").
class Attitude {
public:
    Attitude()
        : type_(Type::Unknown) {}
    explicit Attitude(const Quaternion &q)
        : type_(Type::QuatType), quat_(q) {}

    bool        isKnown()       const { return type_ != Type::Unknown; }
    Quaternion  getQuaternion()  const;
    EulerAngles toSpherical()   const;
    Vec3        rotate(const Vec3 &v) const;

private:
    enum class Type { Unknown, QuatType };
    Type       type_;
    Quaternion quat_;
};

Quaternion sphericalToQuaternion(stfloat ra, stfloat dec, stfloat roll);

// ── Angle utilities ─────────────────────────────────────────────────────────

/// Angle between two arbitrary vectors (normalizes internally).
stfloat angle(const Vec3 &a, const Vec3 &b);

/// Angle between two unit vectors (skip normalization).
stfloat angleUnit(const Vec3 &a, const Vec3 &b);

/// Convert spherical (ra, de in radians) → unit vector.
Vec3 sphericalToSpatial(stfloat ra, stfloat de);

// ── Star types ──────────────────────────────────────────────────────────────

/// A catalog star with a known 3-D unit-sphere position and identifier.
struct CatalogStar {
    Vec3     spatial;    // unit vector on the celestial sphere
    int      magnitude;  // magnitude × 100  (lower = brighter)
    int16_t  name;       // unique catalog identifier

    CatalogStar() = default;
    CatalogStar(stfloat ra, stfloat de, int magnitude, int16_t name);
    CatalogStar(Vec3 spatial, int magnitude, int16_t name);
};

/// A centroid detected in an image (not yet identified).
struct Star {
    Vec2    position;
    stfloat radiusX;
    stfloat radiusY;
    int     magnitude;

    Star() : position{0, 0}, radiusX(0), radiusY(0), magnitude(0) {}
    Star(stfloat x, stfloat y, stfloat rx, stfloat ry, int mag)
        : position{x, y}, radiusX(rx), radiusY(ry), magnitude(mag) {}
    Star(stfloat x, stfloat y, stfloat r)
        : Star(x, y, r, r, 0) {}
};

/// Records which centroid (starIndex) maps to which catalog entry (catalogIndex).
struct StarIdentifier {
    int     starIndex;
    int     catalogIndex;
    stfloat weight;

    StarIdentifier(int si, int ci, stfloat w = (stfloat)1.0)
        : starIndex(si), catalogIndex(ci), weight(w) {}

    bool operator==(const StarIdentifier &o) const {
        return starIndex == o.starIndex && catalogIndex == o.catalogIndex;
    }
};

// ── Convenience typedefs ────────────────────────────────────────────────────

using Catalog         = std::vector<CatalogStar>;
using Stars           = std::vector<Star>;
using StarIdentifiers = std::vector<StarIdentifier>;
