#include "star.hpp"

#include <cmath>
#include <algorithm>

// ── Vec2 ────────────────────────────────────────────────────────────────────

stfloat Vec2::magnitude()   const { return stsqrt(x * x + y * y); }
stfloat Vec2::magnitudeSq() const { return x * x + y * y; }
Vec2    Vec2::normalize()   const { stfloat m = magnitude(); return {x / m, y / m}; }

stfloat Vec2::operator*(const Vec2 &o) const { return x * o.x + y * o.y; }
Vec2    Vec2::operator*(stfloat s)      const { return {x * s, y * s}; }
Vec2    Vec2::operator-(const Vec2 &o) const { return {x - o.x, y - o.y}; }
Vec2    Vec2::operator+(const Vec2 &o) const { return {x + o.x, y + o.y}; }

// ── Vec3 ────────────────────────────────────────────────────────────────────

stfloat Vec3::magnitude()   const { return stsqrt(x * x + y * y + z * z); }
stfloat Vec3::magnitudeSq() const { return x * x + y * y + z * z; }
Vec3    Vec3::normalize()   const { stfloat m = magnitude(); return {x / m, y / m, z / m}; }

stfloat Vec3::operator*(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
Vec3    Vec3::operator*(stfloat s)      const { return {x * s, y * s, z * s}; }
Vec3    Vec3::operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
Vec3    Vec3::operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }

Vec3 Vec3::cross(const Vec3 &o) const {
    return {
        y * o.z - z * o.y,
        z * o.x - x * o.z,
        x * o.y - y * o.x,
    };
}

// ── Angle helpers ───────────────────────────────────────────────────────────

stfloat angle(const Vec3 &a, const Vec3 &b) {
    Vec3 an = a.normalize();
    Vec3 bn = b.normalize();
    return angleUnit(an, bn);
}

stfloat angleUnit(const Vec3 &a, const Vec3 &b) {
    stfloat dot = a * b;
    // clamp to [-1, 1] to avoid NaN from floating-point drift
    dot = std::max((stfloat)-1.0, std::min((stfloat)1.0, dot));
    return stacos(dot);
}

Vec3 sphericalToSpatial(stfloat ra, stfloat de) {
    return {
        stcos(ra) * stcos(de),
        stsin(ra) * stcos(de),
        stsin(de),
    };
}

// ── Mat3 ────────────────────────────────────────────────────────────────────

const Mat3 kIdentityMat3 = {1,0,0, 0,1,0, 0,0,1};

stfloat Mat3::at(int i, int j)  const { return x[3 * i + j]; }

Vec3 Mat3::column(int j) const { return {at(0,j), at(1,j), at(2,j)}; }
Vec3 Mat3::row(int i)    const { return {at(i,0), at(i,1), at(i,2)}; }

Mat3 Mat3::operator+(const Mat3 &o) const {
    Mat3 r;
    for (int n = 0; n < 9; n++) r.x[n] = x[n] + o.x[n];
    return r;
}

Mat3 Mat3::operator*(const Mat3 &o) const {
    Mat3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r.x[3*i+j] = at(i,0)*o.at(0,j) + at(i,1)*o.at(1,j) + at(i,2)*o.at(2,j);
    return r;
}

Vec3 Mat3::operator*(const Vec3 &v) const {
    return {
        v.x*at(0,0) + v.y*at(0,1) + v.z*at(0,2),
        v.x*at(1,0) + v.y*at(1,1) + v.z*at(1,2),
        v.x*at(2,0) + v.y*at(2,1) + v.z*at(2,2),
    };
}

Mat3 Mat3::operator*(stfloat s) const {
    Mat3 r;
    for (int n = 0; n < 9; n++) r.x[n] = x[n] * s;
    return r;
}

Mat3 Mat3::transpose() const {
    return {
        at(0,0), at(1,0), at(2,0),
        at(0,1), at(1,1), at(2,1),
        at(0,2), at(1,2), at(2,2),
    };
}

stfloat Mat3::trace() const { return at(0,0) + at(1,1) + at(2,2); }

stfloat Mat3::det() const {
    return at(0,0) * (at(1,1)*at(2,2) - at(2,1)*at(1,2))
         - at(0,1) * (at(1,0)*at(2,2) - at(2,0)*at(1,2))
         + at(0,2) * (at(1,0)*at(2,1) - at(2,0)*at(1,1));
}

Mat3 Mat3::inverse() const {
    stfloat d = det();
    stfloat s = (stfloat)1.0 / d;
    Mat3 r = {
        at(1,1)*at(2,2) - at(1,2)*at(2,1),  at(0,2)*at(2,1) - at(0,1)*at(2,2),  at(0,1)*at(1,2) - at(0,2)*at(1,1),
        at(1,2)*at(2,0) - at(1,0)*at(2,2),  at(0,0)*at(2,2) - at(0,2)*at(2,0),  at(0,2)*at(1,0) - at(0,0)*at(1,2),
        at(1,0)*at(2,1) - at(1,1)*at(2,0),  at(0,1)*at(2,0) - at(0,0)*at(2,1),  at(0,0)*at(1,1) - at(0,1)*at(1,0),
    };
    return r * s;
}

Mat3 outerProduct(const Vec3 &a, const Vec3 &b) {
    return {
        a.x*b.x, a.x*b.y, a.x*b.z,
        a.y*b.x, a.y*b.y, a.y*b.z,
        a.z*b.x, a.z*b.y, a.z*b.z,
    };
}

// ── Quaternion ──────────────────────────────────────────────────────────────

Quaternion::Quaternion(const Vec3 &axis, stfloat theta)
    : real(stcos(theta / (stfloat)2.0)),
      i(axis.x * stsin(theta / (stfloat)2.0)),
      j(axis.y * stsin(theta / (stfloat)2.0)),
      k(axis.z * stsin(theta / (stfloat)2.0)) {}

Quaternion Quaternion::operator*(const Quaternion &o) const {
    return Quaternion(
        real*o.real - i*o.i - j*o.j - k*o.k,
        real*o.i + o.real*i + j*o.k - k*o.j,
        real*o.j + o.real*j + k*o.i - i*o.k,
        real*o.k + o.real*k + i*o.j - j*o.i);
}

Quaternion Quaternion::conjugate() const {
    return Quaternion(real, -i, -j, -k);
}

Vec3 Quaternion::vector() const { return {i, j, k}; }

Vec3 Quaternion::rotate(const Vec3 &v) const {
    // q * v_pure * q*
    Quaternion vq(0, v.x, v.y, v.z);
    Quaternion result = (*this) * vq * conjugate();
    return result.vector();
}

stfloat Quaternion::angle() const {
    stfloat r = std::max((stfloat)-1.0, std::min((stfloat)1.0, real));
    return stacos(r) * (stfloat)2.0;
}

stfloat Quaternion::smallestAngle() const {
    stfloat a = angle();
    return a > (stfloat)M_PI ? (stfloat)(2.0 * M_PI) - a : a;
}

bool Quaternion::isUnit(stfloat tolerance) const {
    return stfabs(real*real + i*i + j*j + k*k - (stfloat)1.0) < tolerance;
}

Quaternion Quaternion::canonicalize() const {
    return real >= 0 ? *this : Quaternion(-real, -i, -j, -k);
}

EulerAngles Quaternion::toSpherical() const {
    stfloat ra   = statan2((stfloat)2.0 * (-real*k + i*j), (stfloat)1.0 - (stfloat)2.0*(j*j + k*k));
    if (ra < 0) ra += (stfloat)(2.0 * M_PI);
    stfloat de   = -stasin(std::max((stfloat)-1.0,
                           std::min((stfloat)1.0, (stfloat)2.0 * (-real*j - i*k))));
    stfloat roll = -statan2((stfloat)2.0 * (-real*i + j*k), (stfloat)1.0 - (stfloat)2.0*(i*i + j*j));
    if (roll < 0) roll += (stfloat)(2.0 * M_PI);
    return EulerAngles(ra, de, roll);
}

// ── Attitude ────────────────────────────────────────────────────────────────

Quaternion Attitude::getQuaternion() const {
    return quat_;
}

EulerAngles Attitude::toSpherical() const {
    return quat_.toSpherical();
}

Vec3 Attitude::rotate(const Vec3 &v) const {
    return quat_.rotate(v);
}

Quaternion sphericalToQuaternion(stfloat ra, stfloat dec, stfloat roll) {
    Quaternion a(Vec3{0, 0, 1}, ra);
    Quaternion b(Vec3{0, 1, 0}, -dec);
    Quaternion c(Vec3{1, 0, 0}, -roll);
    return (a * b * c).conjugate();
}

// ── CatalogStar constructors ────────────────────────────────────────────────

CatalogStar::CatalogStar(stfloat ra, stfloat de, int magnitude, int16_t name)
    : spatial(sphericalToSpatial(ra, de)), magnitude(magnitude), name(name) {}

CatalogStar::CatalogStar(Vec3 spatial, int magnitude, int16_t name)
    : spatial(spatial), magnitude(magnitude), name(name) {}
