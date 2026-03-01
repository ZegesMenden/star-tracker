#include "camera.hpp"

#include <cassert>

Vec2 Camera::spatialToCamera(const Vec3 &v) const {
    assert(v.x > 0);
    stfloat f = focalLength / v.x;
    return { -v.y * f + xCenter, -v.z * f + yCenter };
}

Vec3 Camera::cameraToSpatial(const Vec2 &v) const {
    stfloat xPixel = -v.x + xCenter;
    stfloat yPixel = -v.y + yCenter;
    return { (stfloat)1.0, xPixel / focalLength, yPixel / focalLength };
}

bool Camera::inSensor(const Vec2 &v) const {
    return v.x >= 0 && v.x <= xResolution
        && v.y >= 0 && v.y <= yResolution;
}

stfloat Camera::fov() const {
    return statan2((stfloat)xResolution / (stfloat)2.0, focalLength) * (stfloat)2.0;
}

stfloat fovToFocalLength(stfloat xFov, stfloat xResolution) {
    return xResolution / (stfloat)2.0 / sttan(xFov / (stfloat)2.0);
}
