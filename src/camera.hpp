#pragma once

#include "star.hpp"

/// Pinhole camera model.  Enough information to project between 3-D spatial
/// vectors and 2-D image coordinates.
///
/// Convention: X is the depth (boresight) axis.  A spatial vector (1, 0, 0) maps
/// to the principal point (xCenter, yCenter) on the sensor.
class Camera {
public:
    Camera(stfloat focalLength, stfloat xCenter, stfloat yCenter,
           int xResolution, int yResolution)
        : focalLength(focalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution) {}

    Camera(stfloat focalLength, int xResolution, int yResolution)
        : Camera(focalLength,
                 (stfloat)xResolution / (stfloat)2.0,
                 (stfloat)yResolution / (stfloat)2.0,
                 xResolution, yResolution) {}

    /// Project a 3-D unit vector onto the 2-D sensor plane.
    Vec2 spatialToCamera(const Vec3 &v) const;

    /// Back-project a 2-D pixel coordinate to a 3-D direction (x-component = 1).
    Vec3 cameraToSpatial(const Vec2 &v) const;

    /// Is the pixel inside the sensor bounds?
    bool inSensor(const Vec2 &v) const;

    int     getXResolution()  const { return xResolution; }
    int     getYResolution()  const { return yResolution; }
    stfloat getFocalLength()  const { return focalLength; }
    stfloat fov()             const;

private:
    stfloat focalLength;
    stfloat xCenter, yCenter;
    int     xResolution, yResolution;
};

stfloat fovToFocalLength(stfloat xFov, stfloat xResolution);
