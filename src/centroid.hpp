#pragma once

#include "stfloat.h"

#include <cstdint>
#include <cstddef>
#include <vector>

struct Centroid {
    stfloat x;
    stfloat y;
    stfloat radiusX;
    stfloat radiusY;
    size_t numPixels;
};

class CentroidAlgorithm {

public:
    virtual ~CentroidAlgorithm() = default;

    virtual std::vector<Centroid> compute(uint8_t *image, size_t imwidth, size_t imheight) = 0;
};

class CenterOfGravityCentroid : public CentroidAlgorithm {

public:
    /// Construct with a configurable sigma threshold multiplier (default: 5).
    explicit CenterOfGravityCentroid(stfloat sigma = 5.0) : sigma_(sigma) {}

    std::vector<Centroid> compute(uint8_t *image, size_t imwidth, size_t imheight) override;

private:
    stfloat sigma_;
    uint8_t threshold(uint8_t *image, size_t imwidth, size_t imheight);
};