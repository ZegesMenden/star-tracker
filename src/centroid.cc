#include "centroid.hpp"

#include <stack>
#include <unordered_set>
#include <vector>

uint8_t CenterOfGravityCentroid::threshold(uint8_t *image, size_t imwidth, size_t imheight) {

    size_t totalPixels = imwidth * imheight;
    unsigned long totalMag = 0;
    stfloat sqTotalMag = 0;

    for (size_t i = 0; i < totalPixels; i++) {
        totalMag += image[i];
        sqTotalMag += (stfloat)image[i] * (stfloat)image[i];
    }

    stfloat mean = (stfloat)totalMag / (stfloat)totalPixels;
    stfloat variance = (sqTotalMag / (stfloat)totalPixels) - (mean * mean);
    stfloat stddev = stsqrt(variance);

    stfloat raw = mean + stddev * sigma_;
    return (raw > (stfloat)255.0) ? 255 : (uint8_t)raw;
}

std::vector<Centroid> CenterOfGravityCentroid::compute(uint8_t *image, size_t imwidth, size_t imheight) {

    std::vector<Centroid> result;
    std::unordered_set<long> visited;

    uint8_t cutoff = threshold(image, imwidth, imheight);

    long totalPixels = (long)(imwidth * imheight);

    for (long i = 0; i < totalPixels; i++) {
        if (image[i] < cutoff || visited.count(i) != 0) {
            continue;
        }

        // Found a new star region — flood-fill iteratively using a stack
        stfloat xCoordMagSum = 0;
        stfloat yCoordMagSum = 0;
        long magSum = 0;

        int xMin = (int)(i % imwidth);
        int xMax = xMin;
        int yMin = (int)(i / imwidth);
        int yMax = yMin;
        bool isValid = true;
        size_t starPixelCount = 0;

        std::stack<long> frontier;
        frontier.push(i);

        while (!frontier.empty()) {
            long idx = frontier.top();
            frontier.pop();

            // Bounds check, brightness check, and duplicate check
            if (idx < 0 || idx >= totalPixels) {
                continue;
            }
            if (image[idx] < cutoff) {
                continue;
            }
            if (visited.count(idx) != 0) {
                continue;
            }

            visited.insert(idx);
            starPixelCount++;

            int x = (int)(idx % imwidth);
            int y = (int)(idx / imwidth);

            // Mark star as invalid if it touches the image border
            if (x == 0 || x == (int)imwidth - 1 || y == 0 || y == (int)imheight - 1) {
                isValid = false;
            }

            // Update bounding box
            if (x > xMax) xMax = x;
            if (x < xMin) xMin = x;
            if (y > yMax) yMax = y;
            if (y < yMin) yMin = y;

            // Accumulate weighted sums for center-of-gravity calculation
            magSum += image[idx];
            xCoordMagSum += (stfloat)x * image[idx];
            yCoordMagSum += (stfloat)y * image[idx];

            // Push 4-connected neighbors onto the stack
            if (x < (int)imwidth - 1)   frontier.push(idx + 1); // right
            if (x > 0)                  frontier.push(idx - 1); // left
            if (y < (int)imheight - 1)  frontier.push(idx + (long)imwidth); // down
            if (y > 0)                  frontier.push(idx - (long)imwidth); // up
        }

        if (isValid && magSum > 0) {
            stfloat xCoord = xCoordMagSum / (stfloat)magSum + (stfloat)0.5;
            stfloat yCoord = yCoordMagSum / (stfloat)magSum + (stfloat)0.5;
            stfloat radiusX = (stfloat)(xMax - xMin + 1) / (stfloat)2.0;
            stfloat radiusY = (stfloat)(yMax - yMin + 1) / (stfloat)2.0;

            Centroid c;
            c.x = xCoord;
            c.y = yCoord;
            c.radiusX = radiusX;
            c.radiusY = radiusY;
            c.numPixels = starPixelCount;
            result.push_back(c);
        }
    }

    return result;
}
