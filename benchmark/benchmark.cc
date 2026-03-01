// benchmark.cc — Run the full star-tracker pipeline on a single image and
// print the resulting attitude (RA, DEC, roll) for comparison against a known
// solution (e.g. from astrometry.net).
//
// Usage:
//   benchmark --image <path>         Path to an image (JPEG/PNG/BMP/…)
//             --database <path>      Path to the serialized MultiDatabase file
//             --catalog <path>       Path to the bright-star-catalog TSV
//             --focal-length <px>    Focal length in pixels
//                                     OR
//             --pixel-size <µm>      Pixel pitch in micrometers
//             --focal-length-mm <mm> Focal length in millimeters
//
// Output (stdout):
//   RA_deg DEC_deg ROLL_deg  NUM_IDENTIFIED  NUM_CENTROIDS

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// stb_image — include the implementation here
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNG
#define STBI_ONLY_BMP
#include "stb_image.h"

// star-tracker headers
#include "centroid.hpp"
#include "star.hpp"
#include "camera.hpp"
#include "database.hpp"
#include "starid.hpp"
#include "attitude.hpp"
#include "io.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void printUsage(const char *argv0) {
    std::fprintf(stderr,
        "Usage: %s --image <path> --database <path> --catalog <bsc.tsv>\n"
        "          --focal-length <px>  [--tolerance <rad>]\n"
        "  OR\n"
        "          --pixel-size <µm> --focal-length-mm <mm>\n"
        "\n"
        "Options:\n"
        "  --algorithm pyramid|tetra3  (default: pyramid)\n"
        "  --pattern-max-error <v>     tetra3 ratio-bin size (default: 0.01)\n"
        "  --pattern-checking-stars N  tetra3 image stars to hash (default: 10)\n"
        "  --catalog-pattern-stars N   tetra3 catalog stars to hash (default: 48)\n"
        "  --pattern-cutoff N          tetra3 max catalog patterns (default: 200000)\n"
        "\n"
        "Output: RA_deg  DEC_deg  ROLL_deg  NUM_ID  NUM_CENTROIDS\n",
        argv0);
}

static std::vector<unsigned char> readBinaryFile(const std::string &path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) {
        std::fprintf(stderr, "Error: cannot open database file: %s\n", path.c_str());
        std::exit(1);
    }
    auto sz = f.tellg();
    f.seekg(0);
    std::vector<unsigned char> buf(sz);
    f.read(reinterpret_cast<char *>(buf.data()), sz);
    return buf;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    std::string imagePath;
    std::string databasePath;
    std::string catalogPath;
    double focalLengthPx    = 0.0;
    double pixelSizeUm      = 0.0;
    double focalLengthMm    = 0.0;
    double tolerance        = 0.001;   // rad  (~0.06°)
    std::string algorithm   = "pyramid";
    int    numFalseStars    = 500;
    double maxMismatchProb  = 0.001;
    long   cutoff           = 1000000;
    double patternMaxError  = 0.01;
    int    patternCheckingStars = 10;
    int    catalogPatternStars = 48;
    long   patternCutoff    = 200000;
    int    maxMagnitude     = 600;
    int    maxStars         = 500;
    int    maxCentroids     = 200;      // keep top N brightest centroids
    int    minStarPixels    = 3;        // discard clusters smaller than this
    double sigma            = 3.0;      // centroid threshold sigma multiplier

    // ── Parse arguments ─────────────────────────────────────────────────
    for (int i = 1; i < argc; ++i) {
        auto arg = [&](const char *name) {
            return std::strcmp(argv[i], name) == 0;
        };
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Missing value for %s\n", argv[i]);
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg("--image"))               imagePath       = next();
        else if (arg("--database"))        databasePath    = next();
        else if (arg("--catalog"))         catalogPath     = next();
        else if (arg("--algorithm"))       algorithm       = next();
        else if (arg("--focal-length"))    focalLengthPx   = std::atof(next());
        else if (arg("--pixel-size"))      pixelSizeUm     = std::atof(next());
        else if (arg("--focal-length-mm")) focalLengthMm   = std::atof(next());
        else if (arg("--tolerance"))       tolerance        = std::atof(next());
        else if (arg("--num-false-stars")) numFalseStars    = std::atoi(next());
        else if (arg("--max-mismatch"))    maxMismatchProb  = std::atof(next());
        else if (arg("--cutoff"))          cutoff           = std::atol(next());
        else if (arg("--pattern-max-error")) patternMaxError = std::atof(next());
        else if (arg("--pattern-checking-stars")) patternCheckingStars = std::atoi(next());
        else if (arg("--catalog-pattern-stars")) catalogPatternStars = std::atoi(next());
        else if (arg("--pattern-cutoff")) patternCutoff = std::atol(next());
        else if (arg("--max-magnitude"))   maxMagnitude     = std::atoi(next());
        else if (arg("--max-stars"))       maxStars         = std::atoi(next());
        else if (arg("--max-centroids"))   maxCentroids     = std::atoi(next());
        else if (arg("--min-star-pixels")) minStarPixels    = std::atoi(next());
        else if (arg("--sigma"))           sigma            = std::atof(next());
        else if (arg("--help") || arg("-h")) { printUsage(argv[0]); return 0; }
        else {
            std::fprintf(stderr, "Unknown option: %s\n", argv[i]);
            printUsage(argv[0]);
            return 1;
        }
    }

    if (imagePath.empty() || databasePath.empty() || catalogPath.empty()) {
        std::fprintf(stderr, "Error: --image, --database, and --catalog are required.\n");
        printUsage(argv[0]);
        return 1;
    }

    if (algorithm != "pyramid" && algorithm != "tetra3") {
        std::fprintf(stderr, "Error: --algorithm must be 'pyramid' or 'tetra3'.\n");
        return 1;
    }

    // ── Compute focal length in pixels ──────────────────────────────────
    if (focalLengthPx <= 0.0) {
        if (pixelSizeUm > 0.0 && focalLengthMm > 0.0) {
            // Convert: focal_length_px = focal_length_mm / (pixel_size_um * 1e-3)
            focalLengthPx = (focalLengthMm * 1000.0) / pixelSizeUm;
        } else {
            std::fprintf(stderr,
                "Error: provide --focal-length <px> or both --pixel-size and --focal-length-mm.\n");
            return 1;
        }
    }

    // ── Load image as grayscale ─────────────────────────────────────────
    int imgW = 0, imgH = 0, imgChannels = 0;
    unsigned char *imgData = stbi_load(imagePath.c_str(), &imgW, &imgH,
                                       &imgChannels, 1 /* force grayscale */);
    if (!imgData) {
        std::fprintf(stderr, "Error: cannot load image: %s (%s)\n",
                     imagePath.c_str(), stbi_failure_reason());
        return 1;
    }

    std::fprintf(stderr, "Image: %dx%d (from %d channels)\n", imgW, imgH, imgChannels);

    // ── Load database ───────────────────────────────────────────────────
    auto dbBlob = readBinaryFile(databasePath);
    MultiDatabase mdb(dbBlob.data());

    const unsigned char *catalogBuf = mdb.subDatabasePointer(kCatalogMagicValue);
    if (!catalogBuf) {
        std::fprintf(stderr, "Error: catalog sub-database not found in %s\n",
                     databasePath.c_str());
        stbi_image_free(imgData);
        return 1;
    }
    DeserializeContext catDes(catalogBuf);
    Catalog catalog = deserializeCatalog(&catDes);
    std::fprintf(stderr, "Catalog: %zu stars\n", catalog.size());

    // ── Build camera model ──────────────────────────────────────────────
    Camera camera((stfloat)focalLengthPx, imgW, imgH);
    std::fprintf(stderr, "Camera: focal=%.1f px, fov=%.2f deg\n",
                 (double)camera.getFocalLength(),
                 (double)camera.fov() * 180.0 / M_PI);

    // ── Step 1: Centroid detection ──────────────────────────────────────
    CenterOfGravityCentroid centroider((stfloat)sigma);
    auto centroids = centroider.compute(imgData, (size_t)imgW, (size_t)imgH);
    stbi_image_free(imgData);

    std::fprintf(stderr, "Centroids (raw): %zu\n", centroids.size());

    // Filter: remove tiny clusters (noise / hot pixels)
    centroids.erase(
        std::remove_if(centroids.begin(), centroids.end(),
            [&](const Centroid &c) { return (int)c.numPixels < minStarPixels; }),
        centroids.end());
    std::fprintf(stderr, "Centroids (>=%d px): %zu\n", minStarPixels, centroids.size());

    // Sort by total pixel count (descending) as a proxy for brightness
    std::sort(centroids.begin(), centroids.end(),
              [](const Centroid &a, const Centroid &b) {
                  return a.numPixels > b.numPixels;
              });

    // Keep only the top N
    if ((int)centroids.size() > maxCentroids) {
        centroids.resize(maxCentroids);
        std::fprintf(stderr, "Centroids (top %d): %zu\n", maxCentroids, centroids.size());
    }

    if (centroids.empty()) {
        std::fprintf(stderr, "No centroids found — cannot proceed.\n");
        std::printf("NaN NaN NaN 0 0\n");
        return 1;
    }

    // Convert Centroids → Stars
    Stars stars;
    stars.reserve(centroids.size());
    for (auto &c : centroids)
        stars.emplace_back(c.x, c.y, c.radiusX, c.radiusY, 0);

    // ── Step 2: Star identification ─────────────────────────────────────
    StarIdentifiers ids;
    if (algorithm == "tetra3") {
        Tetra3StarIdAlgorithm tetra3((stfloat)tolerance,
                                     (stfloat)patternMaxError,
                                     patternCheckingStars,
                                     catalogPatternStars,
                                     patternCutoff);
        ids = tetra3.identify(dbBlob.data(), stars, catalog, camera);
    } else {
        PyramidStarIdAlgorithm pyramid((stfloat)tolerance, numFalseStars,
                                       (stfloat)maxMismatchProb, cutoff);
        ids = pyramid.identify(dbBlob.data(), stars, catalog, camera);
    }

    std::fprintf(stderr, "Identified: %zu / %zu stars\n", ids.size(), stars.size());

    if (ids.size() < 2) {
        std::fprintf(stderr, "Insufficient identifications for attitude estimation.\n");
        std::printf("NaN NaN NaN %zu %zu\n", ids.size(), centroids.size());
        return 1;
    }

    // ── Step 3: Attitude estimation ─────────────────────────────────────
    QuestAlgorithm quest;
    Attitude att = quest.estimate(camera, stars, catalog, ids);

    if (!att.isKnown()) {
        std::fprintf(stderr, "Attitude estimation failed.\n");
        std::printf("NaN NaN NaN %zu %zu\n", ids.size(), centroids.size());
        return 1;
    }

    EulerAngles ea = att.toSpherical();
    double raDeg   = ea.ra   * 180.0 / M_PI;
    double deDeg   = ea.de   * 180.0 / M_PI;
    double rollDeg = ea.roll * 180.0 / M_PI;

    std::fprintf(stderr, "Attitude:  RA=%.4f°  DEC=%.4f°  Roll=%.4f°\n",
                 raDeg, deDeg, rollDeg);

    // Machine-readable output on stdout
    std::printf("%.6f %.6f %.6f %zu %zu\n",
                raDeg, deDeg, rollDeg, ids.size(), centroids.size());

    return 0;
}
