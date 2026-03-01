// generate-stars.cc — Generate synthetic star field images with known attitudes.
//
// Modeled after LOST's GeneratedPipelineInput (lost/src/io.cpp).
// Each star is rendered as a 2-D Gaussian PSF onto a photon buffer,
// then converted to grayscale pixel values with optional noise.
//
// The core generation logic lives in generateStarImage(), which is a
// pure function returning a pixel buffer.  main() is a thin CLI wrapper
// that calls it in a loop and writes PNGs + metadata.
//
// Usage:
//   generate-stars --catalog <bsc.tsv> --output <dir>
//                  --count <N>
//                  [--fov <deg>]  [--width <px>]  [--height <px>]
//                  [--ra <deg>]   [--dec <deg>]   [--roll <deg>]
//                  [--random-attitudes]
//                  [--spread <stddev>]
//                  [--photons <zero-mag-photons>]
//                  [--saturation <photons>]
//                  [--dark-current <0-1>]
//                  [--read-noise <stddev>]
//                  [--shot-noise]
//                  [--max-magnitude <mag100>]
//                  [--max-stars <N>]
//                  [--seed <int>]
//
// Output:
//   For each image, writes:
//     <output>/image_NNN.png      — 8-bit grayscale PNG
//     <output>/image_NNN.txt      — metadata (true RA, DEC, Roll, star count)
//
//   Also writes a summary file <output>/manifest.txt listing all images.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

// stb_image_write for PNG output
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// star-tracker library
#include "camera.hpp"
#include "io.hpp"
#include "star.hpp"

// ── Constants ───────────────────────────────────────────────────────────────

static constexpr int    kMaxBrightness = 255;
static constexpr double kPi            = 3.14159265358979323846;

// ═══════════════════════════════════════════════════════════════════════════
//  Image generation parameters
// ═══════════════════════════════════════════════════════════════════════════

/// All tuneable parameters for synthetic star image generation.
struct StarImageParams {
    double spreadStdDev      = 2.0;       // Gaussian PSF σ (pixels)
    double zeroMagPhotons    = 200000.0;  // photons from a mag-0 star
    double saturationPhotons = 1000.0;    // photon count for pixel = 1.0
    double darkCurrent       = 0.0;       // constant background [0,1]
    double readNoiseStdDev   = 0.01;      // Gaussian read-noise σ
    bool   shotNoise         = false;     // Poisson shot noise
    double exposureTime      = 0.5;       // exposure time (arbitrary units)
};

/// Result returned by generateStarImage().
struct StarImageResult {
    std::vector<unsigned char> pixels;    // row-major 8-bit grayscale
    int width;
    int height;
    int numStarsRendered;                 // stars that fell inside the FOV
};

// ═══════════════════════════════════════════════════════════════════════════
//  generateStarImage  —  pure image-generation function
// ═══════════════════════════════════════════════════════════════════════════

/// Render a synthetic star field image for the given attitude.
///
/// @param catalog   Narrowed star catalog (positions + magnitudes).
/// @param camera    Pinhole camera model (focal length, width, height).
/// @param attitude  Orientation of the camera (body→inertial rotation).
/// @param params    Sensor / rendering parameters.
/// @param rng       Random engine (mutated for noise; pass by reference).
/// @return          The rendered 8-bit grayscale image and metadata.
StarImageResult generateStarImage(const Catalog           &catalog,
                                  const Camera            &camera,
                                  const Attitude          &attitude,
                                  const StarImageParams   &params,
                                  std::default_random_engine &rng) {

    const int imgWidth  = camera.getXResolution();
    const int imgHeight = camera.getYResolution();

    // ── Magnitude → relative brightness ─────────────────────────────────
    // From LOST: brightness = 10^(-mag/250).  Magnitude is stored as mag×100.
    auto magToBrightness = [](int mag100) -> double {
        return std::pow(10.0, -mag100 / 250.0);
    };

    // ── Gaussian PSF helper ─────────────────────────────────────────────
    auto staticPixelBrightness = [&](double px, double py,
                                     double starX, double starY,
                                     double peakBrightness) -> double {
        double dx = starX - px;
        double dy = starY - py;
        double d2 = dx * dx + dy * dy;
        return peakBrightness * params.exposureTime
             * std::exp(-d2 / (2.0 * params.spreadStdDev * params.spreadStdDev));
    };

    // ── Precompute peak photon density ──────────────────────────────────
    double zeroMagPeakDensity =
        params.zeroMagPhotons
        / (2.0 * kPi * params.spreadStdDev * params.spreadStdDev);

    // ── Project catalog stars onto sensor ────────────────────────────────
    struct RenderedStar {
        double x, y;           // sub-pixel position
        double peakBrightness; // peak photon density per time unit
        double radius;         // rendering radius in pixels
    };
    std::vector<RenderedStar> renderedStars;

    for (const auto &catStar : catalog) {
        Vec3 rotated = attitude.rotate(catStar.spatial);
        if (rotated.x <= 0) continue;

        Vec2 camCoords = camera.spatialToCamera(rotated);
        if (!camera.inSensor(camCoords)) continue;

        double brightness = zeroMagPeakDensity * magToBrightness(catStar.magnitude);

        // Rendering radius from LOST:
        //   ceil(sqrt(-log(threshold / peak / exposure) × 2π σ²))
        double interestingThreshold = 0.05;
        double logArg = interestingThreshold / (brightness * params.exposureTime);
        double radius = 1.0;
        if (logArg > 0.0 && logArg < 1.0) {
            radius = std::ceil(std::sqrt(-std::log(logArg) * 2.0 * kPi
                                         * params.spreadStdDev * params.spreadStdDev));
        }

        renderedStars.push_back({(double)camCoords.x, (double)camCoords.y,
                                 brightness, radius});
    }

    // ── Render photon buffer ────────────────────────────────────────────
    std::vector<double> photons(imgWidth * imgHeight, 0.0);

    for (const auto &star : renderedStars) {
        int xMin = std::max(0,            (int)(star.x - star.radius));
        int xMax = std::min(imgWidth - 1,  (int)(star.x + star.radius));
        int yMin = std::max(0,             (int)(star.y - star.radius));
        int yMax = std::min(imgHeight - 1, (int)(star.y + star.radius));

        for (int yp = yMin; yp <= yMax; ++yp) {
            for (int xp = xMin; xp <= xMax; ++xp) {
                double px = xp + 0.5;
                double py = yp + 0.5;
                photons[xp + yp * imgWidth] +=
                    staticPixelBrightness(px, py, star.x, star.y,
                                          star.peakBrightness);
            }
        }
    }

    // ── Convert photons → pixel brightness ──────────────────────────────
    std::vector<unsigned char> imageData(imgWidth * imgHeight);
    std::normal_distribution<double> readNoiseDist(0.0, params.readNoiseStdDev);

    for (int i = 0; i < imgWidth * imgHeight; ++i) {
        double brightness = 0.0;

        // Dark current (constant background)
        brightness += params.darkCurrent;

        // Read noise (Gaussian)
        brightness += readNoiseDist(rng);

        // Quantize photons (with optional shot noise)
        long quantized;
        if (params.shotNoise && photons[i] > 0) {
            std::poisson_distribution<long> shotDist(photons[i]);
            quantized = shotDist(rng);
        } else {
            quantized = std::lround(photons[i]);
        }

        // Convert photons to brightness via saturation level
        brightness += (double)quantized / params.saturationPhotons;

        // Clamp to [0, 1] and scale to [0, 255]
        brightness = std::max(0.0, std::min(1.0, brightness));
        imageData[i] = (unsigned char)std::floor(brightness * kMaxBrightness);
    }

    return {std::move(imageData), imgWidth, imgHeight,
            (int)renderedStars.size()};
}

// ═══════════════════════════════════════════════════════════════════════════
//  CLI helpers
// ═══════════════════════════════════════════════════════════════════════════

static void printUsage(const char *argv0) {
    std::fprintf(stderr,
        "Usage: %s --catalog <bsc.tsv> --output <dir> --count <N>\n"
        "          [--fov <deg>] [--width <px>] [--height <px>]\n"
        "          [--ra <deg>] [--dec <deg>] [--roll <deg>]\n"
        "          [--random-attitudes]\n"
        "          [--spread <stddev>]     Star spread (default: 2.0)\n"
        "          [--photons <N>]         Zero-mag total photons (default: 200000)\n"
        "          [--saturation <N>]      Saturation photons (default: 1000)\n"
        "          [--dark-current <0-1>]  Dark current level (default: 0.0)\n"
        "          [--read-noise <stddev>] Read noise stddev (default: 0.01)\n"
        "          [--shot-noise]          Enable shot noise\n"
        "          [--exposure <time>]     Exposure time (default: 0.5)\n"
        "          [--max-magnitude <M>]   Max magnitude*100 (default: 700)\n"
        "          [--max-stars <N>]       Max catalog stars (default: 5000)\n"
        "          [--seed <int>]          Random seed (default: time-based)\n",
        argv0);
}

static void mkdirp(const std::string &path) {
    struct stat st{};
    if (stat(path.c_str(), &st) == 0) return;
    mkdir(path.c_str(), 0755);
}

// ═══════════════════════════════════════════════════════════════════════════
//  main  —  CLI wrapper: parse args, call generateStarImage(), write files
// ═══════════════════════════════════════════════════════════════════════════

int main(int argc, char *argv[]) {
    // Defaults
    std::string catalogPath;
    std::string outputDir;
    int    count           = 1;
    double fovDeg          = 15.0;
    int    imgWidth        = 512;
    int    imgHeight       = 512;
    double raDeg           = 0.0;
    double decDeg          = 0.0;
    double rollDeg         = 0.0;
    bool   randomAttitudes = false;
    int    maxMagnitude    = 700;
    int    maxStars        = 5000;
    int    seed            = -1;  // -1 = time-based

    StarImageParams params;

    // ── Parse arguments ─────────────────────────────────────────────────
    for (int i = 1; i < argc; ++i) {
        auto arg = [&](const char *name) { return std::strcmp(argv[i], name) == 0; };
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) { std::fprintf(stderr, "Missing value for %s\n", argv[i]); std::exit(1); }
            return argv[++i];
        };

        if      (arg("--catalog"))           catalogPath              = next();
        else if (arg("--output"))            outputDir                = next();
        else if (arg("--count"))             count                    = std::atoi(next());
        else if (arg("--fov"))               fovDeg                   = std::atof(next());
        else if (arg("--width"))             imgWidth                 = std::atoi(next());
        else if (arg("--height"))            imgHeight                = std::atoi(next());
        else if (arg("--ra"))                raDeg                    = std::atof(next());
        else if (arg("--dec"))               decDeg                   = std::atof(next());
        else if (arg("--roll"))              rollDeg                  = std::atof(next());
        else if (arg("--random-attitudes"))  randomAttitudes          = true;
        else if (arg("--spread"))            params.spreadStdDev      = std::atof(next());
        else if (arg("--photons"))           params.zeroMagPhotons    = std::atof(next());
        else if (arg("--saturation"))        params.saturationPhotons = std::atof(next());
        else if (arg("--dark-current"))      params.darkCurrent       = std::atof(next());
        else if (arg("--read-noise"))        params.readNoiseStdDev   = std::atof(next());
        else if (arg("--shot-noise"))        params.shotNoise         = true;
        else if (arg("--exposure"))          params.exposureTime      = std::atof(next());
        else if (arg("--max-magnitude"))     maxMagnitude             = std::atoi(next());
        else if (arg("--max-stars"))         maxStars                 = std::atoi(next());
        else if (arg("--seed"))              seed                     = std::atoi(next());
        else if (arg("--help") || arg("-h")) { printUsage(argv[0]); return 0; }
        else { std::fprintf(stderr, "Unknown option: %s\n", argv[i]); printUsage(argv[0]); return 1; }
    }

    if (catalogPath.empty() || outputDir.empty()) {
        std::fprintf(stderr, "Error: --catalog and --output are required.\n");
        printUsage(argv[0]);
        return 1;
    }

    // ── Read and narrow catalog ─────────────────────────────────────────
    Catalog fullCatalog = bscParse(catalogPath);
    std::fprintf(stderr, "Full catalog: %zu stars\n", fullCatalog.size());

    Catalog catalog = narrowCatalog(fullCatalog, maxMagnitude, maxStars, 0.0);
    std::fprintf(stderr, "Narrowed catalog: %zu stars (maxMag=%d, maxStars=%d)\n",
                 catalog.size(), maxMagnitude, maxStars);

    // ── Set up camera ───────────────────────────────────────────────────
    double fovRad = fovDeg * kPi / 180.0;
    double focalLength = (double)imgWidth / (2.0 * std::tan(fovRad / 2.0));
    Camera camera((stfloat)focalLength, imgWidth, imgHeight);
    std::fprintf(stderr, "Camera: %dx%d, fov=%.1f°, fl=%.1f px\n",
                 imgWidth, imgHeight, fovDeg, focalLength);

    // ── RNG ─────────────────────────────────────────────────────────────
    if (seed < 0) seed = (int)std::time(nullptr);
    std::default_random_engine rng(seed);
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    std::fprintf(stderr, "Seed: %d\n", seed);

    // ── Create output directory ─────────────────────────────────────────
    mkdirp(outputDir);

    std::string manifestPath = outputDir + "/manifest.txt";
    std::ofstream manifest(manifestPath);
    if (!manifest) {
        std::fprintf(stderr, "Error: cannot open %s for writing\n", manifestPath.c_str());
        return 1;
    }
    manifest << "# image_file  ra_deg  dec_deg  roll_deg  num_stars  fov_deg  width  height  focal_length_px\n";

    // ── Generate images ─────────────────────────────────────────────────
    for (int imgIdx = 0; imgIdx < count; ++imgIdx) {
        // Determine attitude
        double curRa, curDec, curRoll;
        if (randomAttitudes) {
            curRa   = uniformDist(rng) * 360.0;
            curDec  = std::asin(uniformDist(rng) * 2.0 - 1.0) * 180.0 / kPi;
            curRoll = uniformDist(rng) * 360.0;
        } else {
            curRa   = raDeg;
            curDec  = decDeg;
            curRoll = rollDeg;
        }

        double raRad   = curRa   * kPi / 180.0;
        double decRad  = curDec  * kPi / 180.0;
        double rollRad = curRoll * kPi / 180.0;

        Quaternion q = sphericalToQuaternion((stfloat)raRad, (stfloat)decRad, (stfloat)rollRad);
        Attitude attitude(q);

        // ── Generate the image ──────────────────────────────────────────
        StarImageResult result =
            generateStarImage(catalog, camera, attitude, params, rng);

        std::fprintf(stderr, "[%d/%d] RA=%.2f° DEC=%.2f° Roll=%.2f°  %d stars in FOV\n",
                     imgIdx + 1, count, curRa, curDec, curRoll,
                     result.numStarsRendered);

        // ── Write PNG ───────────────────────────────────────────────────
        std::ostringstream baseName;
        baseName << "image_" << std::setw(3) << std::setfill('0') << imgIdx;

        std::string pngPath = outputDir + "/" + baseName.str() + ".png";
        if (!stbi_write_png(pngPath.c_str(), result.width, result.height, 1,
                            result.pixels.data(), result.width)) {
            std::fprintf(stderr, "  Error writing PNG: %s\n", pngPath.c_str());
            continue;
        }

        // ── Write metadata file ─────────────────────────────────────────
        std::string metaPath = outputDir + "/" + baseName.str() + ".txt";
        std::ofstream meta(metaPath);
        meta << std::fixed << std::setprecision(6);
        meta << "ra_deg " << curRa << "\n"
             << "dec_deg " << curDec << "\n"
             << "roll_deg " << curRoll << "\n"
             << "fov_deg " << fovDeg << "\n"
             << "width " << result.width << "\n"
             << "height " << result.height << "\n"
             << "focal_length_px " << focalLength << "\n"
             << "num_stars " << result.numStarsRendered << "\n"
             << "seed " << seed << "\n"
             << "spread_stddev " << params.spreadStdDev << "\n"
             << "zero_mag_photons " << params.zeroMagPhotons << "\n"
             << "saturation_photons " << params.saturationPhotons << "\n"
             << "dark_current " << params.darkCurrent << "\n"
             << "read_noise_stddev " << params.readNoiseStdDev << "\n"
             << "shot_noise " << (params.shotNoise ? "true" : "false") << "\n"
             << "exposure_time " << params.exposureTime << "\n";
        meta.close();

        // ── Append to manifest ──────────────────────────────────────────
        manifest << std::fixed << std::setprecision(6);
        manifest << baseName.str() << ".png  "
                 << curRa << "  " << curDec << "  " << curRoll << "  "
                 << result.numStarsRendered << "  "
                 << fovDeg << "  " << result.width << "  " << result.height << "  "
                 << focalLength << "\n";
    }

    manifest.close();
    std::fprintf(stderr, "Done. Wrote %d images to %s/\n", count, outputDir.c_str());
    std::fprintf(stderr, "Manifest: %s\n", manifestPath.c_str());
    return 0;
}
