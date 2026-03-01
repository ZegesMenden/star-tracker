/// generate-db — Build a serialized star database from the Bright Star Catalog.
///
/// Usage:
///   generate-db [options]
///
/// Options:
///   --catalog PATH           Path to bright-star-catalog.tsv  (required)
///   --output  PATH           Output database file             (default: stdout)
///   --format  binary|header  Output format                    (default: binary)
///   --max-stars N            Keep at most N brightest stars    (default: 5000)
///   --max-magnitude M        Max magnitude × 100              (default: 600)
///   --min-separation DEG     Remove pairs closer than DEG°    (default: 0.08)
///   --kvector-min-distance D Minimum pair distance in degrees (default: 0.5)
///   --kvector-max-distance D Maximum pair distance in degrees (default: 15.0)
///   --kvector-distance-bins N Number of K-Vector bins         (default: 10000)
///   --guard NAME             Header include-guard name        (default: STARTRACKER_DATABASE_DATA_HPP)
///   --namespace NAME         C++ namespace for header output  (default: startracker)
///   --variable NAME          Array variable name              (default: kDatabase)
///   -h, --help               Print help and exit

#include "database.hpp"
#include "io.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

static constexpr stfloat kDegToRad = (stfloat)(M_PI / 180.0);

// ── Option defaults ─────────────────────────────────────────────────────────

enum class Format { Binary, Header };

struct Options {
    std::string catalogPath;
    std::string outputPath;                           // empty → stdout
    Format      format             = Format::Binary;
    int         maxStars           = 5000;
    int         maxMagnitude       = 600;             // = mag 6.00
    stfloat     minSeparation      = (stfloat)0.08;   // degrees
    stfloat     kvectorMinDistance  = (stfloat)0.5;    // degrees
    stfloat     kvectorMaxDistance  = (stfloat)15.0;   // degrees
    long        kvectorDistanceBins = 10000;

    // Header-mode options
    std::string guard     = "STARTRACKER_DATABASE_DATA_HPP";
    std::string ns        = "startracker";
    std::string variable  = "kDatabase";
};

static void printUsage(const char *prog) {
    fprintf(stderr,
        "Usage: %s --catalog PATH [options]\n"
        "\n"
        "Options:\n"
        "  --catalog PATH            Path to bright-star-catalog.tsv  (required)\n"
        "  --output  PATH            Output database file             (default: stdout)\n"
        "  --format  binary|header   Output format                    (default: binary)\n"
        "  --max-stars N             Keep at most N brightest stars    (default: 5000)\n"
        "  --max-magnitude M         Max magnitude × 100              (default: 600)\n"
        "  --min-separation DEG      Remove pairs closer than DEG°    (default: 0.08)\n"
        "  --kvector-min-distance D  Minimum pair distance in degrees (default: 0.5)\n"
        "  --kvector-max-distance D  Maximum pair distance in degrees (default: 15.0)\n"
        "  --kvector-distance-bins N Number of K-Vector bins          (default: 10000)\n"
        "\n"
        "Header-mode options (only with --format header):\n"
        "  --guard NAME              Include-guard macro name         (default: STARTRACKER_DATABASE_DATA_HPP)\n"
        "  --namespace NAME          C++ namespace                    (default: startracker)\n"
        "  --variable NAME           Array variable name              (default: kDatabase)\n"
        "  -h, --help                Print this message\n",
        prog);
}

static bool parseArgs(int argc, char *argv[], Options &opts) {
    for (int i = 1; i < argc; i++) {
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: %s requires an argument\n", argv[i]);
                return nullptr;
            }
            return argv[++i];
        };

        if (strcmp(argv[i], "--catalog") == 0) {
            const char *v = next(); if (!v) return false;
            opts.catalogPath = v;
        } else if (strcmp(argv[i], "--output") == 0) {
            const char *v = next(); if (!v) return false;
            opts.outputPath = v;
        } else if (strcmp(argv[i], "--format") == 0) {
            const char *v = next(); if (!v) return false;
            if (strcmp(v, "binary") == 0)      opts.format = Format::Binary;
            else if (strcmp(v, "header") == 0)  opts.format = Format::Header;
            else {
                fprintf(stderr, "error: --format must be 'binary' or 'header'\n");
                return false;
            }
        } else if (strcmp(argv[i], "--guard") == 0) {
            const char *v = next(); if (!v) return false;
            opts.guard = v;
        } else if (strcmp(argv[i], "--namespace") == 0) {
            const char *v = next(); if (!v) return false;
            opts.ns = v;
        } else if (strcmp(argv[i], "--variable") == 0) {
            const char *v = next(); if (!v) return false;
            opts.variable = v;
        } else if (strcmp(argv[i], "--max-stars") == 0) {
            const char *v = next(); if (!v) return false;
            opts.maxStars = atoi(v);
        } else if (strcmp(argv[i], "--max-magnitude") == 0) {
            const char *v = next(); if (!v) return false;
            opts.maxMagnitude = atoi(v);
        } else if (strcmp(argv[i], "--min-separation") == 0) {
            const char *v = next(); if (!v) return false;
            opts.minSeparation = (stfloat)atof(v);
        } else if (strcmp(argv[i], "--kvector-min-distance") == 0) {
            const char *v = next(); if (!v) return false;
            opts.kvectorMinDistance = (stfloat)atof(v);
        } else if (strcmp(argv[i], "--kvector-max-distance") == 0) {
            const char *v = next(); if (!v) return false;
            opts.kvectorMaxDistance = (stfloat)atof(v);
        } else if (strcmp(argv[i], "--kvector-distance-bins") == 0) {
            const char *v = next(); if (!v) return false;
            opts.kvectorDistanceBins = atol(v);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            exit(0);
        } else {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return false;
        }
    }
    return true;
}

// ── Main ────────────────────────────────────────────────────────────────────

int main(int argc, char *argv[]) {
    printf("stfloat size: %zu\n", sizeof(stfloat));
    Options opts;
    if (!parseArgs(argc, argv, opts)) {
        printUsage(argv[0]);
        return 1;
    }
    if (opts.catalogPath.empty()) {
        fprintf(stderr, "error: --catalog is required\n");
        printUsage(argv[0]);
        return 1;
    }

    // ── 1.  Read the Bright Star Catalog ────────────────────────────────────

    fprintf(stderr, "Reading catalog: %s\n", opts.catalogPath.c_str());
    Catalog fullCatalog = bscParse(opts.catalogPath);
    if (fullCatalog.empty()) {
        fprintf(stderr, "error: catalog is empty or unreadable\n");
        return 1;
    }
    fprintf(stderr, "  Full catalog: %zu stars\n", fullCatalog.size());

    // ── 2.  Narrow the catalog ──────────────────────────────────────────────

    Catalog catalog = narrowCatalog(fullCatalog,
                                    opts.maxMagnitude,
                                    opts.maxStars,
                                    opts.minSeparation * kDegToRad);
    fprintf(stderr, "  Narrowed catalog: %zu stars  "
            "(maxMag=%d, maxStars=%d, minSep=%.2f°)\n",
            catalog.size(), opts.maxMagnitude, opts.maxStars,
            (double)opts.minSeparation);

    // ── 3.  Build sub-databases ─────────────────────────────────────────────

    MultiDatabaseDescriptor dbEntries;

    //  3a.  Serialized catalog
    {
        SerializeContext ser;
        serializeCatalog(&ser, catalog);
        dbEntries.emplace_back(kCatalogMagicValue, std::move(ser.buffer));
        fprintf(stderr, "  Catalog sub-database: %zu bytes\n",
                dbEntries.back().bytes.size());
    }

    //  3b.  PairDistance K-Vector database
    {
        stfloat minDist = opts.kvectorMinDistance * kDegToRad;
        stfloat maxDist = opts.kvectorMaxDistance * kDegToRad;
        long    numBins = opts.kvectorDistanceBins;

        fprintf(stderr, "  Building K-Vector  (%.2f°–%.2f°, %ld bins) ...\n",
                (double)opts.kvectorMinDistance, (double)opts.kvectorMaxDistance,
                numBins);

        SerializeContext ser;
        serializePairDistanceKVector(&ser, catalog, minDist, maxDist, numBins);
        dbEntries.emplace_back(PairDistanceKVectorDatabase::kMagicValue,
                               std::move(ser.buffer));
        fprintf(stderr, "  K-Vector sub-database: %zu bytes\n",
                dbEntries.back().bytes.size());
    }

    // ── 4.  Combine into a single MultiDatabase blob ────────────────────────

    SerializeContext multi;
    serializeMultiDatabase(&multi, dbEntries);
    fprintf(stderr, "  Total database size: %zu bytes\n", multi.buffer.size());

    // ── 5.  Write output ────────────────────────────────────────────────────

    if (opts.format == Format::Binary) {
        // Raw binary output
        if (opts.outputPath.empty() || opts.outputPath == "-") {
            fwrite(multi.buffer.data(), 1, multi.buffer.size(), stdout);
        } else {
            std::ofstream out(opts.outputPath, std::ios::binary);
            if (!out) {
                fprintf(stderr, "error: could not open output file: %s\n",
                        opts.outputPath.c_str());
                return 1;
            }
            out.write(reinterpret_cast<const char *>(multi.buffer.data()),
                      (std::streamsize)multi.buffer.size());
            fprintf(stderr, "  Wrote binary: %s\n", opts.outputPath.c_str());
        }
    } else {
        // C++ header output
        FILE *out = stdout;
        if (!opts.outputPath.empty() && opts.outputPath != "-") {
            out = fopen(opts.outputPath.c_str(), "w");
            if (!out) {
                fprintf(stderr, "error: could not open output file: %s\n",
                        opts.outputPath.c_str());
                return 1;
            }
        }

        const auto &buf = multi.buffer;

        fprintf(out, "// Auto-generated by generate-db — do not edit.\n");
        fprintf(out, "//\n");
        fprintf(out, "// Stars: %zu   Database size: %zu bytes\n",
                catalog.size(), buf.size());
        fprintf(out, "//\n");
        fprintf(out, "// Options:\n");
        fprintf(out, "//   max-stars:             %d\n", opts.maxStars);
        fprintf(out, "//   max-magnitude:         %d  (%.2f)\n",
                opts.maxMagnitude, opts.maxMagnitude / 100.0);
        fprintf(out, "//   min-separation:        %.4f deg\n",
                (double)opts.minSeparation);
        fprintf(out, "//   kvector-min-distance:  %.4f deg\n",
                (double)opts.kvectorMinDistance);
        fprintf(out, "//   kvector-max-distance:  %.4f deg\n",
                (double)opts.kvectorMaxDistance);
        fprintf(out, "//   kvector-distance-bins: %ld\n",
                opts.kvectorDistanceBins);
        fprintf(out, "\n");
        fprintf(out, "#ifndef %s\n", opts.guard.c_str());
        fprintf(out, "#define %s\n\n", opts.guard.c_str());
        fprintf(out, "#include <cstddef>\n");
        fprintf(out, "#include <cstdint>\n\n");

        if (!opts.ns.empty())
            fprintf(out, "namespace %s {\n\n", opts.ns.c_str());

        fprintf(out, "/// Number of catalog stars in the database.\n");
        fprintf(out, "inline constexpr int %sNumStars = %zu;\n\n",
                opts.variable.c_str(), catalog.size());

        fprintf(out, "/// Total size of the serialized database in bytes.\n");
        fprintf(out, "inline constexpr size_t %sSize = %zuU;\n\n",
                opts.variable.c_str(), buf.size());

        fprintf(out, "/// Serialized MultiDatabase blob.\n");
        fprintf(out, "/// Use: MultiDatabase db(%s);\n",
                opts.variable.c_str());
        fprintf(out, "alignas(8) inline constexpr unsigned char %s[] = {\n",
                opts.variable.c_str());

        for (size_t i = 0; i < buf.size(); i++) {
            if (i % 16 == 0) fprintf(out, "    ");
            fprintf(out, "0x%02X,", buf[i]);
            if (i % 16 == 15 || i == buf.size() - 1)
                fprintf(out, "\n");
            else
                fprintf(out, " ");
        }

        fprintf(out, "};\n");

        if (!opts.ns.empty())
            fprintf(out, "\n} // namespace %s\n", opts.ns.c_str());

        fprintf(out, "\n#endif // %s\n", opts.guard.c_str());

        if (out != stdout) {
            fclose(out);
            fprintf(stderr, "  Wrote header: %s\n", opts.outputPath.c_str());
        }
    }

    return 0;
}
