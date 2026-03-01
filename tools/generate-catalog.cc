/// generate-catalog — Emit a narrowed star catalog as a C++ header file.
///
/// This allows an embedded target to carry the catalog in flash/ROM without
/// needing the BSC TSV file or any file-system access at runtime.
///
/// Usage:
///   generate-catalog --catalog <bsc.tsv> [options]
///
/// Options:
///   --catalog PATH           Path to bright-star-catalog.tsv  (required)
///   --output  PATH           Output header file               (default: stdout)
///   --max-stars N            Keep at most N brightest stars    (default: 5000)
///   --max-magnitude M        Max magnitude × 100              (default: 600)
///   --min-separation DEG     Remove pairs closer than DEG°    (default: 0.08)
///   --guard NAME             Header include-guard macro        (default: STARTRACKER_CATALOG_DATA_HPP)
///   --namespace NAME         C++ namespace                     (default: startracker)
///   --variable NAME          Array variable name               (default: kCatalog)
///   -h, --help               Print help and exit
///
/// Output:
///   A self-contained C++ header with:
///     - A constexpr array of {x, y, z, magnitude, name} structs
///     - A helper function makeCatalog() that returns a std::vector<CatalogStar>

#include "io.hpp"
#include "star.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

static constexpr stfloat kDegToRad = (stfloat)(M_PI / 180.0);

// ── Option defaults ─────────────────────────────────────────────────────────

struct Options {
    std::string catalogPath;
    std::string outputPath;                           // empty → stdout
    int         maxStars           = 5000;
    int         maxMagnitude       = 600;             // = mag 6.00
    stfloat     minSeparation      = (stfloat)0.08;   // degrees

    std::string guard     = "STARTRACKER_CATALOG_DATA_HPP";
    std::string ns        = "startracker";
    std::string variable  = "kCatalog";
};

static void printUsage(const char *prog) {
    fprintf(stderr,
        "Usage: %s --catalog PATH [options]\n"
        "\n"
        "Options:\n"
        "  --catalog PATH            Path to bright-star-catalog.tsv   (required)\n"
        "  --output  PATH            Output header file                (default: stdout)\n"
        "  --max-stars N             Keep at most N brightest stars     (default: 5000)\n"
        "  --max-magnitude M         Max magnitude × 100               (default: 600)\n"
        "  --min-separation DEG      Remove pairs closer than DEG°     (default: 0.08)\n"
        "\n"
        "Header options:\n"
        "  --guard NAME              Include-guard macro name           (default: STARTRACKER_CATALOG_DATA_HPP)\n"
        "  --namespace NAME          C++ namespace                      (default: startracker)\n"
        "  --variable NAME           Array variable name                (default: kCatalog)\n"
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
        } else if (strcmp(argv[i], "--max-stars") == 0) {
            const char *v = next(); if (!v) return false;
            opts.maxStars = atoi(v);
        } else if (strcmp(argv[i], "--max-magnitude") == 0) {
            const char *v = next(); if (!v) return false;
            opts.maxMagnitude = atoi(v);
        } else if (strcmp(argv[i], "--min-separation") == 0) {
            const char *v = next(); if (!v) return false;
            opts.minSeparation = (stfloat)atof(v);
        } else if (strcmp(argv[i], "--guard") == 0) {
            const char *v = next(); if (!v) return false;
            opts.guard = v;
        } else if (strcmp(argv[i], "--namespace") == 0) {
            const char *v = next(); if (!v) return false;
            opts.ns = v;
        } else if (strcmp(argv[i], "--variable") == 0) {
            const char *v = next(); if (!v) return false;
            opts.variable = v;
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

    // ── 3.  Write the header ────────────────────────────────────────────────

    FILE *out = stdout;
    if (!opts.outputPath.empty() && opts.outputPath != "-") {
        out = fopen(opts.outputPath.c_str(), "w");
        if (!out) {
            fprintf(stderr, "error: could not open output file: %s\n",
                    opts.outputPath.c_str());
            return 1;
        }
    }

    const size_t N = catalog.size();

    // ── Preamble ────────────────────────────────────────────────────────

    fprintf(out, "// Auto-generated by generate-catalog — do not edit.\n");
    fprintf(out, "//\n");
    fprintf(out, "// Stars: %zu\n", N);
    fprintf(out, "//\n");
    fprintf(out, "// Options:\n");
    fprintf(out, "//   max-stars:        %d\n", opts.maxStars);
    fprintf(out, "//   max-magnitude:    %d  (%.2f)\n",
            opts.maxMagnitude, opts.maxMagnitude / 100.0);
    fprintf(out, "//   min-separation:   %.4f deg\n",
            (double)opts.minSeparation);
    fprintf(out, "\n");
    fprintf(out, "#ifndef %s\n", opts.guard.c_str());
    fprintf(out, "#define %s\n\n", opts.guard.c_str());
    fprintf(out, "#include \"star.hpp\"\n\n");
    fprintf(out, "#include <cstddef>\n");
    fprintf(out, "#include <cstdint>\n");
    fprintf(out, "#include <vector>\n\n");

    if (!opts.ns.empty())
        fprintf(out, "namespace %s {\n\n", opts.ns.c_str());

    // ── Star count ──────────────────────────────────────────────────────

    fprintf(out, "/// Number of stars in the embedded catalog.\n");
    fprintf(out, "inline constexpr int %sNumStars = %zu;\n\n",
            opts.variable.c_str(), N);

    // ── POD entry struct (constexpr-friendly) ───────────────────────────

    fprintf(out, "/// Plain-data entry for constexpr storage.\n");
    fprintf(out, "struct %sEntry {\n", opts.variable.c_str());
    fprintf(out, "    double x, y, z;    // spatial unit vector\n");
    fprintf(out, "    int    magnitude;   // magnitude × 100\n");
    fprintf(out, "    int16_t name;       // BSC catalog name/index\n");
    fprintf(out, "};\n\n");

    // ── Data array ──────────────────────────────────────────────────────

    fprintf(out, "/// Embedded catalog data.  Each entry is a star with its\n");
    fprintf(out, "/// unit-sphere direction, magnitude, and BSC name.\n");
    fprintf(out, "inline constexpr %sEntry %sData[] = {\n",
            opts.variable.c_str(), opts.variable.c_str());

    for (size_t i = 0; i < N; i++) {
        const CatalogStar &s = catalog[i];
        fprintf(out, "    {%+.17e, %+.17e, %+.17e, %4d, %5d},\n",
                (double)s.spatial.x,
                (double)s.spatial.y,
                (double)s.spatial.z,
                s.magnitude,
                (int)s.name);
    }

    fprintf(out, "};\n\n");

    // ── Helper: makeCatalog() ───────────────────────────────────────────

    fprintf(out, "/// Build a Catalog (std::vector<CatalogStar>) from the\n");
    fprintf(out, "/// embedded constexpr data.  Call once at startup.\n");
    fprintf(out, "inline Catalog %sMake() {\n", opts.variable.c_str());
    fprintf(out, "    Catalog c;\n");
    fprintf(out, "    c.reserve(%sNumStars);\n", opts.variable.c_str());
    fprintf(out, "    for (int i = 0; i < %sNumStars; i++) {\n",
            opts.variable.c_str());
    fprintf(out, "        const auto &e = %sData[i];\n", opts.variable.c_str());
    fprintf(out, "        c.push_back(CatalogStar(\n");
    fprintf(out, "            Vec3{(stfloat)e.x, (stfloat)e.y, (stfloat)e.z},\n");
    fprintf(out, "            e.magnitude, e.name));\n");
    fprintf(out, "    }\n");
    fprintf(out, "    return c;\n");
    fprintf(out, "}\n");

    // ── Footer ──────────────────────────────────────────────────────────

    if (!opts.ns.empty())
        fprintf(out, "\n} // namespace %s\n", opts.ns.c_str());

    fprintf(out, "\n#endif // %s\n", opts.guard.c_str());

    if (out != stdout) {
        fclose(out);
        fprintf(stderr, "  Wrote header: %s (%zu stars)\n",
                opts.outputPath.c_str(), N);
    }

    return 0;
}
