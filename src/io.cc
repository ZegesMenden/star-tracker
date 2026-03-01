#include "io.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
//  BSC catalog parsing
// ═══════════════════════════════════════════════════════════════════════════

static constexpr stfloat kDegToRad = (stfloat)(M_PI / 180.0);

Catalog bscParse(const std::string &path) {
    Catalog result;

    FILE *file = fopen(path.c_str(), "r");
    if (!file) {
        fprintf(stderr, "error: could not open catalog file: %s\n", path.c_str());
        return result;
    }

    stfloat raj2000, dej2000;
    int magnitudeHigh, magnitudeLow;
    int name;
    char weird;

#ifdef STARTRACK_USE_FLOAT
    const char *format = "%f|%f|%d|%c|%d.%d";
#else
    const char *format = "%lf|%lf|%d|%c|%d.%d";
#endif

    while (EOF != fscanf(file, format,
                         &raj2000, &dej2000,
                         &name, &weird,
                         &magnitudeHigh, &magnitudeLow)) {
        int magnitude = magnitudeHigh * 100
                      + (magnitudeHigh < 0 ? -magnitudeLow : magnitudeLow);
        result.push_back(CatalogStar(raj2000 * kDegToRad,
                                     dej2000 * kDegToRad,
                                     magnitude,
                                     (int16_t)name));
    }

    fclose(file);
    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Catalog narrowing (filter by magnitude, separation, max count)
// ═══════════════════════════════════════════════════════════════════════════

Catalog narrowCatalog(const Catalog &catalog,
                      int maxMagnitude,
                      int maxStars,
                      stfloat minSeparation) {
    // 1.  Filter by magnitude.
    Catalog result;
    for (const auto &s : catalog) {
        if (s.magnitude <= maxMagnitude) {
            result.push_back(s);
        }
    }

    // 2.  Remove stars that are too close to another star.
    if (minSeparation > 0) {
        std::vector<bool> tooClose(result.size(), false);
        for (size_t i = 0; i < result.size(); i++) {
            if (tooClose[i]) continue;
            for (size_t j = i + 1; j < result.size(); j++) {
                if (tooClose[j]) continue;
                stfloat d = angleUnit(result[i].spatial, result[j].spatial);
                if (d < minSeparation) {
                    tooClose[i] = true;
                    tooClose[j] = true;
                }
            }
        }
        Catalog filtered;
        for (size_t i = 0; i < result.size(); i++) {
            if (!tooClose[i]) filtered.push_back(result[i]);
        }
        result = std::move(filtered);
    }

    // 3.  Sort by brightness (lower magnitude = brighter) and truncate.
    std::sort(result.begin(), result.end(),
              [](const CatalogStar &a, const CatalogStar &b) {
                  return a.magnitude < b.magnitude;
              });

    if (maxStars > 0 && (int)result.size() > maxStars) {
        result.resize(maxStars);
    }

    return result;
}
