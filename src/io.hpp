#pragma once

#include "star.hpp"

#include <string>

// ═══════════════════════════════════════════════════════════════════════════
//  Catalog I/O
// ═══════════════════════════════════════════════════════════════════════════

/// Parse the Yale Bright Star Catalog from a pipe-delimited TSV file.
///
/// Expected format per line:
///   RA_deg|DEC_deg|NAME|FLAG|MAG_high.MAG_low
///
/// Example:  001.291250|+45.229167|   1| | 6.70
Catalog bscParse(const std::string &path);

/// Filter a catalog down to the brightest, well-separated stars.
///
/// @param catalog         The full catalog.
/// @param maxMagnitude    Stars dimmer than this (magnitude × 100) are removed.
///                        Lower values = stricter.  Use e.g. 600 for mag ≤ 6.00.
/// @param maxStars        Keep at most this many of the brightest stars.
/// @param minSeparation   Remove both stars in any pair separated by less than
///                        this angle (radians).  0 to skip.
Catalog narrowCatalog(const Catalog &catalog,
                      int maxMagnitude,
                      int maxStars,
                      stfloat minSeparation);
