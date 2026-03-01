// Precision-adaptive tolerances for tests that must pass in both float and double modes.
// Float has ~7 significant digits (epsilon ~1.19e-7), double ~15 digits (epsilon ~2.22e-16).

#pragma once
#include "stfloat.h"
#include <type_traits>

// For checking values expected to be exactly zero (e.g. sin(0), cos(pi/2)):
static constexpr double NEAR_ZERO = std::is_same<stfloat, float>::value ? 1e-5 : 1e-15;

// General-purpose tolerance for computed results:
static constexpr double LOOSE_EPS = std::is_same<stfloat, float>::value ? 1e-4 : 1e-10;

// Tighter tolerance for simple operations (magnitudes, dot products):
static constexpr double TIGHT_EPS = std::is_same<stfloat, float>::value ? 1e-4 : 1e-12;
