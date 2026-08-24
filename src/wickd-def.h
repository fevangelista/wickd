#pragma once

#include <cstddef>
#include <vector>

#include <bitset>

/// Maximum number of orbital spaces supported by graph-matrix storage.
inline constexpr std::size_t max_orbital_spaces = 8;

#define DEBUG_PRINT 0

#if DEBUG_PRINT == 1
#define WPRINT(code)                                                           \
  { code }

#else
#define WPRINT(code)                                                           \
  {}
#endif

/// Rational numbers
#include "helpers/rational.h"
using scalar_t = rational;
