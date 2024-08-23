#ifndef SAPLING_ESTD_H_
#define SAPLING_ESTD_H_

#include <algorithm>
#include <charconv>
#include <numeric>
#include <stdexcept>
#include <string>

#include "absl/strings/str_format.h"

// estd contains extensions to std that probably should have been there
namespace estd {

// Debug support (we try very hard not to use conditional compilation)
#ifndef NDEBUG
inline constexpr bool is_debug_enabled = true;
#else
inline constexpr bool is_debug_enabled = false;
#endif

}  // namespace estd

#endif // SAPLING_ESTD_H_
