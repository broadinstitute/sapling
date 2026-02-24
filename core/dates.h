#ifndef SAPLING_DATES_H_
#define SAPLING_DATES_H_

#include <string>
#include <string_view>

#include "absl/time/time.h"

namespace sapling {

// Parse an ISO date string (YYYY-MM-DD) and return the time in years relative to the given epoch.
// Throws std::invalid_argument on malformed input.
auto parse_iso_date(std::string_view iso_date_str, const absl::Time& epoch) -> double;

// Convert a time in years relative to the given epoch to an ISO date string (YYYY-MM-DD).
auto to_iso_date(double t, const absl::Time& epoch) -> std::string;

}  // namespace sapling

#endif // SAPLING_DATES_H_
