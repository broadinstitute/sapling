#include "dates.h"

#include <stdexcept>

#include "absl/strings/str_format.h"
#include "absl/time/time.h"

namespace sapling {

auto parse_iso_date(std::string_view iso_date_str, const absl::Time& epoch) -> double {
  auto d = absl::Time{};
  auto err = std::string{};
  if (not absl::ParseTime("%Y-%m-%d", iso_date_str, &d, &err)) {
    throw std::invalid_argument(absl::StrFormat(
        "Badly formatted ISO date: %s (error: %s)", iso_date_str, err));
  }
  return absl::ToDoubleHours(d - epoch) / 24.0 / 365.0;
}

auto to_iso_date(double t, const absl::Time& epoch) -> std::string {
  auto abs_time = absl::Time{epoch + absl::Hours(t * 24.0 * 365.0 + 1e-5)};  // +1e-5 to avoid roundoff
  return absl::FormatTime("%Y-%m-%d", abs_time, absl::UTCTimeZone());
}

}  // namespace sapling
