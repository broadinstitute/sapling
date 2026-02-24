#ifndef SAPLING_TIP_FILE_H_
#define SAPLING_TIP_FILE_H_

#include <istream>
#include <string>
#include <utility>
#include <vector>

#include "absl/time/time.h"

namespace sapling {

// Parse a tip file, returning a vector of (name, time) pairs.
// Each non-empty line should have the format "Name|YYYY-MM-DD".
// Times are in years relative to the epoch t0.
// Throws on malformed input or if no tips are found.
auto parse_tip_file(std::istream& is, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>>;

// Convenience overload that opens a file by name.
// Throws if the file cannot be opened.
auto parse_tip_file(const std::string& filename, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>>;

}  // namespace sapling

#endif // SAPLING_TIP_FILE_H_
