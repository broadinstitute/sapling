#include "tip_file.h"

#include <fstream>
#include <stdexcept>

#include "absl/strings/str_format.h"

#include "dates.h"

namespace sapling {

auto parse_tip_file(std::istream& is, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>> {
  auto result = std::vector<std::pair<std::string, double>>{};

  auto line = std::string{};
  auto line_num = 0;
  while (std::getline(is, line)) {
    ++line_num;

    // Skip empty lines
    if (line.empty()) { continue; }

    // Split on the last '|'
    auto pos = line.rfind('|');
    if (pos == std::string::npos) {
      throw std::invalid_argument(absl::StrFormat(
          "Line %d: expected 'Name|YYYY-MM-DD', but no '|' found: %s", line_num, line));
    }

    auto name = line.substr(0, pos);
    auto date_str = std::string_view{line}.substr(pos + 1);

    auto t = parse_iso_date(date_str, t0);  // throws on bad date
    result.emplace_back(std::move(name), t);
  }

  if (result.empty()) {
    throw std::invalid_argument("Tip file is empty: no tips found");
  }

  return result;
}

auto parse_tip_file(const std::string& filename, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>> {
  auto is = std::ifstream{filename};
  if (not is) {
    throw std::invalid_argument(absl::StrFormat(
        "Could not open tip file '%s'", filename));
  }
  return parse_tip_file(is, t0);
}

}  // namespace sapling
