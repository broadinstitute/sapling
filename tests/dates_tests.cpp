#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "dates.h"

namespace sapling {

// Helper: parse a fixed epoch for all tests
static auto test_epoch() -> absl::Time {
  auto t0 = absl::Time{};
  auto err = std::string{};
  auto ok = absl::ParseTime("%Y-%m-%d", "2020-01-01", &t0, &err);
  if (not ok) { throw std::runtime_error("test_epoch: " + err); }
  return t0;
}

TEST(Dates_test, parse_iso_date_at_epoch) {
  auto t0 = test_epoch();
  EXPECT_THAT(parse_iso_date("2020-01-01", t0), testing::DoubleNear(0.0, 1e-6));
}

TEST(Dates_test, parse_iso_date_one_year_later) {
  auto t0 = test_epoch();
  // 2021-01-01 is approximately 1.0 years from 2020-01-01
  // (366 days / 365.0 because 2020 is a leap year)
  auto t = parse_iso_date("2021-01-01", t0);
  EXPECT_THAT(t, testing::DoubleNear(366.0 / 365.0, 1e-6));
}

TEST(Dates_test, parse_iso_date_before_epoch) {
  auto t0 = test_epoch();
  auto t = parse_iso_date("2019-01-01", t0);
  EXPECT_LT(t, 0.0);
  // 2019-01-01 to 2020-01-01 is 365 days
  EXPECT_THAT(t, testing::DoubleNear(-365.0 / 365.0, 1e-6));
}

TEST(Dates_test, parse_iso_date_bad_format) {
  auto t0 = test_epoch();
  EXPECT_THROW(parse_iso_date("not-a-date", t0), std::invalid_argument);
  EXPECT_THROW(parse_iso_date("2024/01/01", t0), std::invalid_argument);
  EXPECT_THROW(parse_iso_date("", t0), std::invalid_argument);
}

TEST(Dates_test, to_iso_date_at_epoch) {
  auto t0 = test_epoch();
  EXPECT_EQ(to_iso_date(0.0, t0), "2020-01-01");
}

TEST(Dates_test, to_iso_date_known_date) {
  auto t0 = test_epoch();
  // 2020-07-01 is 182 days from 2020-01-01 (2020 is a leap year)
  auto t = 182.0 / 365.0;
  EXPECT_EQ(to_iso_date(t, t0), "2020-07-01");
}

TEST(Dates_test, roundtrip) {
  auto t0 = test_epoch();

  // Parse, then convert back -- should be stable
  auto dates = std::vector<std::string>{"2019-06-15", "2020-01-01", "2024-07-31", "2025-12-25"};
  for (const auto& date : dates) {
    auto t = parse_iso_date(date, t0);
    EXPECT_EQ(to_iso_date(t, t0), date) << "Roundtrip failed for " << date;
  }
}

}  // namespace sapling
