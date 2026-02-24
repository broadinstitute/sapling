#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <sstream>

#include "tip_file.h"

namespace sapling {

// Helper: parse a fixed epoch for all tests
static auto test_epoch() -> absl::Time {
  auto t0 = absl::Time{};
  auto err = std::string{};
  auto ok = absl::ParseTime("%Y-%m-%d", "2024-07-31", &t0, &err);
  if (not ok) { throw std::runtime_error("test_epoch: " + err); }
  return t0;
}

TEST(Tip_file_test, valid_multiple_tips) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"Alpha|2024-01-15\nBeta|2024-03-22\nGamma|2024-07-01\n"};

  auto tips = parse_tip_file(is, t0);

  ASSERT_EQ(std::ssize(tips), 3);
  EXPECT_EQ(tips[0].first, "Alpha|2024-01-15");
  EXPECT_EQ(tips[1].first, "Beta|2024-03-22");
  EXPECT_EQ(tips[2].first, "Gamma|2024-07-01");
  // All dates are before the epoch, so times should be negative
  EXPECT_LT(tips[0].second, tips[1].second);
  EXPECT_LT(tips[1].second, tips[2].second);
  EXPECT_LT(tips[2].second, 0.0);
}

TEST(Tip_file_test, single_tip) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"OnlyTip|2024-07-31\n"};

  auto tips = parse_tip_file(is, t0);

  ASSERT_EQ(std::ssize(tips), 1);
  EXPECT_EQ(tips[0].first, "OnlyTip|2024-07-31");
  EXPECT_THAT(tips[0].second, testing::DoubleNear(0.0, 1e-6));
}

TEST(Tip_file_test, empty_lines_skipped) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"\nAlpha|2024-01-15\n\nBeta|2024-03-22\n\n"};

  auto tips = parse_tip_file(is, t0);

  ASSERT_EQ(std::ssize(tips), 2);
  EXPECT_EQ(tips[0].first, "Alpha|2024-01-15");
  EXPECT_EQ(tips[1].first, "Beta|2024-03-22");
}

TEST(Tip_file_test, missing_pipe) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"NoPipeHere\n"};

  EXPECT_THROW(parse_tip_file(is, t0), std::invalid_argument);
}

TEST(Tip_file_test, bad_date) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"Foo|not-a-date\n"};

  EXPECT_THROW(parse_tip_file(is, t0), std::invalid_argument);
}

TEST(Tip_file_test, empty_file) {
  auto t0 = test_epoch();
  auto is = std::istringstream{""};

  EXPECT_THROW(parse_tip_file(is, t0), std::invalid_argument);
}

TEST(Tip_file_test, only_blank_lines) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"\n\n\n"};

  EXPECT_THROW(parse_tip_file(is, t0), std::invalid_argument);
}

TEST(Tip_file_test, name_with_multiple_pipes) {
  auto t0 = test_epoch();
  // Split on the *last* pipe to extract the date; full line is preserved as the name
  auto is = std::istringstream{"A|B|2024-07-31\n"};

  auto tips = parse_tip_file(is, t0);

  ASSERT_EQ(std::ssize(tips), 1);
  EXPECT_EQ(tips[0].first, "A|B|2024-07-31");
  EXPECT_THAT(tips[0].second, testing::DoubleNear(0.0, 1e-6));
}

TEST(Tip_file_test, name_with_spaces) {
  auto t0 = test_epoch();
  auto is = std::istringstream{"Sample One|2024-07-31\n"};

  auto tips = parse_tip_file(is, t0);

  ASSERT_EQ(std::ssize(tips), 1);
  EXPECT_EQ(tips[0].first, "Sample One|2024-07-31");
}

}  // namespace sapling
