#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <cmath>
#include <random>
#include <numbers>

#include "pop_model.h"

namespace sapling {

TEST(Pop_model_test, const_pop_model_invalid) {
  // Don't support nonpositive populations
  EXPECT_THROW((Const_pop_model{0.0}), std::invalid_argument);
  EXPECT_THROW((Const_pop_model{-1.0}), std::invalid_argument);
}

TEST(Pop_model_test, const_pop_model_normal) {
  auto pop = 20.0;   // N_e * rho, in days

  auto pop_model = Const_pop_model{pop};  // 20.0 = N_e * rho, in days

  EXPECT_EQ(pop_model.pop(), pop);
  EXPECT_EQ(pop_model.pop_at_time(0.0), pop);
  EXPECT_EQ(pop_model.pop_at_time(5.0), pop);
  EXPECT_EQ(pop_model.intensity_at_time(0.0), 0.0 / pop);
  EXPECT_EQ(pop_model.intensity_at_time(5.0), 5.0 / pop);
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time(4.2)), testing::DoubleNear(4.2, 1e-6));
}

TEST(Pop_model_test, const_pop_model_print) {
  auto pop = 20.0;
  auto pop_model = Const_pop_model{pop};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat("Const_pop_model{pop=%g}", pop)));
}

TEST(Pop_model_test, exp_pop_model_invalid) {
  // Don't support nonpositive populations
  EXPECT_THROW((Exp_pop_model{0.0, 0.0, 1.0}), std::invalid_argument);
  EXPECT_THROW((Exp_pop_model{0.0, -1.0, 1.0}), std::invalid_argument);
}

TEST(Pop_model_test, exp_pop_model_normal) {
  auto n0 = 10.0;                // N_e(t=0) * rho, in days
  auto g = std::numbers::ln2;    // e^(ln(2)*t) = 2^t, i.e., doubles every day

  auto pop_model = Exp_pop_model{0.0, n0, g};

  EXPECT_EQ(pop_model.pop_at_t0(), n0);
  EXPECT_EQ(pop_model.growth_rate(), g);

  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(20.0, 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(-2.0), testing::DoubleNear(2.5, 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(-3.0), testing::DoubleNear(1.25, 1e-6));

  EXPECT_THAT(pop_model.cum_pop_at_time(0.0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.cum_pop_at_time(1.0), testing::DoubleNear((1/g) * (20.0 - 10.0), 1e-6));
  EXPECT_THAT(pop_model.cum_pop_at_time(-2.0), testing::DoubleNear((1/g) * (2.5 - 10.0), 1e-6));
  EXPECT_THAT(pop_model.cum_pop_at_time(-3.0), testing::DoubleNear(
      (1/g) * (1.25 - 10.0),
      1e-6));

  EXPECT_THAT(pop_model.intensity_at_time(0.0), testing::DoubleNear(
      0.0, 1e-6));
  EXPECT_THAT(pop_model.intensity_at_time(1.0), testing::DoubleNear(
      (1/n0) * (1/g) * (1.0 - 0.5), 1e-6));
  EXPECT_THAT(pop_model.intensity_at_time(2.0), testing::DoubleNear(
      (1/n0) * (1/g) * (1.0 - 0.25), 1e-6));
  EXPECT_THAT(pop_model.intensity_at_time(-2.0), testing::DoubleNear(
      (1/n0) * (1/g) * (1.0 - 4.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_at_time(-3.0), testing::DoubleNear(
      (1/n0) * (1/g) * (1.0 - 8.0), 1e-6));
  
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time( 0.0)), testing::DoubleNear( 0.0, 1e-6));
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time( 1.0)), testing::DoubleNear( 1.0, 1e-6));
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time( 2.0)), testing::DoubleNear( 2.0, 1e-6));
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time(-2.0)), testing::DoubleNear(-2.0, 1e-6));
  EXPECT_THAT(pop_model.inverse_intensity(pop_model.intensity_at_time(-3.0)), testing::DoubleNear(-3.0, 1e-6));
}

TEST(Pop_model_test, exp_pop_model_print) {
  auto t0 = 3.25;
  auto n0 = 10.0;
  auto g = 2.0;

  auto pop_model = Exp_pop_model{t0, n0, g};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat(
      "Exp_pop_model{t0=%g, n0=%g, g=%g}", t0, n0, g)));
}

}  // namespace sapling
