#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <cmath>
#include <random>
#include <numbers>

#include "coal_sim.h"
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

TEST(Pop_model_test, const_pop_model_inverse_cum_pop_roundtrip) {
  auto pop_model = Const_pop_model{20.0};

  for (auto t : {-3.0, -1.0, 0.0, 1.0, 5.0, 100.0}) {
    EXPECT_THAT(pop_model.inverse_cum_pop(pop_model.cum_pop_at_time(t)),
                testing::DoubleNear(t, 1e-10));
  }
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

TEST(Pop_model_test, exp_pop_model_inverse_cum_pop_roundtrip) {
  auto n0 = 10.0;
  auto g = std::numbers::ln2;
  auto pop_model = Exp_pop_model{0.0, n0, g};

  for (auto t : {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0}) {
    EXPECT_THAT(pop_model.inverse_cum_pop(pop_model.cum_pop_at_time(t)),
                testing::DoubleNear(t, 1e-10));
  }
}

TEST(Pop_model_test, exp_pop_model_inverse_cum_pop_zero_growth) {
  auto n0 = 10.0;
  auto pop_model = Exp_pop_model{0.0, n0, 0.0};

  for (auto t : {-3.0, 0.0, 1.0, 5.0}) {
    EXPECT_THAT(pop_model.inverse_cum_pop(pop_model.cum_pop_at_time(t)),
                testing::DoubleNear(t, 1e-10));
  }
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

// Skygrid tests
// =============

TEST(Pop_model_test, skygrid_construction_invalid) {
  // Need at least 2 knots
  EXPECT_THROW(
      (Skygrid_pop_model{{0.0}, {1.0}, Skygrid_pop_model::Type::k_staircase}),
      std::invalid_argument);

  // Mismatched sizes
  EXPECT_THROW(
      (Skygrid_pop_model{{0.0, 1.0}, {1.0}, Skygrid_pop_model::Type::k_staircase}),
      std::invalid_argument);

  // Non-increasing knot times
  EXPECT_THROW(
      (Skygrid_pop_model{{1.0, 0.0}, {1.0, 2.0}, Skygrid_pop_model::Type::k_staircase}),
      std::invalid_argument);

  // Equal knot times
  EXPECT_THROW(
      (Skygrid_pop_model{{0.0, 0.0}, {1.0, 2.0}, Skygrid_pop_model::Type::k_staircase}),
      std::invalid_argument);
}

TEST(Pop_model_test, skygrid_staircase_pop_at_time) {
  // 3 knots at t=0, 1, 2; gamma = ln(10), ln(20), ln(30)
  auto g0 = std::log(10.0);
  auto g1 = std::log(20.0);
  auto g2 = std::log(30.0);
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0}, {g0, g1, g2}, Skygrid_pop_model::Type::k_staircase};

  // Interval 0 (t <= 0): N = exp(gamma_0) = 10
  EXPECT_THAT(pop_model.pop_at_time(-1.0), testing::DoubleNear(10.0, 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(0.0), testing::DoubleNear(10.0, 1e-10));

  // Interval 1 (0 < t <= 1): N = exp(gamma_1) = 20
  EXPECT_THAT(pop_model.pop_at_time(0.5), testing::DoubleNear(20.0, 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(20.0, 1e-10));

  // Interval 2 (1 < t <= 2): N = exp(gamma_2) = 30
  EXPECT_THAT(pop_model.pop_at_time(1.5), testing::DoubleNear(30.0, 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(2.0), testing::DoubleNear(30.0, 1e-10));

  // Interval 3 (t > 2): N = exp(gamma_2) = 30
  EXPECT_THAT(pop_model.pop_at_time(3.0), testing::DoubleNear(30.0, 1e-10));
}

TEST(Pop_model_test, skygrid_log_linear_pop_at_time) {
  // 3 knots at t=0, 1, 2; gamma = ln(10), ln(20), ln(20)
  auto g0 = std::log(10.0);
  auto g1 = std::log(20.0);
  auto g2 = std::log(20.0);
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0}, {g0, g1, g2}, Skygrid_pop_model::Type::k_log_linear};

  // Interval 0 (t <= 0): N = exp(gamma_0) = 10
  EXPECT_THAT(pop_model.pop_at_time(-1.0), testing::DoubleNear(10.0, 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(0.0), testing::DoubleNear(10.0, 1e-10));

  // Interval 1 (0 < t <= 1): log N(t) = (1-c)*g0 + c*g1, c = t
  // At t=0.5: log N = 0.5*ln(10) + 0.5*ln(20)
  EXPECT_THAT(pop_model.pop_at_time(0.5), testing::DoubleNear(std::exp(0.5*std::log(10.0) + 0.5*std::log(20.0)), 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(20.0, 1e-10));

  // Interval 2 (1 < t <= 2): gamma_1 == gamma_2, so N = 20 throughout
  EXPECT_THAT(pop_model.pop_at_time(1.5), testing::DoubleNear(20.0, 1e-10));
  EXPECT_THAT(pop_model.pop_at_time(2.0), testing::DoubleNear(20.0, 1e-10));

  // Interval 3 (t > 2): N = exp(gamma_2) = 20
  EXPECT_THAT(pop_model.pop_at_time(3.0), testing::DoubleNear(20.0, 1e-10));
}

TEST(Pop_model_test, skygrid_staircase_intensity_roundtrip) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_staircase};

  for (auto t : {-2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0}) {
    auto I = pop_model.intensity_at_time(t);
    EXPECT_THAT(pop_model.inverse_intensity(I), testing::DoubleNear(t, 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_log_linear_intensity_roundtrip) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_log_linear};

  for (auto t : {-2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0}) {
    auto I = pop_model.intensity_at_time(t);
    EXPECT_THAT(pop_model.inverse_intensity(I), testing::DoubleNear(t, 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_staircase_cum_pop_roundtrip) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_staircase};

  for (auto t : {-2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0}) {
    auto P = pop_model.cum_pop_at_time(t);
    EXPECT_THAT(pop_model.inverse_cum_pop(P), testing::DoubleNear(t, 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_log_linear_cum_pop_roundtrip) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_log_linear};

  for (auto t : {-2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0}) {
    auto P = pop_model.cum_pop_at_time(t);
    EXPECT_THAT(pop_model.inverse_cum_pop(P), testing::DoubleNear(t, 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_staircase_intensity_numerical_derivative) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_staircase};

  auto dt = 1e-6;
  // Avoid testing exactly at knot boundaries where the derivative is discontinuous
  for (auto t : {-2.0, -0.5, 0.5, 1.5, 2.5, 4.0}) {
    auto dI = pop_model.intensity_at_time(t + dt) - pop_model.intensity_at_time(t);
    auto expected = dt / pop_model.pop_at_time(t);
    EXPECT_THAT(dI, testing::DoubleNear(expected, 1e-8))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_log_linear_intensity_numerical_derivative) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_log_linear};

  auto dt = 1e-6;
  for (auto t : {-2.0, -0.5, 0.5, 1.5, 2.5, 4.0}) {
    auto dI = pop_model.intensity_at_time(t + dt) - pop_model.intensity_at_time(t);
    auto expected = dt / pop_model.pop_at_time(t);
    EXPECT_THAT(dI, testing::DoubleNear(expected, 1e-8))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_staircase_cum_pop_numerical_derivative) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_staircase};

  auto dt = 1e-6;
  for (auto t : {-2.0, -0.5, 0.5, 1.5, 2.5, 4.0}) {
    auto dP = pop_model.cum_pop_at_time(t + dt) - pop_model.cum_pop_at_time(t);
    auto expected = dt * pop_model.pop_at_time(t);
    EXPECT_THAT(dP, testing::DoubleNear(expected, 1e-8))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_log_linear_cum_pop_numerical_derivative) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0, 3.0},
      {std::log(10.0), std::log(20.0), std::log(5.0), std::log(50.0)},
      Skygrid_pop_model::Type::k_log_linear};

  auto dt = 1e-6;
  for (auto t : {-2.0, -0.5, 0.5, 1.5, 2.5, 4.0}) {
    auto dP = pop_model.cum_pop_at_time(t + dt) - pop_model.cum_pop_at_time(t);
    auto expected = dt * pop_model.pop_at_time(t);
    EXPECT_THAT(dP, testing::DoubleNear(expected, 1e-8))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_uniform_gamma_matches_const_pop) {
  // Uniform gamma should behave like a constant population model
  auto N = 42.0;
  auto g = std::log(N);
  auto const_model = Const_pop_model{N};
  auto skygrid_model = Skygrid_pop_model{
      {0.0, 1.0, 2.0}, {g, g, g}, Skygrid_pop_model::Type::k_staircase};

  for (auto t : {-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0}) {
    EXPECT_THAT(skygrid_model.pop_at_time(t),
                testing::DoubleNear(const_model.pop_at_time(t), 1e-10))
        << "at t=" << t;
  }

  // Intensity and cum_pop differ by a constant offset (different reference points),
  // but differences should match
  for (auto t : {-1.0, 0.5, 1.0, 1.5, 2.0, 3.0}) {
    auto dI_skygrid = skygrid_model.intensity_at_time(t) - skygrid_model.intensity_at_time(0.0);
    auto dI_const = const_model.intensity_at_time(t) - const_model.intensity_at_time(0.0);
    EXPECT_THAT(dI_skygrid, testing::DoubleNear(dI_const, 1e-10))
        << "at t=" << t;

    auto dP_skygrid = skygrid_model.cum_pop_at_time(t) - skygrid_model.cum_pop_at_time(0.0);
    auto dP_const = const_model.cum_pop_at_time(t) - const_model.cum_pop_at_time(0.0);
    EXPECT_THAT(dP_skygrid, testing::DoubleNear(dP_const, 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_single_interval_log_linear_matches_exp_pop) {
  // Single interval log-linear Skygrid should match Exp_pop_model
  auto n0 = 10.0;
  auto g = 2.0;
  auto x0 = 0.0;
  auto x1 = 1.0;
  auto gamma_0 = std::log(n0 * std::exp(g * (x0 - x1)));  // N at x0
  auto gamma_1 = std::log(n0);                              // N at x1 (reference point)

  auto exp_model = Exp_pop_model{x1, n0, g};
  auto skygrid_model = Skygrid_pop_model{
      {x0, x1}, {gamma_0, gamma_1}, Skygrid_pop_model::Type::k_log_linear};

  // Should match within the interval
  for (auto t : {0.0, 0.25, 0.5, 0.75, 1.0}) {
    EXPECT_THAT(skygrid_model.pop_at_time(t),
                testing::DoubleNear(exp_model.pop_at_time(t), 1e-10))
        << "at t=" << t;
  }
}

TEST(Pop_model_test, skygrid_with_coal_sim) {
  auto pop_model = Skygrid_pop_model{
      {-1.0, 0.0, 1.0},
      {std::log(10.0), std::log(20.0), std::log(15.0)},
      Skygrid_pop_model::Type::k_staircase};

  auto tip_times = std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

  auto rng = std::mt19937{42};
  auto tree = coal_sim(pop_model, tip_times, rng);
  assert_phylo_tree_integrity(tree);
}

TEST(Pop_model_test, skygrid_log_linear_with_coal_sim) {
  auto pop_model = Skygrid_pop_model{
      {-1.0, 0.0, 1.0},
      {std::log(10.0), std::log(20.0), std::log(15.0)},
      Skygrid_pop_model::Type::k_log_linear};

  auto tip_times = std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

  auto rng = std::mt19937{42};
  auto tree = coal_sim(pop_model, tip_times, rng);
  assert_phylo_tree_integrity(tree);
}

TEST(Pop_model_test, skygrid_print) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0}, {3.0, 4.0}, Skygrid_pop_model::Type::k_staircase};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::HasSubstr("Skygrid_pop_model"));
  EXPECT_THAT(ss.str(), testing::HasSubstr("k_staircase"));
}

TEST(Pop_model_test, skygrid_log_linear_print) {
  auto pop_model = Skygrid_pop_model{
      {0.0, 1.0}, {3.0, 4.0}, Skygrid_pop_model::Type::k_log_linear};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::HasSubstr("k_log_linear"));
}

}  // namespace sapling
