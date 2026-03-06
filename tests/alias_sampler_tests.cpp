#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cmath>
#include <random>

#include "alias_sampler.h"

namespace sapling {

TEST(Alias_sampler_test, single_element) {
  auto sampler = Alias_sampler{{5.0}};
  EXPECT_EQ(sampler.size(), 1);

  auto rng = std::mt19937{42};
  for (auto i = 0; i < 100; ++i) {
    EXPECT_EQ(sampler.sample(rng), 0);
  }
}

TEST(Alias_sampler_test, single_nonzero_weight) {
  auto sampler = Alias_sampler{{0.0, 0.0, 7.0, 0.0}};

  auto rng = std::mt19937{42};
  for (auto i = 0; i < 100; ++i) {
    EXPECT_EQ(sampler.sample(rng), 2);
  }
}

TEST(Alias_sampler_test, uniform_weights) {
  auto n = 5;
  auto weights = std::vector<double>(n, 1.0);
  auto sampler = Alias_sampler{weights};
  EXPECT_EQ(sampler.size(), n);

  auto rng = std::mt19937{42};
  auto counts = std::vector<int>(n, 0);
  auto num_samples = 50000;
  for (auto i = 0; i < num_samples; ++i) {
    auto s = sampler.sample(rng);
    ASSERT_GE(s, 0);
    ASSERT_LT(s, n);
    ++counts[s];
  }

  auto p = 1.0 / n;
  auto expected = num_samples * p;
  // 3 standard deviations of a Binomial(num_samples, p) distribution
  auto tolerance = 3.0 * std::sqrt(num_samples * p * (1.0 - p));
  for (auto i = 0; i < n; ++i) {
    EXPECT_NEAR(counts[i], expected, tolerance);
  }
}

TEST(Alias_sampler_test, known_weights) {
  // Weights 1:2:3:4 => probabilities 0.1, 0.2, 0.3, 0.4
  auto weights = std::vector<double>{1.0, 2.0, 3.0, 4.0};
  auto sampler = Alias_sampler{weights};

  auto rng = std::mt19937{42};
  auto counts = std::vector<int>(4, 0);
  auto num_samples = 100000;
  for (auto i = 0; i < num_samples; ++i) {
    ++counts[sampler.sample(rng)];
  }

  auto total_weight = 10.0;
  for (auto i = 0; i < 4; ++i) {
    auto p = weights[i] / total_weight;
    auto expected = num_samples * p;
    // 3 standard deviations of a Binomial(num_samples, p) distribution
    auto tolerance = 3.0 * std::sqrt(num_samples * p * (1.0 - p));
    EXPECT_NEAR(counts[i], expected, tolerance);
  }
}

TEST(Alias_sampler_test, extreme_ratio) {
  // One very large weight and many small ones
  auto weights = std::vector<double>{1000.0, 1.0, 1.0, 1.0};
  auto sampler = Alias_sampler{weights};

  auto rng = std::mt19937{42};
  auto counts = std::vector<int>(4, 0);
  auto num_samples = 100000;
  for (auto i = 0; i < num_samples; ++i) {
    ++counts[sampler.sample(rng)];
  }

  auto total_weight = 1003.0;
  for (auto i = 0; i < 4; ++i) {
    auto p = weights[i] / total_weight;
    auto expected = num_samples * p;
    // 3 standard deviations of a Binomial(num_samples, p) distribution
    auto tolerance = 3.0 * std::sqrt(num_samples * p * (1.0 - p));
    EXPECT_NEAR(counts[i], expected, tolerance);
  }
}

}  // namespace sapling
