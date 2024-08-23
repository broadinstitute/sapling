#include "pop_model.h"

#include <algorithm>
#include <cmath>

#include "absl/strings/str_format.h"

#include "distributions.h"

namespace sapling {

Const_pop_model::Const_pop_model(double pop) : pop_{pop} {
  if (pop <= 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Population size should be positive (not %f)", pop));
  }
}

auto Const_pop_model::pop_at_time(double) const -> double { return pop_; }
auto Const_pop_model::cum_pop_at_time(double t) const -> double { return t * pop_; }
auto Const_pop_model::intensity_at_time(double t) const -> double { return t / pop_; }
auto Const_pop_model::inverse_intensity(double I) const -> double { return I * pop_; }
auto Const_pop_model::sample(double min_t, double max_t, absl::BitGenRef prng) const -> double {
  return absl::Uniform(absl::IntervalClosedOpen, prng, min_t, max_t);
}

Exp_pop_model::Exp_pop_model(double t0, double pop_at_t0, double growth_rate)
    : t0_{t0}, pop_at_t0_{pop_at_t0}, growth_rate_{growth_rate} {
  if (pop_at_t0 <= 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Initial population size should be positive (not %f)", pop_at_t0));
  }
}

auto Exp_pop_model::pop_at_time(double t) const -> double {
  return pop_at_t0_ * std::exp((t - t0_) * growth_rate_);
}

auto Exp_pop_model::cum_pop_at_time(double t) const -> double {
  // int_t0^t dt' n0 exp[g (t' - t0)] = n0/g [exp[g (t - t0)] - 1]
  if (growth_rate_ == 0.0) {
    return t * pop_at_t0_;
  } else {
    return (std::exp(growth_rate_ * (t - t0_)) - 1) * pop_at_t0_ / growth_rate_;
  }
}

auto Exp_pop_model::intensity_at_time(double t) const -> double {
  // int_t0^t dt' 1/n0 exp[-g (t - t0))] = 1/(g n0) [1 - exp[-g (t-t0)]]
  if (growth_rate_ == 0.0) {
    return t / pop_at_t0_;
  } else {
    return (1 - std::exp(-growth_rate_ * (t-t0_))) / (growth_rate_ * pop_at_t0_);
  }
}

auto Exp_pop_model::inverse_intensity(double I) const -> double {
  //    I = 1/(g n0) [1 - exp[-g (t-t0)]]
  // => t = t0 - ln(1 - I*(g n0))/g
  if (growth_rate_ == 0.0) {
    return I * pop_at_t0_;
  } else {
    return t0_ - std::log1p(-I * growth_rate_ * pop_at_t0_) / growth_rate_;
  }
}

auto Exp_pop_model::sample(double min_t, double max_t, absl::BitGenRef prng) const -> double {
  return Bounded_exponential_distribution{growth_rate_, min_t, max_t}(prng);
}


}  // namespace sapling
