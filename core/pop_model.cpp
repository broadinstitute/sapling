#include "pop_model.h"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "absl/log/check.h"
#include "absl/strings/str_format.h"

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
auto Const_pop_model::inverse_cum_pop(double P) const -> double { return P / pop_; }

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

auto Exp_pop_model::inverse_cum_pop(double P) const -> double {
  // cum_pop_at_time(t) = (exp(g*(t-t0)) - 1) * n0/g  [g != 0]
  //                    = t * n0                        [g == 0]
  // Inverse: t = t0 + log1p(P * g / n0) / g           [g != 0]
  //          t = P / n0                                [g == 0]
  if (growth_rate_ == 0.0) {
    return P / pop_at_t0_;
  } else {
    return t0_ + std::log1p(P * growth_rate_ / pop_at_t0_) / growth_rate_;
  }
}


// Skygrid_pop_model
// ==================

Skygrid_pop_model::Skygrid_pop_model(
    std::vector<double> x,
    std::vector<double> gamma,
    Type type)
    : x_{std::move(x)}, gamma_{std::move(gamma)}, type_{type} {

  if (std::ssize(x_) < 2) {
    throw std::invalid_argument("Skygrid_pop_model needs at least two knots");
  }
  if (std::ssize(x_) != std::ssize(gamma_)) {
    throw std::invalid_argument(absl::StrFormat(
        "Skygrid_pop_model needs knot times array (x) and knot values array (gamma) to "
        "have equal size, but there are %d x's and %d gamma's",
        std::ssize(x_), std::ssize(gamma_)));
  }

  // Strictly increasing knot times
  for (auto i = 0; i < std::ssize(x_) - 1; ++i) {
    if (not (x_[i] < x_[i+1])) {
      throw std::invalid_argument(absl::StrFormat(
          "Skygrid_pop_model needs strictly increasing knot times, but "
          "x[%d] = %g >= x[%d] = %g",
          i, x_[i], i+1, x_[i+1]));
    }
  }

  // Precompute cumulative intensity and cumulative pop at knot boundaries.
  // cum_intensity_[0] = 0, cum_intensity_[k] = int_{x_0}^{x_k} 1/N(t) dt
  // cum_pop_[0] = 0, cum_pop_[k] = int_{x_0}^{x_k} N(t) dt
  auto num_knots = std::ssize(x_);
  cum_intensity_.resize(num_knots, 0.0);
  cum_pop_.resize(num_knots, 0.0);

  for (auto k = 1; k < num_knots; ++k) {
    auto w = x_[k] - x_[k-1];
    CHECK_GT(w, 0.0);

    auto interval_intensity = 0.0;
    auto interval_pop = 0.0;

    switch (type_) {
      case Type::k_staircase: {
        // N(t) = exp(gamma_k) throughout interval k
        auto N_k = std::exp(gamma_[k]);
        CHECK_GT(N_k, 0.0);
        interval_intensity = w / N_k;
        interval_pop = N_k * w;
        break;
      }

      case Type::k_log_linear: {
        // N(t) = n0 * exp(g * (t - x_k))  where n0 = exp(gamma_k)
        auto n0 = std::exp(gamma_[k]);
        auto g = (gamma_[k] - gamma_[k-1]) / w;
        CHECK_GT(n0, 0.0);

        if (g == 0.0) {
          interval_intensity = w / n0;
          interval_pop = n0 * w;
        } else {
          // intensity: (1/(g*n0)) * expm1(g*w)
          interval_intensity = std::expm1(g * w) / (g * n0);
          // pop: -(n0/g) * expm1(-g*w)
          interval_pop = -(n0 / g) * std::expm1(-g * w);
        }
        break;
      }
    }

    CHECK_GT(interval_intensity, 0.0);
    CHECK_GT(interval_pop, 0.0);

    cum_intensity_[k] = cum_intensity_[k-1] + interval_intensity;
    cum_pop_[k] = cum_pop_[k-1] + interval_pop;
  }
}

auto Skygrid_pop_model::interval_containing_t(double t) const -> int {
  auto it = std::ranges::lower_bound(x_, t);  // lower_bound means *(it-1) < t <= *it
  return static_cast<int>(std::ranges::distance(x_.begin(), it));
}

auto Skygrid_pop_model::pop_at_time(double t) const -> double {
  auto k = interval_containing_t(t);
  auto m = M();

  switch (type_) {
    case Type::k_staircase:
      if      (k == 0)  { return std::exp(gamma(0)); }
      else if (k <= m)  { return std::exp(gamma(k)); }
      else              { return std::exp(gamma(m)); }

    case Type::k_log_linear:
      if      (k == 0)  { return std::exp(gamma(0)); }
      else if (k <= m)  {
        auto c = (t - x(k-1)) / (x(k) - x(k-1));
        CHECK_GE(c, 0.0);
        CHECK_LE(c, 1.0);
        return std::exp((1-c)*gamma(k-1) + c*gamma(k));
      }
      else              { return std::exp(gamma(m)); }
  }
  CHECK(false) << "unrecognized type " << static_cast<int>(type_);
}

auto Skygrid_pop_model::intensity_at_time(double t) const -> double {
  auto k = interval_containing_t(t);
  auto m = M();

  // Interval 0: t <= x_0, constant N = exp(gamma_0)
  if (k == 0) {
    return (t - x(0)) / std::exp(gamma(0));
  }

  // Interval M+1: t > x_M, constant N = exp(gamma_M)
  if (k > m) {
    return cum_intensity_[m] + (t - x(m)) / std::exp(gamma(m));
  }

  // Interior interval k (1 <= k <= M)
  CHECK_GE(k, 1);
  CHECK_LE(k, m);
  auto delta = t - x(k-1);
  CHECK_GE(delta, 0.0);

  switch (type_) {
    case Type::k_staircase: {
      return cum_intensity_[k-1] + delta / std::exp(gamma(k));
    }

    case Type::k_log_linear: {
      auto n0 = std::exp(gamma(k));
      auto w = x(k) - x(k-1);
      auto g = (gamma(k) - gamma(k-1)) / w;

      if (g == 0.0) {
        return cum_intensity_[k-1] + delta / n0;
      } else {
        // (1/(g*n0)) * exp(-g*(t - x_k)) * expm1(g*delta)
        return cum_intensity_[k-1] +
            (1.0 / (g * n0)) * std::exp(-g * (t - x(k))) * std::expm1(g * delta);
      }
    }
  }
  CHECK(false) << "unrecognized type " << static_cast<int>(type_);
}

auto Skygrid_pop_model::inverse_intensity(double I) const -> double {
  auto m = M();

  // Region at or below x_0: I <= 0
  if (I <= 0.0) {
    return x(0) + I * std::exp(gamma(0));
  }

  // Region at or above x_M: I >= cum_intensity_[M]
  if (I >= cum_intensity_[m]) {
    return x(m) + (I - cum_intensity_[m]) * std::exp(gamma(m));
  }

  // Interior: 0 < I < cum_intensity_[M]
  // Binary search to find smallest k such that cum_intensity_[k] > I
  auto it = std::upper_bound(cum_intensity_.begin(), cum_intensity_.end(), I);
  CHECK(it != cum_intensity_.begin());
  CHECK(it != cum_intensity_.end());
  auto k = static_cast<int>(std::distance(cum_intensity_.begin(), it));
  CHECK_GE(k, 1);
  CHECK_LE(k, m);

  auto remaining_I = I - cum_intensity_[k-1];
  CHECK_GE(remaining_I, 0.0);

  switch (type_) {
    case Type::k_staircase: {
      auto result = x(k-1) + remaining_I * std::exp(gamma(k));
      CHECK_GE(result, x(k-1));
      CHECK_LE(result, x(k));
      return result;
    }

    case Type::k_log_linear: {
      auto w = x(k) - x(k-1);
      auto g = (gamma(k) - gamma(k-1)) / w;
      auto N_left = std::exp(gamma(k-1));

      auto delta = 0.0;
      if (g == 0.0) {
        delta = remaining_I * N_left;
      } else {
        // delta = -log1p(-remaining_I * g * exp(gamma_{k-1})) / g
        delta = -std::log1p(-remaining_I * g * N_left) / g;
      }
      CHECK_GE(delta, -1e-10);
      CHECK_LE(delta, w + 1e-10);
      return x(k-1) + delta;
    }
  }
  CHECK(false) << "unrecognized type " << static_cast<int>(type_);
}

auto Skygrid_pop_model::cum_pop_at_time(double t) const -> double {
  auto k = interval_containing_t(t);
  auto m = M();

  // Interval 0: t <= x_0, constant N = exp(gamma_0)
  if (k == 0) {
    return std::exp(gamma(0)) * (t - x(0));
  }

  // Interval M+1: t > x_M, constant N = exp(gamma_M)
  if (k > m) {
    return cum_pop_[m] + std::exp(gamma(m)) * (t - x(m));
  }

  // Interior interval k (1 <= k <= M)
  CHECK_GE(k, 1);
  CHECK_LE(k, m);
  auto delta = t - x(k-1);
  CHECK_GE(delta, 0.0);

  switch (type_) {
    case Type::k_staircase: {
      return cum_pop_[k-1] + std::exp(gamma(k)) * delta;
    }

    case Type::k_log_linear: {
      auto w = x(k) - x(k-1);
      auto g = (gamma(k) - gamma(k-1)) / w;
      auto N_left = std::exp(gamma(k-1));

      if (g == 0.0) {
        return cum_pop_[k-1] + N_left * delta;
      } else {
        // (exp(gamma_{k-1})/g) * expm1(g*delta)
        return cum_pop_[k-1] + (N_left / g) * std::expm1(g * delta);
      }
    }
  }
  CHECK(false) << "unrecognized type " << static_cast<int>(type_);
}

auto Skygrid_pop_model::inverse_cum_pop(double P) const -> double {
  auto m = M();

  // Region at or below x_0: P <= 0
  if (P <= 0.0) {
    return x(0) + P / std::exp(gamma(0));
  }

  // Region at or above x_M: P >= cum_pop_[M]
  if (P >= cum_pop_[m]) {
    return x(m) + (P - cum_pop_[m]) / std::exp(gamma(m));
  }

  // Interior: 0 < P < cum_pop_[M]
  // Binary search to find smallest k such that cum_pop_[k] > P
  auto it = std::upper_bound(cum_pop_.begin(), cum_pop_.end(), P);
  CHECK(it != cum_pop_.begin());
  CHECK(it != cum_pop_.end());
  auto k = static_cast<int>(std::distance(cum_pop_.begin(), it));
  CHECK_GE(k, 1);
  CHECK_LE(k, m);

  auto remaining = P - cum_pop_[k-1];
  CHECK_GE(remaining, 0.0);

  switch (type_) {
    case Type::k_staircase: {
      auto result = x(k-1) + remaining / std::exp(gamma(k));
      CHECK_GE(result, x(k-1));
      CHECK_LE(result, x(k));
      return result;
    }

    case Type::k_log_linear: {
      auto w = x(k) - x(k-1);
      auto g = (gamma(k) - gamma(k-1)) / w;
      auto N_left = std::exp(gamma(k-1));

      auto delta = 0.0;
      if (g == 0.0) {
        delta = remaining / N_left;
      } else {
        // delta = log1p(remaining * g / exp(gamma_{k-1})) / g
        delta = std::log1p(remaining * g / N_left) / g;
      }
      CHECK_GE(delta, -1e-10);
      CHECK_LE(delta, w + 1e-10);
      return x(k-1) + delta;
    }
  }
  CHECK(false) << "unrecognized type " << static_cast<int>(type_);
}

auto Skygrid_pop_model::print_to(std::ostream& os) const -> void {
  os << "Skygrid_pop_model{";
  os << absl::StreamFormat("type=%s, ",
                           (type_ == Type::k_staircase ? "k_staircase" :
                            type_ == Type::k_log_linear ? "k_log_linear" :
                            "???"));
  os << "ln N(t)=[";
  for (auto k = 0; k <= M(); ++k) {
    if (k != 0) { os << "; "; }
    os << absl::StreamFormat("at t=%g, %g", x(k), gamma(k));
  }
  os << "]}";
}

}  // namespace sapling
