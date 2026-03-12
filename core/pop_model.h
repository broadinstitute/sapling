#ifndef SAPLING_POP_MODEL_H_
#define SAPLING_POP_MODEL_H_

#include <ostream>
#include <vector>

#include "absl/strings/str_format.h"

namespace sapling {

class Pop_model {
 public:
  virtual ~Pop_model() = default;

  // N(t) * rho.  Time is measured in calendar units and moves forwards in time from an arbitrary epoch
  // The constant rho is the generation time.
  // Should be positive, but can be zero in special cases (e.g., N(t) = {1 - t for t <= 1, 0 otherwise}).
  // Should handle any value of t, not just nonnegative ones
  virtual auto pop_at_time(double t) const -> double = 0;

  // \int_c^t dt' N(t') * rho  [lower bound unspecified, but always the same]
  virtual auto cum_pop_at_time(double t) const -> double = 0;

  // I(t) := \int_c^t dt' 1/N(t') [lower bound unspecified, but always the same]
  // Can be +/-INFINITY.
  virtual auto intensity_at_time(double t) const -> double = 0;

  // I^{-1}(I) := t such that I(t) = I
  // Can be +/-INFINITY.
  virtual auto inverse_intensity(double I) const -> double = 0;

  // Inverse of cum_pop_at_time: find t such that cum_pop_at_time(t) = P
  virtual auto inverse_cum_pop(double P) const -> double = 0;

  // Debug printing
  friend auto operator<<(std::ostream& os, const Pop_model& pop_model) -> std::ostream& {
    pop_model.print_to(os);
    return os;
  }

 private:
  // Debug printing details
  virtual auto print_to(std::ostream& os) const -> void = 0;
};

class Const_pop_model : public Pop_model {
 public:
  explicit Const_pop_model(double pop);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto pop() const -> double { return pop_; }

  auto pop_at_time(double t) const -> double override;
  auto cum_pop_at_time(double t) const -> double override;
  auto intensity_at_time(double t) const -> double override;
  auto inverse_intensity(double I) const -> double override;
  auto inverse_cum_pop(double P) const -> double override;

 private:
  double pop_;

  auto print_to(std::ostream& os) const -> void override { os << absl::StreamFormat("Const_pop_model{pop=%g}", pop()); }
};

class Exp_pop_model : public Pop_model {
 public:
  explicit Exp_pop_model(double t0, double pop_at_t0, double growth_rate);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto t0() const -> double { return t0_; }
  auto pop_at_t0() const -> double { return pop_at_t0_; }
  auto growth_rate() const -> double { return growth_rate_; }

  auto pop_at_time(double t) const -> double override;
  auto cum_pop_at_time(double t) const -> double override;
  auto intensity_at_time(double t) const -> double override;
  auto inverse_intensity(double I) const -> double override;
  auto inverse_cum_pop(double P) const -> double override;

 private:
  double t0_;
  double pop_at_t0_;
  double growth_rate_;

  auto print_to(std::ostream& os) const -> void override {
    os << absl::StreamFormat("Exp_pop_model{t0=%g, n0=%g, g=%g}",
                             t0(), pop_at_t0(), growth_rate());
  }
};

// A Skygrid model has a population defined piecewise over a finite
// number of predefined intervals:
//
//
//                      |         |       |           |       |
//                      |         |       |           |       |
//         Int.0        |  Int.1  | Int.2 |    ...    | Int.M |   Int.{M+1}
//                      |         |       |           |       |
//                      |         |       |           |       |
//    -INFINITY ...-----+---------+-------+----....---+-------+---... +INFINITY  --> t
//                     x_0       x_1     x_2       x_{M-1}   x_M
//
// The time range (-INFINITY, +INFINITY) is partitioned into M+2 intervals,
// the first and last of which are open-ended.  The boundaries between consecutive intervals
// are called "knots", and are specified at x_0 < x_1 < ... < x_M.  The interval
// for a given time t is defined as:
//
//               /   0,                     t <= x_0;
// Interval(t) = |   k,           x_{k-1} < t <= x_k;    (1 <= k <= M)
//               \   M+1,             x_M < t.
//
// The log-population-size log(N(t)) is specified at each of the knots:
//
//     log N(t_k) =: gamma_k     0 <= k <= M.
//
// We support two variants for interpolating the population size at all other times.
//
// In the traditional Skygrid model (Gill et al 2012, BEAST's "gmrfSkyGridLikelihood",
// k_staircase below), the population is constant over the length of a single interval,
// but can vary from interval to interval:
//
//  Staircase
//  ---------
//                 / gamma_0,            t <= x_0;       // Interval 0
//      log N(t) = | gamma_k,  x_{k-1} < t <= x_k;       // Interval 1 <= k <= M
//                 \ gamma_M,      x_M < t.              // Interval M+1
//
// A somewhat more natural choice (k_log_linear below) is to have the population
// grow exponentially within a single interval, be continuous across intervals,
// but allowing different growth rates (possibly negative) in different intervals.
// We parametrize this curve as follows:
//
//  Log-linear
//  ----------
//                 / gamma_0,                                  t <= x_0;  // Interval 0
//      log N(t) = | (1-c) gamma_{k-1} + c gamma_k,  x_{k-1} < t <= x_k;  // Interval 1 <= k <= M
//                 \ gamma_M,                            x_M < t.         // Interval M+1
//
//         [N.B.: for intervals 1 <= k <= M, we have
//                  t = (1-c) x_{k-1} + c x_k
//               => c = (t - x_{k-1}) / (x_k - x_{k-1}) ]
//
class Skygrid_pop_model : public Pop_model {
 public:
  enum class Type {

    // Staircase = the traditional Skygrid model from Gill et al 2012 and BEAST.
    // We allow arbitrary knot times satisfying x_0 < x_1 < ... < x_M.
    // The main notational discrepancy is that, unlike in Gill et al 2012 or in BEAST,
    // Sapling's internal time axis increases towards the future and has a fixed epoch:
    //
    //  Here                 Gill et al 2012
    //  ----                 ---------------
    //  t                    T-t
    //  M                    M
    //  T                    0
    //  K                    K
    //  x_k                  x_{M-k}
    //  gamma_k              gamma_{M+1-k}
    //  exp(gamma_k)         theta_{M+1-k}

    k_staircase = 1,

    // Log-linear = A continuous population curve s.t. log(N(t)) is continuous
    // and linearly interpolates between the values at the knots.

    k_log_linear = 2
  };

  Skygrid_pop_model(
      std::vector<double> x,      // x[k] = time of knot k
      std::vector<double> gamma,  // gamma[k] = log(N(x[k]))
      Type type);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto x() const -> const std::vector<double>& { return x_; }
  auto gamma() const -> const std::vector<double>& { return gamma_; }
  auto type() const -> Type { return type_; }

  // Derived quantities in convenient notation (see notation above)
  auto M() const -> int { return std::ssize(x_) - 1; }
  auto x(int k) const -> double { return x_.at(k); }
  auto gamma(int k) const -> double { return gamma_.at(k); }

  // Pop_model interface
  auto pop_at_time(double t) const -> double override;
  auto cum_pop_at_time(double t) const -> double override;
  auto intensity_at_time(double t) const -> double override;
  auto inverse_intensity(double I) const -> double override;
  auto inverse_cum_pop(double P) const -> double override;

  // Interval(t) in the top-level comment for Skygrid_pop_model
  auto interval_containing_t(double t) const -> int;

 private:
  std::vector<double> x_;
  std::vector<double> gamma_;
  Type type_;

  // Precomputed at knot boundaries (size = M+1, indices 0..M)
  // cum_intensity_[0] = 0, cum_intensity_[k] = integral of 1/N(t) from x_0 to x_k
  std::vector<double> cum_intensity_;
  // cum_pop_[0] = 0, cum_pop_[k] = integral of N(t) from x_0 to x_k
  std::vector<double> cum_pop_;

  auto print_to(std::ostream& os) const -> void override;
};

}  // namespace sapling

#endif // SAPLING_POP_MODEL_H_
