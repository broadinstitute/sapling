#ifndef SAPLING_POP_MODEL_H_
#define SAPLING_POP_MODEL_H_

#include <ostream>

#include "absl/random/bit_gen_ref.h"
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

  // Sample the time of a random individual between min_t and max_t
  virtual auto sample(double min_t, double max_t, absl::BitGenRef prng) const -> double = 0;

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
  auto sample(double min_t, double max_t, absl::BitGenRef prng) const -> double override;

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
  auto sample(double min_t, double max_t, absl::BitGenRef prng) const -> double override;

 private:
  double t0_;
  double pop_at_t0_;
  double growth_rate_;
  
  auto print_to(std::ostream& os) const -> void override {
    os << absl::StreamFormat("Exp_pop_model{t0=%g, n0=%g, g=%g}",
                             t0(), pop_at_t0(), growth_rate());
  }
};

}  // namespace sapling

#endif // SAPLING_POP_MODEL_H_
