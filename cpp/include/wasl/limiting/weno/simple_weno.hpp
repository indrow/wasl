#pragma once

#include <cstddef>
#include <vector>

namespace wasl {

inline double weno_smoothness(const double *coefficients, std::size_t mode_count) noexcept {
  double value = 0.0;
  for (std::size_t mode = 1; mode < mode_count; ++mode)
    value += coefficients[mode] * coefficients[mode];
  return value;
}

// Blend the target and transformed neighbor candidates. candidates is packed
// as [candidate][mode]. The target receives half the ideal weight and the
// immediate neighbors share the other half.
inline void simple_weno(const double *target, const double *candidates, std::size_t candidate_count,
                        std::size_t mode_count, double *limited, double epsilon = 1e-12) {
  if (candidate_count == 0) {
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      limited[mode] = target[mode];
    return;
  }

  std::vector<double> weights(candidate_count + 1);
  const double target_denominator = epsilon + weno_smoothness(target, mode_count);
  weights[0] = 0.5 / (target_denominator * target_denominator);
  double weight_sum = weights[0];
  const double neighbor_ideal_weight = 0.5 / static_cast<double>(candidate_count);
  for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
    const double *values = candidates + candidate * mode_count;
    const double denominator = epsilon + weno_smoothness(values, mode_count);
    weights[candidate + 1] = neighbor_ideal_weight / (denominator * denominator);
    weight_sum += weights[candidate + 1];
  }

  limited[0] = target[0];
  for (std::size_t mode = 1; mode < mode_count; ++mode) {
    double value = weights[0] * target[mode];
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate)
      value += weights[candidate + 1] * candidates[candidate * mode_count + mode];
    limited[mode] = value / weight_sum;
  }
}

} // namespace wasl
