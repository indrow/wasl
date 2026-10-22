#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/basis/reference_element.hpp"

namespace wasl {

// Persson-style modal decay sensor on one triangle polynomial. The returned
// ratio is the energy in modes of exactly degree p divided by total modal
// energy. P0 has no resolvable decay and is never marked.
inline double modal_decay_ratio(const double *coefficients, int polynomial_order) noexcept {
  if (polynomial_order < 1)
    return 0.0;

  const std::size_t mode_count = triangle_mode_count(polynomial_order);
  const std::size_t highest_degree_begin = triangle_mode_count(polynomial_order - 1);
  double total_energy = 0.0;
  double highest_degree_energy = 0.0;
  for (std::size_t mode = 0; mode < mode_count; ++mode) {
    const double energy = coefficients[mode] * coefficients[mode];
    total_energy += energy;
    if (mode >= highest_degree_begin)
      highest_degree_energy += energy;
  }
  return total_energy > 0.0 ? highest_degree_energy / total_energy : 0.0;
}

inline bool modal_decay_troubled(const double *coefficients, int polynomial_order,
                                 double threshold) noexcept {
  return modal_decay_ratio(coefficients, polynomial_order) > threshold;
}

inline std::vector<std::uint32_t> compact_troubled_elements(const std::uint8_t *flags,
                                                            std::size_t element_count) {
  std::vector<std::uint32_t> troubled;
  troubled.reserve(element_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    if (flags[element] != 0)
      troubled.push_back(static_cast<std::uint32_t>(element));
  }
  return troubled;
}

} // namespace wasl
