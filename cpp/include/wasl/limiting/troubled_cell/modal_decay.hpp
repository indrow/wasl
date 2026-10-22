// Persson-style modal-decay troubled-cell sensor.
//
// The sensor is the energy in the highest-degree modes divided by the total
// modal energy. P0 has no resolvable decay and is never marked.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/basis/mode_count.hpp"

namespace wasl {

struct ModalEnergy {
  double total = 0.0;
  double highest_degree = 0.0;
};

// Energies of mode_count coefficients spaced stride apart; modes from
// highest_degree_begin onward have the highest polynomial degree.
inline ModalEnergy modal_energy(const double *coefficients, std::size_t stride,
                                std::size_t mode_count, std::size_t highest_degree_begin) noexcept {
  ModalEnergy energy;
  for (std::size_t mode = 0; mode < mode_count; ++mode) {
    const double value = coefficients[mode * stride];
    energy.total += value * value;
    if (mode >= highest_degree_begin)
      energy.highest_degree += value * value;
  }
  return energy;
}

// Sensor ratio of one contiguous triangle polynomial.
inline double modal_decay_ratio(const double *coefficients, int polynomial_order) noexcept {
  if (polynomial_order < 1)
    return 0.0;
  const auto energy = modal_energy(coefficients, 1, triangle_mode_count(polynomial_order),
                                   triangle_mode_count(polynomial_order - 1));
  return energy.total > 0.0 ? energy.highest_degree / energy.total : 0.0;
}

inline bool modal_decay_troubled(const double *coefficients, int polynomial_order,
                                 double threshold) noexcept {
  return modal_decay_ratio(coefficients, polynomial_order) > threshold;
}

// Indices of flagged elements in ascending order.
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
