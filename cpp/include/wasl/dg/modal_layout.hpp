// Structure-of-arrays modal coefficient layout (NUMERICAL_CONVENTIONS.md):
//   U[variable][mode][element],
//   index = ((variable * mode_count) + mode) * element_count + element.
// Euler residuals, limiters, and checkpoints use this layout. The scalar
// solver path keeps U[element][mode] and converts at checkpoint boundaries.

#pragma once

#include <cmath>
#include <cstddef>

#include "wasl/equations/euler/euler.hpp"

namespace wasl {

struct ModalLayout {
  std::size_t mode_count;
  std::size_t element_count;

  std::size_t index(std::size_t variable, std::size_t mode, std::size_t element) const noexcept {
    return (variable * mode_count + mode) * element_count + element;
  }
};

// Euler state of one element at a point with the given basis values.
inline EulerConserved evaluate_euler_state(const ModalLayout &layout, const double *coefficients,
                                           std::size_t element, const double *basis) noexcept {
  EulerConserved state{0.0, 0.0, 0.0, 0.0};
  for (int variable = 0; variable < Euler::variable_count; ++variable) {
    double &component = euler_component(state, variable);
    for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
      component += coefficients[layout.index(static_cast<std::size_t>(variable), mode, element)] *
                   basis[mode];
  }
  return state;
}

// Cell average of one element: mode zero is phi_0 = 1/sqrt(2).
inline EulerConserved euler_cell_average(const ModalLayout &layout, const double *coefficients,
                                         std::size_t element) noexcept {
  const double scale = 1.0 / std::sqrt(2.0);
  return {scale * coefficients[layout.index(0, 0, element)],
          scale * coefficients[layout.index(1, 0, element)],
          scale * coefficients[layout.index(2, 0, element)],
          scale * coefficients[layout.index(3, 0, element)]};
}

// Multiply every coefficient of one element by factor, for example to apply the
// orthonormal mass inverse 1/det(J).
inline void scale_element(const ModalLayout &layout, std::size_t variable_count, double *values,
                          std::size_t element, double factor) noexcept {
  for (std::size_t variable = 0; variable < variable_count; ++variable)
    for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
      values[layout.index(variable, mode, element)] *= factor;
}

} // namespace wasl
