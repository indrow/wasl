// Point evaluation of a modal expansion u(x) = sum_i c_i phi_i(x).

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"

namespace wasl {

// Sum c_i b_i over mode_count modes, accumulated in mode order.
inline double modal_sum(const double *coefficients, const double *basis,
                        std::size_t mode_count) noexcept {
  double value = 0.0;
  for (std::size_t mode = 0; mode < mode_count; ++mode)
    value += coefficients[mode] * basis[mode];
  return value;
}

// Value of one line element's Legendre expansion at reference point xi.
inline double line_expansion_value(const double *modes, int polynomial_order, double xi) {
  const std::size_t mode_count = line_mode_count(polynomial_order);
  std::vector<double> basis(mode_count);
  legendre_basis(xi, polynomial_order, basis.data());
  return modal_sum(modes, basis.data(), mode_count);
}

} // namespace wasl
