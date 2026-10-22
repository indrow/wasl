// Number of modal coefficients of a complete polynomial space of order p.

#pragma once

#include <cstddef>

namespace wasl {

constexpr std::size_t line_mode_count(std::size_t polynomial_order) { return polynomial_order + 1; }

constexpr std::size_t triangle_mode_count(std::size_t polynomial_order) {
  return (polynomial_order + 1) * (polynomial_order + 2) / 2;
}

} // namespace wasl
