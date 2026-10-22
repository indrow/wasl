// A one-dimensional quadrature rule on the reference line xi in [-1, 1].

#pragma once

#include <vector>

namespace wasl {

struct QuadratureRule {
  std::vector<double> points;
  std::vector<double> weights;

  int size() const noexcept { return static_cast<int>(points.size()); }
};

} // namespace wasl
