// Face quadrature of the reference line element.
//
// A 1D face is a single point, so its weight is one. The left face sits at
// xi = -1 (outward normal -1) and the right face at xi = +1 (outward normal +1).
// The rule shares one face-evaluation shape with the triangle edge rules.

#pragma once

#include <vector>

#include "wasl/quadrature/quadrature_rule.hpp"

namespace wasl {

inline QuadratureRule endpoint_faces() {
  return QuadratureRule{
      std::vector<double>{-1.0, 1.0},
      std::vector<double>{1.0, 1.0},
  };
}

} // namespace wasl
