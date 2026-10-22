// Outflow (zero-gradient) boundary policy for 1D scalar advection.
//
// Produces an ExteriorFunction for ssp_rk3_step that copies the owner element's
// own interior trace at each boundary face, so no gradient is imposed and the
// wave leaves the domain freely. For the right boundary (outward_normal +1) the
// exterior is the owner's right-endpoint trace; for the left boundary
// (outward_normal -1) it is the owner's left-endpoint trace. Interior faces use
// the right element's left-endpoint trace, as the residual operator expects.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

inline EulerConserved euler_outflow_state(const EulerConserved &interior) noexcept {
  return interior;
}

// Build an outflow (zero-gradient) exterior-trace callback. Both boundary faces
// copy their owner's own interior trace.
inline ExteriorFunction make_outflow_exterior(const LineMesh &mesh, int polynomial_order) {
  return [&mesh, polynomial_order](const std::vector<double> &coefficients) -> std::vector<double> {
    const int mode_count = static_cast<int>(line_mode_count(polynomial_order));
    std::vector<double> basis(mode_count);
    std::vector<double> exterior(mesh.faces.size(), 0.0);
    for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
      const auto right = mesh.faces[f].right_element;
      if (right != invalid_element) {
        // Interior face: right element's trace at xi = -1.
        legendre_basis(-1.0, polynomial_order, basis.data());
        const double *modes = coefficients.data() + static_cast<std::size_t>(right) * mode_count;
        double trace = 0.0;
        for (int j = 0; j < mode_count; ++j) {
          trace += modes[j] * basis[j];
        }
        exterior[f] = trace;
      } else {
        // Boundary face: copy the owner's own trace at the owned endpoint.
        const auto owner = mesh.faces[f].left_element;
        const double xi = (mesh.faces[f].outward_normal < 0.0) ? -1.0 : 1.0;
        legendre_basis(xi, polynomial_order, basis.data());
        const double *modes = coefficients.data() + static_cast<std::size_t>(owner) * mode_count;
        double trace = 0.0;
        for (int j = 0; j < mode_count; ++j) {
          trace += modes[j] * basis[j];
        }
        exterior[f] = trace;
      }
    }
    return exterior;
  };
}

} // namespace wasl
