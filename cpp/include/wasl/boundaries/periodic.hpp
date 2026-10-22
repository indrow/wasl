// Periodic boundary policy for 1D scalar advection.
//
// Produces an ExteriorFunction for ssp_rk3_step that treats the domain as a
// ring: the left boundary face wraps to the last element's right-endpoint
// trace, and the right boundary face wraps to the first element's left-endpoint
// trace. Interior faces use the right element's left-endpoint trace, as the
// residual operator expects. With periodic boundaries the total mass is
// conserved and a smooth pulse returns to its initial state after one full
// revolution.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

// Build a periodic exterior-trace callback for the given mesh and polynomial
// order. The returned function recomputes traces from whatever coefficients the
// stepper hands it, so it is valid at every RK stage.
inline ExteriorFunction make_periodic_exterior(const LineMesh &mesh, int polynomial_order) {
  return [&mesh, polynomial_order](const std::vector<double> &coefficients) -> std::vector<double> {
    const int mode_count = static_cast<int>(line_mode_count(polynomial_order));
    const std::size_t n = mesh.num_elements();
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
      } else if (mesh.faces[f].outward_normal < 0.0) {
        // Left boundary wraps to the last element's right trace (xi = +1).
        legendre_basis(1.0, polynomial_order, basis.data());
        const double *modes = coefficients.data() + (n - 1) * mode_count;
        double trace = 0.0;
        for (int j = 0; j < mode_count; ++j) {
          trace += modes[j] * basis[j];
        }
        exterior[f] = trace;
      } else {
        // Right boundary wraps to the first element's left trace (xi = -1).
        legendre_basis(-1.0, polynomial_order, basis.data());
        const double *modes = coefficients.data();
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
