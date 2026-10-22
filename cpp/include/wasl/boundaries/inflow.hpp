// Inflow boundary policy for 1D scalar advection.
//
// Produces an ExteriorFunction for ssp_rk3_step that prescribes a fixed state
// at both boundary faces. The interior faces use the right element's
// left-endpoint trace, as the residual operator expects. For a one-sided inflow
// the caller pairs this with an outflow policy on the other boundary by
// constructing a combined callback; the common two-sided inflow case is covered
// directly here.
//
// At an outflow boundary (wind blowing outward) the prescribed value is unused
// because the upwind flux selects the interior trace, so prescribing inflow on
// both ends is harmless when only one end is physically inflow.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

inline EulerConserved euler_inflow_state(const EulerConserved &prescribed) noexcept {
  return prescribed;
}

// Build an inflow exterior-trace callback. inflow_value is the fixed state
// imposed at both boundary faces.
inline ExteriorFunction make_inflow_exterior(const LineMesh &mesh, int polynomial_order,
                                             double inflow_value) {
  return [&mesh, polynomial_order,
          inflow_value](const std::vector<double> &coefficients) -> std::vector<double> {
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
        // Boundary face: prescribe the fixed inflow state.
        exterior[f] = inflow_value;
      }
    }
    return exterior;
  };
}

} // namespace wasl
