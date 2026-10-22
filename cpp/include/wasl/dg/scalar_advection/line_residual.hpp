// Modal DG residual for linear scalar advection on a 1D line mesh.
//
// The semi-discrete weak DG equation for element K and mode i is
//   J dc_i/dt = volume_i + surface_i
// where J = dx/d xi is the constant element Jacobian, the basis is orthonormal
// (mass matrix = J I, so the mass inverse is the scalar 1/J), and
//   volume_i  = sum_q w_q F(U_h(xi_q)) phi_i'(xi_q)
//   surface_i = -phi_i(+1) F*_R + phi_i(-1) F*_L.
// The Jacobian cancels inside the volume integral (dx = J d xi against
// d/dx = (1/J) d/d xi). F*_R and F*_L are the +x upwind numerical fluxes at the
// element's right and left faces. The residual applies 1/J, so it returns
// dc_i/dt and a forward-Euler update is c_i += dt * rhs_i.
//
// Coefficients and rhs use the scalar solver layout U[element][mode].
//
// Each interior face subtracts phi(+1) F* from its left owner and adds
// phi(-1) F* to its right owner, so periodic reassembly conserves mass to
// round-off. The left boundary face carries outward normal -1 and is owned at
// xi = -1; every other face carries +1 and is owned at xi = +1. The caller
// supplies boundary exterior traces through the same face loop, so periodic,
// inflow, and outflow conditions are uniform.

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/modal_expansion.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/quadrature/gauss_legendre.hpp"

namespace wasl {

// The volume integrand a U_h phi_i' has degree 2p - 1, which the smallest
// Gauss rule with 2n - 1 >= 2p - 1 integrates exactly: n = max(1, p).
inline int scalar_advection_volume_quadrature_points(int polynomial_order) noexcept {
  return polynomial_order > 0 ? polynomial_order : 1;
}

// exterior_state[f] is the state outside face f in its canonical direction:
// the right element's trace for an interior face, or the boundary value
// (periodic neighbor trace, inflow data, or outflow copy) for a boundary face.
// The residual preserves a constant state exactly.
inline void scalar_advection_residual(const LineMesh &mesh, const ScalarAdvection &eq,
                                      int polynomial_order, const double *coefficients,
                                      const double *exterior_state, double *rhs) {
  const std::size_t element_count = mesh.num_elements();
  const std::size_t mode_count = line_mode_count(polynomial_order);
  const QuadratureRule volume_rule =
      gauss_legendre(scalar_advection_volume_quadrature_points(polynomial_order));
  std::fill(rhs, rhs + element_count * mode_count, 0.0);

  std::vector<double> basis(mode_count), gradient(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const double *modes = coefficients + element * mode_count;
    double *element_rhs = rhs + element * mode_count;
    for (int q = 0; q < volume_rule.size(); ++q) {
      const double xi = volume_rule.points[q];
      legendre_basis(xi, polynomial_order, basis.data());
      legendre_basis_gradient(xi, polynomial_order, gradient.data());
      const double flux = scalar_advection_flux(eq, modal_sum(modes, basis.data(), mode_count));
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        element_rhs[mode] += volume_rule.weights[q] * flux * gradient[mode];
    }
  }

  // Endpoint basis values are shared by every element.
  std::vector<double> basis_at_left(mode_count), basis_at_right(mode_count);
  legendre_basis(-1.0, polynomial_order, basis_at_left.data());
  legendre_basis(1.0, polynomial_order, basis_at_right.data());

  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const LineFace &face = mesh.faces[face_index];
    const double *left_modes =
        coefficients + static_cast<std::size_t>(face.left_element) * mode_count;
    double *left_rhs = rhs + static_cast<std::size_t>(face.left_element) * mode_count;

    if (face.outward_normal > 0.0) {
      // The owner sees this face at xi = +1: its trace is the upwind left state.
      const double interior = modal_sum(left_modes, basis_at_right.data(), mode_count);
      const double flux = upwind_flux(eq, interior, exterior_state[face_index]);
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        left_rhs[mode] -= basis_at_right[mode] * flux;
      if (!face.is_boundary()) {
        // The right owner sees the same flux at xi = -1 with opposite sign, so
        // the shared-face contributions cancel.
        double *right_rhs = rhs + static_cast<std::size_t>(face.right_element) * mode_count;
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          right_rhs[mode] += basis_at_left[mode] * flux;
      }
    } else {
      // Left boundary: the owner sees it at xi = -1 and the exterior is the
      // upwind left state.
      const double interior = modal_sum(left_modes, basis_at_left.data(), mode_count);
      const double flux = upwind_flux(eq, exterior_state[face_index], interior);
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        left_rhs[mode] += basis_at_left[mode] * flux;
    }
  }

  for (std::size_t element = 0; element < element_count; ++element) {
    const double inverse_jacobian = 1.0 / element_map(mesh, element).jacobian();
    double *element_rhs = rhs + element * mode_count;
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      element_rhs[mode] *= inverse_jacobian;
  }
}

} // namespace wasl
