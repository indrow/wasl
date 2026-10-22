// Modal DG residual for linear scalar advection on a 1D line mesh.
//
// The semi-discrete weak DG equation for element K and mode i is
//   J dc_i/dt = volume_i + surface_i
// where J = dx/d xi is the constant element Jacobian, the basis is orthonormal
// (mass matrix = J I, so the mass inverse is the scalar 1/J), and
//   volume_i  = sum_q w_q F(U_h(xi_q)) phi_i'(xi_q)
//   surface_i = -phi_i(+1) F*_R + phi_i(-1) F*_L
// The Jacobian cancels inside the volume integral (dx = J d xi against d/dx =
// (1/J) d/d xi), so it does not appear in volume_i. F*_R and F*_L are the +x
// upwind numerical fluxes at the element's right and left faces.
//
// This function returns dc_i/dt directly, i.e. it applies the 1/J mass inverse,
// so a forward Euler update is c_i += dt * rhs_i. This is the quantity SSP-RK3
// composes.
//
// State layout matches NUMERICAL_CONVENTIONS.md for one variable, indexed
//   U[element][mode]  =>  index = element * mode_count + mode
// A scalar equation has variable_count = 1, so the variable axis is dropped.
//
// Face accumulation is shared-endpoint: each interior face subtracts phi(+1) F*
// from its left owner and adds phi(-1) F* to its right owner, so a periodic
// reassembly conserves total mass to round-off. The left boundary face carries
// outward normal -1 and is an inflow face at xi = -1; the right boundary and all
// interior faces carry outward normal +1. Boundary faces have invalid_element on
// the side they do not own; the caller supplies their exterior trace through the
// same face loop, so periodic, inflow, and outflow conditions are uniform.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line_geometry.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/quadrature/quadrature.hpp"

namespace wasl {

// Volume-rule selection. The volume integrand is F(U_h) * phi_i', where for
// linear advection F(U_h) = a * U_h is degree p and phi_i' is degree p - 1, so
// the product is degree 2p - 1. The smallest n-point Gauss rule with
// 2n - 1 >= 2p - 1 has n = max(1, p) points.
inline int scalar_advection_volume_quadrature_points(int polynomial_order) noexcept {
  return polynomial_order > 0 ? polynomial_order : 1;
}

// Evaluate the state at a reference point by projecting modal coefficients.
// modes points to the mode_count coefficients of one element.
inline double evaluate_state(const double *modes, int polynomial_order, double xi) {
  const int mode_count = static_cast<int>(line_mode_count(polynomial_order));
  std::vector<double> basis(mode_count);
  legendre_basis(xi, polynomial_order, basis.data());
  double value = 0.0;
  for (int i = 0; i < mode_count; ++i) {
    value += modes[i] * basis[i];
  }
  return value;
}

// Compute the semi-discrete right-hand side dc_i/dt for every element of a line
// mesh, given the interior modal coefficients and the exterior trace state at
// each face. coefficients and rhs index U[element][mode].
//
// exterior_state[f] is the state seen from outside face f in its canonical
// direction (from the face's left element toward its right element). For an
// interior face it is the right element's interior trace at that endpoint; the
// upwind flux decides which side is used. For a boundary face the caller sets it
// to the physical boundary value: periodic neighbor trace, inflow data, or an
// outflow copy (the outflow value is unused when wind blows outward).
//
// On return, rhs holds dc_i/dt for each (element, mode). The residual preserves a
// constant state exactly: the volume integral of a constant flux against a
// gradient is exactly cancelled by the matching constant surface fluxes.
inline void scalar_advection_residual(const LineMesh &mesh, const ScalarAdvection &eq,
                                      int polynomial_order, const double *coefficients,
                                      const double *exterior_state, double *rhs) {
  const std::size_t element_count = mesh.num_elements();
  const int mode_count = static_cast<int>(line_mode_count(polynomial_order));

  const int quad_points = scalar_advection_volume_quadrature_points(polynomial_order);
  const QuadratureRule volume_rule = gauss_legendre(quad_points);

  for (std::size_t i = 0; i < element_count * mode_count; ++i) {
    rhs[i] = 0.0;
  }

  // Endpoint basis values are identical for every element (reference basis at
  // xi = -1, +1), so evaluate them once.
  std::vector<double> basis_at_left(mode_count);  // phi_i(-1)
  std::vector<double> basis_at_right(mode_count); // phi_i(+1)
  legendre_basis(-1.0, polynomial_order, basis_at_left.data());
  legendre_basis(1.0, polynomial_order, basis_at_right.data());

  for (std::size_t e = 0; e < element_count; ++e) {
    const double *modes = coefficients + e * mode_count;
    double *element_rhs = rhs + e * mode_count;

    // Volume term: sum_q w_q F(U_h(xi_q)) phi_i'(xi_q). The Jacobian from dx and
    // the 1/J from d/dx cancel for the affine line element, so the physical
    // Jacobian does not appear here.
    std::vector<double> grad(mode_count);
    std::vector<double> basis(mode_count);
    for (int q = 0; q < volume_rule.size(); ++q) {
      const double xi = volume_rule.points[q];
      const double w = volume_rule.weights[q];
      legendre_basis(xi, polynomial_order, basis.data());
      const double state = evaluate_state(modes, polynomial_order, xi);
      const double flux = scalar_advection_flux(eq, state);
      legendre_basis_gradient(xi, polynomial_order, grad.data());
      for (int i = 0; i < mode_count; ++i) {
        element_rhs[i] += w * flux * grad[i];
      }
    }
  }

  // Surface term. Each face sits between left_element and (when valid)
  // right_element. The stored outward_normal distinguishes the left boundary
  // face (normal -1, owned at its left end xi = -1, an inflow face) from every
  // other face (normal +1, owned at its right end xi = +1).
  for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
    const LineFace &face = mesh.faces[f];
    const std::int32_t left = face.left_element;
    const std::int32_t right = face.right_element;
    const double *left_modes = coefficients + static_cast<std::size_t>(left) * mode_count;
    double *left_rhs = rhs + static_cast<std::size_t>(left) * mode_count;

    if (face.outward_normal > 0.0) {
      // Owner sees this face at its right end (xi = +1, outflow). The +x upwind
      // flux takes the interior trace as the left state and the exterior as the
      // right state.
      double state_interior = 0.0;
      for (int i = 0; i < mode_count; ++i) {
        state_interior += left_modes[i] * basis_at_right[i];
      }
      const double flux = upwind_flux(eq, state_interior, exterior_state[f]);
      for (int i = 0; i < mode_count; ++i) {
        left_rhs[i] -= basis_at_right[i] * flux;
      }
      if (right != invalid_element) {
        // Right owner sees this face at its left end (xi = -1, inflow). The same
        // +x flux feeds it; the sign flip makes the shared-face contribution
        // cancel across owners.
        double *right_rhs = rhs + static_cast<std::size_t>(right) * mode_count;
        for (int i = 0; i < mode_count; ++i) {
          right_rhs[i] += basis_at_left[i] * flux;
        }
      }
    } else {
      // Left boundary face (outward normal -1): owner sees it at its left end
      // (xi = -1, inflow). The exterior is the left state and the interior trace
      // is the right state of the +x upwind flux.
      double state_interior = 0.0;
      for (int i = 0; i < mode_count; ++i) {
        state_interior += left_modes[i] * basis_at_left[i];
      }
      const double flux = upwind_flux(eq, exterior_state[f], state_interior);
      for (int i = 0; i < mode_count; ++i) {
        left_rhs[i] += basis_at_left[i] * flux;
      }
    }
  }

  // Apply the orthonormal mass inverse (1/J) so rhs is dc_i/dt. J is constant
  // per element.
  for (std::size_t e = 0; e < element_count; ++e) {
    const double inv_jac = 1.0 / element_map(mesh, e).jacobian();
    double *element_rhs = rhs + e * mode_count;
    for (int i = 0; i < mode_count; ++i) {
      element_rhs[i] *= inv_jac;
    }
  }
}

} // namespace wasl
