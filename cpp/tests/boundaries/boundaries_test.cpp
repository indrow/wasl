// Tests for the 1D scalar-advection boundary policies.
//
// Each policy produces an ExteriorFunction whose per-face exterior trace must
// match the documented contract: periodic wraps boundary faces to the opposite
// end of the domain, inflow prescribes a fixed value, and outflow copies the
// owner's own interior trace. Interior faces are identical across all three
// policies (only boundary faces differ).

#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/inflow.hpp"
#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/boundaries/outflow.hpp"
#include "wasl/boundaries/periodic.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line/line_geometry.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/time/cfl_timestep.hpp"
#include "wasl/time/line_advection_step.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Build a uniform 4-element mesh on [0, 1]. Faces: 0=left boundary, 1-3=interior,
// 4=right boundary.
wasl::LineMesh make_uniform_mesh(double length, int n) {
  wasl::LineMesh mesh;
  const double h = length / n;
  mesh.elements.reserve(n);
  for (int e = 0; e < n; ++e) {
    mesh.elements.push_back(wasl::LineElement{e, e + 1, e * h, (e + 1) * h, 1});
  }
  mesh.faces.reserve(n + 1);
  mesh.faces.push_back(wasl::LineFace{0, wasl::invalid_element, 0.0, -1.0, wasl::invalid_physical});
  for (int e = 0; e + 1 < n; ++e) {
    mesh.faces.push_back(wasl::LineFace{e, e + 1, (e + 1) * h, 1.0, wasl::invalid_physical});
  }
  mesh.faces.push_back(
      wasl::LineFace{n - 1, wasl::invalid_element, length, 1.0, wasl::invalid_physical});
  return mesh;
}

// Interior trace of an element at a reference endpoint.
double element_trace(const std::vector<double> &coefficients, int polynomial_order, int element,
                     double xi) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  std::vector<double> basis(mode_count);
  wasl::legendre_basis(xi, polynomial_order, basis.data());
  const double *modes = coefficients.data() + static_cast<std::size_t>(element) * mode_count;
  double trace = 0.0;
  for (int j = 0; j < mode_count; ++j) {
    trace += modes[j] * basis[j];
  }
  return trace;
}

// Periodic: left boundary wraps to the last element's right trace; right
// boundary wraps to the first element's left trace.
bool check_periodic_wrapping() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  constexpr int p = 1;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  // Distinct coefficients so traces are distinguishable.
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = static_cast<double>(e) + 1.0;
    coefficients[e * mode_count + 1] = 0.3 * static_cast<double>(e);
  }

  const auto exterior_fn = wasl::make_periodic_exterior(mesh, p);
  const auto exterior = exterior_fn(coefficients);

  // Left boundary (face 0) must equal the last element's right trace.
  const double last_right = element_trace(coefficients, p, 3, 1.0);
  if (!approx(exterior[0], last_right))
    return false;
  // Right boundary (face 4) must equal the first element's left trace.
  const double first_left = element_trace(coefficients, p, 0, -1.0);
  if (!approx(exterior[4], first_left))
    return false;
  // Interior face 1 must equal element 1's left trace.
  const double e1_left = element_trace(coefficients, p, 1, -1.0);
  if (!approx(exterior[1], e1_left))
    return false;
  return true;
}

// Inflow: boundary faces carry the prescribed value regardless of interior state.
bool check_inflow_prescription() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  constexpr int p = 2;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));
  const double inflow = 5.0;

  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = static_cast<double>(e) * std::sqrt(2.0);
    coefficients[e * mode_count + 2] = 0.1;
  }

  const auto exterior_fn = wasl::make_inflow_exterior(mesh, p, inflow);
  const auto exterior = exterior_fn(coefficients);

  // Both boundary faces must equal the prescribed inflow.
  if (!approx(exterior[0], inflow))
    return false;
  if (!approx(exterior[4], inflow))
    return false;
  // Interior face 2 must equal element 2's left trace (not the inflow).
  const double e2_left = element_trace(coefficients, p, 2, -1.0);
  if (!approx(exterior[2], e2_left))
    return false;
  return true;
}

// Outflow: boundary faces copy the owner's own interior trace (zero-gradient).
bool check_outflow_zero_gradient() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  constexpr int p = 1;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = static_cast<double>(e) + 1.0;
    coefficients[e * mode_count + 1] = 0.25;
  }

  const auto exterior_fn = wasl::make_outflow_exterior(mesh, p);
  const auto exterior = exterior_fn(coefficients);

  // Left boundary (face 0, owner element 0) copies element 0's left trace.
  const double e0_left = element_trace(coefficients, p, 0, -1.0);
  if (!approx(exterior[0], e0_left))
    return false;
  // Right boundary (face 4, owner element 3) copies element 3's right trace.
  const double e3_right = element_trace(coefficients, p, 3, 1.0);
  if (!approx(exterior[4], e3_right))
    return false;
  // Interior face 3 must equal element 3's left trace.
  const double e3_left = element_trace(coefficients, p, 3, -1.0);
  if (!approx(exterior[3], e3_left))
    return false;
  return true;
}

// All three policies agree on interior faces; only boundary faces differ.
bool check_interior_faces_match() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  constexpr int p = 2;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = std::cos(static_cast<double>(e));
    coefficients[e * mode_count + 1] = 0.4;
    coefficients[e * mode_count + 2] = 0.1 * static_cast<double>(e);
  }

  const auto periodic = wasl::make_periodic_exterior(mesh, p)(coefficients);
  const auto inflow = wasl::make_inflow_exterior(mesh, p, 9.0)(coefficients);
  const auto outflow = wasl::make_outflow_exterior(mesh, p)(coefficients);

  // Faces 1, 2, 3 are interior and must match across all policies.
  for (std::size_t f = 1; f <= 3; ++f) {
    if (!approx(periodic[f], inflow[f]))
      return false;
    if (!approx(periodic[f], outflow[f]))
      return false;
  }
  return true;
}

// Periodic conservation: after one SSP-RK3 step with periodic boundaries, the
// total mode-zero mass (volume-weighted) is unchanged.
bool check_periodic_conservation_one_step() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  const wasl::ScalarAdvection eq{1.0};
  constexpr int p = 1;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));
  const double dt = wasl::advection_cfl_timestep(mesh, eq, p, 0.2);

  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = std::cos(static_cast<double>(e)) * std::sqrt(2.0);
    coefficients[e * mode_count + 1] = 0.3;
  }

  auto mass = [&]() {
    double total = 0.0;
    for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
      const wasl::LineMap map = wasl::element_map(mesh, e);
      total += map.jacobian() * coefficients[e * mode_count + 0];
    }
    return total;
  };
  const double mass_before = mass();

  const auto exterior_fn = wasl::make_periodic_exterior(mesh, p);
  wasl::SspRk3Workspace ws;
  wasl::ssp_rk3_step(mesh, eq, p, dt, exterior_fn, coefficients, ws);

  return approx(mass(), mass_before, 1e-12);
}

// Mixed ends: the left face takes the inflow value, the right face its owner's
// own trace at xi = +1, and interior faces are unchanged.
bool check_mixed_inflow_outflow() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  constexpr int p = 2;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e)
    for (int j = 0; j < mode_count; ++j)
      coefficients[e * mode_count + j] = 0.5 + 0.1 * static_cast<double>(e) - 0.2 * j;

  const wasl::LineBoundaryConditions conditions{{wasl::BoundaryKind::inflow, 2.5},
                                                {wasl::BoundaryKind::outflow, 0.0}};
  const auto mixed = wasl::make_line_boundary_exterior(mesh, p, conditions)(coefficients);
  const auto periodic = wasl::make_periodic_exterior(mesh, p)(coefficients);
  const auto last = static_cast<int>(mesh.num_elements()) - 1;
  if (mixed.front() != 2.5)
    return false;
  if (!approx(mixed.back(), element_trace(coefficients, p, last, 1.0)))
    return false;
  for (std::size_t f = 1; f + 1 < mesh.faces.size(); ++f)
    if (mixed[f] != periodic[f])
      return false;
  return true;
}

} // namespace

int main() {
  if (!check_periodic_wrapping())
    return 1;
  if (!check_inflow_prescription())
    return 1;
  if (!check_outflow_zero_gradient())
    return 1;
  if (!check_interior_faces_match())
    return 1;
  if (!check_periodic_conservation_one_step())
    return 1;
  if (!check_mixed_inflow_outflow())
    return 1;
  return 0;
}
