// Tests for the modal DG scalar-advection residual on a 1D line mesh.
//
// These are the conservation, consistency, and constant-state checks required by
// a numerical change. Meshes are built directly so the tests do not depend on a
// Gmsh binary. Properties checked:
//   - Constant state gives zero residual (free-stream / constant preservation).
//   - P0 reduces exactly to first-order upwind finite volume.
//   - Shared interior endpoints show equal and opposite flux contributions.
//   - A periodic reassembly conserves total mass to round-off.

#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/dg/residual_operator.hpp"
#include "wasl/mesh/line_mesh.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Build a uniform line mesh on [0, L] with n elements. Faces are interior pairs
// plus two boundary faces at x = 0 and x = L. This mirrors the canonical mesh
// produced by build_line_mesh but avoids the MSH parsing dependency.
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

// Interior trace of a single element at a reference endpoint.
double element_trace(const std::vector<double> &coefficients, int polynomial_order, int element,
                     double xi) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  std::vector<double> basis(mode_count);
  wasl::legendre_basis(xi, polynomial_order, basis.data());
  const double *modes = coefficients.data() + static_cast<std::size_t>(element) * mode_count;
  double trace = 0.0;
  for (int i = 0; i < mode_count; ++i) {
    trace += modes[i] * basis[i];
  }
  return trace;
}

// Build the exterior trace vector for given boundary values. For interior faces
// the exterior is the right element's trace at its left endpoint (xi = -1). The
// left boundary face (normal -1) takes left_boundary; the right boundary face
// (normal +1) takes right_boundary. The exterior is the state seen from outside
// the face in its canonical left-to-right direction.
std::vector<double> build_exterior(const wasl::LineMesh &mesh, int polynomial_order,
                                   const std::vector<double> &coefficients, double left_boundary,
                                   double right_boundary) {
  std::vector<double> exterior(mesh.faces.size(), 0.0);
  for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
    const auto right = mesh.faces[f].right_element;
    if (right == wasl::invalid_element) {
      exterior[f] = (mesh.faces[f].outward_normal < 0.0) ? left_boundary : right_boundary;
    } else {
      exterior[f] = element_trace(coefficients, polynomial_order, right, -1.0);
    }
  }
  return exterior;
}

// Periodic exterior: the left boundary wraps to the rightmost element's right
// trace, and the right boundary wraps to the leftmost element's left trace, so
// the domain is a ring regardless of wavespeed direction.
std::vector<double> build_periodic_exterior(const wasl::LineMesh &mesh, int polynomial_order,
                                            const std::vector<double> &coefficients) {
  const std::size_t n = mesh.num_elements();
  const double wrap_left =
      element_trace(coefficients, polynomial_order, static_cast<int>(n - 1), 1.0);
  const double wrap_right = element_trace(coefficients, polynomial_order, 0, -1.0);
  return build_exterior(mesh, polynomial_order, coefficients, wrap_left, wrap_right);
}

// Constant-state preservation: a uniform field is a steady solution of advection
// with matching boundary data, so the residual must be zero at every mode and
// element. This is the free-stream preservation invariant.
bool check_constant_state_zero_residual() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  const wasl::ScalarAdvection eq{1.3};
  constexpr int p = 2;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  const double c = 0.7;
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    // A constant field U = c is the single mode-zero value c / phi_0, and phi_0
    // = 1/sqrt(2), so the modal coefficient is c * sqrt(2).
    coefficients[e * mode_count] = c * std::sqrt(2.0);
  }
  // Matching boundary data on both ends so the constant state is steady.
  const auto exterior = build_exterior(mesh, p, coefficients, c, c);

  std::vector<double> rhs(mesh.num_elements() * mode_count, 0.0);
  wasl::scalar_advection_residual(mesh, eq, p, coefficients.data(), exterior.data(), rhs.data());

  for (double value : rhs) {
    if (!approx(value, 0.0))
      return false;
  }
  return true;
}

// P0 reduction to first-order upwind finite volume: with one mode per element and
// a positive wavespeed, dc_0/dt = sqrt(2) * a (U_{e-1} - U_e) / h, so the
// underlying field rate is dU_e/dt = a (U_{e-1} - U_e) / h.
bool check_p0_upwind_finite_volume() {
  const int n = 4;
  const double length = 1.0;
  const wasl::LineMesh mesh = make_uniform_mesh(length, n);
  const wasl::ScalarAdvection eq{1.0};
  constexpr int p = 0;
  const double h = length / n;

  // Distinct cell averages.
  const double values[] = {1.0, 2.0, 4.0, 8.0};
  std::vector<double> coefficients(mesh.num_elements(), 0.0);
  for (int e = 0; e < n; ++e) {
    coefficients[e] = values[e] * std::sqrt(2.0); // mode 0 carries avg * sqrt(2)
  }
  // Left boundary feeds U = 0 (inflow); right boundary copies the last cell
  // (zero-gradient outflow, whose value is unused since wind blows outward).
  const auto exterior = build_exterior(mesh, p, coefficients, 0.0, values[n - 1]);

  std::vector<double> rhs(mesh.num_elements(), 0.0);
  wasl::scalar_advection_residual(mesh, eq, p, coefficients.data(), exterior.data(), rhs.data());

  // Expected dc_0/dt = sqrt(2) * a (U_{e-1} - U_e) / h.
  const double expected[] = {
      (0.0 - values[0]) / h,
      (values[0] - values[1]) / h,
      (values[1] - values[2]) / h,
      (values[2] - values[3]) / h,
  };
  for (int e = 0; e < n; ++e) {
    if (!approx(rhs[e], expected[e] * std::sqrt(2.0), 1e-12))
      return false;
  }
  return true;
}

// Shared-face flux cancellation: each interior face subtracts phi(+1) F* from
// its left owner and adds phi(-1) F* to its right owner. Mode zero has
// phi(+1) = phi(-1) = 1/sqrt(2), so the mode-zero contributions cancel across
// owners. With periodic boundaries the whole-domain mode-zero sum is therefore
// zero to round-off, the discrete conservation statement.
bool check_shared_face_cancellation() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 3);
  const wasl::ScalarAdvection eq{1.0};
  constexpr int p = 1;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count + 0] = 1.0; // mode 0
    coefficients[e * mode_count + 1] = 0.4; // mode 1
  }
  const auto exterior = build_periodic_exterior(mesh, p, coefficients);

  std::vector<double> rhs(mesh.num_elements() * mode_count, 0.0);
  wasl::scalar_advection_residual(mesh, eq, p, coefficients.data(), exterior.data(), rhs.data());

  double sum_mode0 = 0.0;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    sum_mode0 += rhs[e * mode_count + 0];
  }
  return approx(sum_mode0, 0.0, 1e-12);
}

// Conservation under periodic reassembly: with periodic boundaries, the total
// integral of U over the domain is constant, so the volume-weighted sum of
// mode-zero residuals must vanish. Use a nonconstant field so the check is
// nontrivial.
bool check_periodic_conservation() {
  const int n = 5;
  const double length = 2.0;
  const wasl::LineMesh mesh = make_uniform_mesh(length, n);
  const wasl::ScalarAdvection eq{0.7};
  constexpr int p = 2;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));
  const double h = length / n;

  // A smooth field via its modal moments; exact shape is irrelevant, only that
  // it varies so the residual is nontrivial.
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (int e = 0; e < n; ++e) {
    const double centre = (e + 0.5) * h;
    coefficients[e * mode_count + 0] = std::cos(centre) * std::sqrt(2.0);
    coefficients[e * mode_count + 1] = 0.0;
    coefficients[e * mode_count + 2] = 0.1 * std::sin(centre);
  }
  const auto exterior = build_periodic_exterior(mesh, p, coefficients);

  std::vector<double> rhs(mesh.num_elements() * mode_count, 0.0);
  wasl::scalar_advection_residual(mesh, eq, p, coefficients.data(), exterior.data(), rhs.data());

  // Mode 0 carries average * sqrt(2); rhs[e,0] = dc_0/dt. The integral rate
  // d/dt integral U dx scales as sum_e J_e * rhs[e,0], which conservation
  // requires to vanish (the constant sqrt(2) factor does not affect the zero).
  double dintegral = 0.0;
  for (int e = 0; e < n; ++e) {
    const wasl::LineMap map = wasl::element_map(mesh, e);
    dintegral += map.jacobian() * rhs[e * mode_count + 0];
  }
  return approx(dintegral, 0.0, 1e-10);
}

// Free-stream preservation at negative wavespeed: a constant field with matching
// periodic data must still give zero residual when wind blows from the right.
bool check_negative_wavespeed_constant_state() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  const wasl::ScalarAdvection eq{-0.9};
  constexpr int p = 1;
  const int mode_count = static_cast<int>(wasl::line_mode_count(p));

  const double c = 2.0;
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    coefficients[e * mode_count] = c * std::sqrt(2.0);
  }
  const auto exterior = build_exterior(mesh, p, coefficients, c, c);

  std::vector<double> rhs(mesh.num_elements() * mode_count, 0.0);
  wasl::scalar_advection_residual(mesh, eq, p, coefficients.data(), exterior.data(), rhs.data());

  for (double value : rhs) {
    if (!approx(value, 0.0))
      return false;
  }
  return true;
}

} // namespace

int main() {
  if (!check_constant_state_zero_residual())
    return 1;
  if (!check_p0_upwind_finite_volume())
    return 1;
  if (!check_shared_face_cancellation())
    return 1;
  if (!check_periodic_conservation())
    return 1;
  if (!check_negative_wavespeed_constant_state())
    return 1;
  return 0;
}
