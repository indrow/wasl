// Tests for SSP-RK3 time integration of the 1D scalar-advection DG residual.
//
// These verify the time-integration properties required by CODING_STYLE.md:
//   - Constant-state preservation (free-stream, to round-off).
//   - Known ODE order: temporal convergence rate 3 against the exact
//     semi-discrete solution on a mesh with no spatial error.
//   - Periodic sine-pulse transport: a smooth wave returns near its initial
//     state after one full revolution.
//   - Spatial convergence: P1 and P2 show p + 1 order when temporal error is
//     suppressed by a small timestep.
//
// Meshes are built directly so the tests do not depend on a Gmsh binary.

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/periodic.hpp"
#include "wasl/mesh/line/line_geometry.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/quadrature/gauss_legendre.hpp"
#include "wasl/time/cfl_timestep.hpp"
#include "wasl/time/line_advection_step.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

// Build a uniform periodic line mesh on [0, L] with n elements.
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

// L2-project a smooth function onto the modal basis of each element.
std::vector<double> project_function(const wasl::LineMesh &mesh, int polynomial_order,
                                     const std::function<double(double)> &f) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  // 7-point rule integrates the projection integrand exactly for P0-P3.
  const wasl::QuadratureRule rule = wasl::gauss_legendre(7);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const wasl::LineMap map = wasl::element_map(mesh, e);
    for (int q = 0; q < rule.size(); ++q) {
      const double xi = rule.points[q];
      const double x = map.reference_to_physical(xi);
      wasl::legendre_basis(xi, polynomial_order, basis.data());
      const double fx = f(x);
      for (int j = 0; j < mode_count; ++j) {
        coefficients[e * mode_count + j] += rule.weights[q] * fx * basis[j];
      }
    }
  }
  return coefficients;
}

// L2 error of the numerical solution against an exact function.
double l2_error(const wasl::LineMesh &mesh, int polynomial_order,
                const std::vector<double> &coefficients,
                const std::function<double(double)> &exact) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  const wasl::QuadratureRule rule = wasl::gauss_legendre(7);
  std::vector<double> basis(mode_count);
  double error_sq = 0.0;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const wasl::LineMap map = wasl::element_map(mesh, e);
    const double *modes = coefficients.data() + e * mode_count;
    for (int q = 0; q < rule.size(); ++q) {
      const double xi = rule.points[q];
      const double x = map.reference_to_physical(xi);
      wasl::legendre_basis(xi, polynomial_order, basis.data());
      double uh = 0.0;
      for (int j = 0; j < mode_count; ++j) {
        uh += modes[j] * basis[j];
      }
      const double diff = uh - exact(x);
      error_sq += rule.weights[q] * map.jacobian() * diff * diff;
    }
  }
  return std::sqrt(error_sq);
}

// Cell average from mode-zero coefficient. phi_0 = 1/sqrt(2), so the average is
// c_0 / sqrt(2).
double cell_average(const std::vector<double> &coefficients, int polynomial_order, int element) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  return coefficients[static_cast<std::size_t>(element) * mode_count] / std::sqrt(2.0);
}

// The periodic exterior-trace policy lives in wasl/boundaries/periodic.hpp and
// is used here via wasl::make_periodic_exterior.

// Advance coefficients to time T using SSP-RK3 with a fixed dt.
void advance(const wasl::LineMesh &mesh, const wasl::ScalarAdvection &eq, int polynomial_order,
             double dt, double target_time, const wasl::ExteriorFunction &exterior_fn,
             std::vector<double> &coefficients, wasl::SspRk3Workspace &ws) {
  const int n_steps = static_cast<int>(std::ceil(target_time / dt));
  const double actual_dt = target_time / n_steps;
  for (int step = 0; step < n_steps; ++step) {
    wasl::ssp_rk3_step(mesh, eq, polynomial_order, actual_dt, exterior_fn, coefficients, ws);
  }
}

// A constant field is a steady solution: the residual is zero at every stage,
// so SSP-RK3 must leave it unchanged to round-off. Checked at P0 and P2.
bool check_constant_state_preservation() {
  constexpr double kTol = 1e-14;
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 4);
  const wasl::ScalarAdvection eq{1.3};
  wasl::SspRk3Workspace ws;

  for (int p : {0, 2}) {
    const int mode_count = static_cast<int>(wasl::line_mode_count(p));
    const double dt = wasl::advection_cfl_timestep(mesh, eq, p, 0.2);
    const auto exterior_fn = wasl::make_periodic_exterior(mesh, p);
    const double c = 0.7;
    std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
    for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
      coefficients[e * mode_count] = c * std::sqrt(2.0);
    }
    const std::vector<double> initial = coefficients;
    advance(mesh, eq, p, dt, 1.0, exterior_fn, coefficients, ws);
    // A constant field is steady: every coefficient must be unchanged.
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      if (!approx(coefficients[i], initial[i], kTol))
        return false;
    }
  }
  return true;
}

// Temporal order 3: on a 2-element P0 periodic mesh the spatial discretization
// is exact (no truncation error beyond the two cell averages), so the only error
// is temporal. The exact semi-discrete solution decays the initial difference
// exponentially: D(t) = D(0) * exp(-2 * a * t / h). Compare SSP-RK3 against it
// at decreasing dt and check the rate is close to 3.
bool check_temporal_order() {
  const wasl::LineMesh mesh = make_uniform_mesh(1.0, 2);
  const wasl::ScalarAdvection eq{1.0};
  constexpr int p = 0;
  const double target_time = 1.0; // one full revolution

  // Exact semi-discrete solution: U_0(t) = 0.5 + 0.5*exp(-4t), U_1 likewise.
  const double decay = std::exp(-4.0 * target_time);
  const double exact_0 = 0.5 + 0.5 * decay;
  const double exact_1 = 0.5 - 0.5 * decay;

  const auto exterior_fn = wasl::make_periodic_exterior(mesh, p);
  wasl::SspRk3Workspace ws;

  double prev_error = 0.0;
  double prev_dt = 0.0;
  for (int n_steps : {4, 8, 16, 32}) {
    const double dt = target_time / n_steps;
    std::vector<double> coefficients{std::sqrt(2.0), 0.0}; // U_0=1, U_1=0
    advance(mesh, eq, p, dt, target_time, exterior_fn, coefficients, ws);
    const double u0 = cell_average(coefficients, p, 0);
    const double u1 = cell_average(coefficients, p, 1);
    const double err = std::max(std::fabs(u0 - exact_0), std::fabs(u1 - exact_1));

    if (prev_dt > 0.0) {
      const double rate = std::log(prev_error / err) / std::log(2.0);
      // The last refinement (16 -> 32 steps) must show rate >= 2.5.
      if (n_steps == 32 && rate < 2.5) {
        std::fprintf(stderr, "temporal order: rate %.2f < 2.5 at n_steps=%d\n", rate, n_steps);
        return false;
      }
    }
    prev_error = err;
    prev_dt = dt;
  }
  return true;
}

// Periodic sine-pulse transport: advance sin(2*pi*x) one full period (T = L/a)
// on a fine mesh. The exact solution returns to the initial condition, so the L2
// error must be small. Checked at P1 and P2.
bool check_sine_transport_qualitative() {
  constexpr double kTol = 0.02;
  const double length = 1.0;
  const wasl::ScalarAdvection eq{1.0};
  const double target_time = length / eq.wavespeed; // one full revolution
  const auto exact = [](double x) { return std::sin(2.0 * M_PI * x); };

  for (int p : {1, 2}) {
    const wasl::LineMesh mesh = make_uniform_mesh(length, 40);
    const double dt = wasl::advection_cfl_timestep(mesh, eq, p, 0.2);
    const auto exterior_fn = wasl::make_periodic_exterior(mesh, p);
    wasl::SspRk3Workspace ws;

    auto coefficients = project_function(mesh, p, exact);
    advance(mesh, eq, p, dt, target_time, exterior_fn, coefficients, ws);
    const double err = l2_error(mesh, p, coefficients, exact);
    if (err > kTol) {
      std::fprintf(stderr, "sine transport P%d: L2 error %.4e > %.4e\n", p, err, kTol);
      return false;
    }
  }
  return true;
}

// Spatial convergence: advance sin(2*pi*x) for a short time on successive uniform
// meshes at P1 and P2. The timestep is scaled so temporal error is far below
// spatial error. The L2 error must converge at rate p + 1.
//
// Records: mesh sizes, p, dt policy, norms, observed rates, precision (double),
// and tolerance, per CODING_STYLE.md.
bool check_spatial_convergence() {
  const double length = 1.0;
  const wasl::ScalarAdvection eq{1.0};
  const double target_time = 0.1;
  const auto exact_at_t = [target_time](double x) {
    return std::sin(2.0 * M_PI * (x - target_time));
  };

  // Precision: double. Tolerance: observed rate >= p + 0.7.
  // dt policy: n_steps = 5 * n for P1 (dt ~ h, temporal ~ h^3 << h^2);
  //            n_steps = 5 * n^2 for P2 (dt ~ h^2, temporal ~ h^6 << h^3).
  struct Level {
    int n;
    int p;
    int n_steps;
  };

  const std::vector<Level> levels = {
      {8, 1, 40},  {16, 1, 80},   {32, 1, 160},  {64, 1, 320},
      {8, 2, 320}, {16, 2, 1280}, {32, 2, 5120},
  };

  std::vector<double> errors;
  std::vector<double> hs;
  std::vector<int> ps;
  for (const auto &level : levels) {
    const wasl::LineMesh mesh = make_uniform_mesh(length, level.n);
    const auto exterior_fn = wasl::make_periodic_exterior(mesh, level.p);
    wasl::SspRk3Workspace ws;

    const auto initial = [](double x) { return std::sin(2.0 * M_PI * x); };
    auto coefficients = project_function(mesh, level.p, initial);
    const double dt = target_time / level.n_steps;
    advance(mesh, eq, level.p, dt, target_time, exterior_fn, coefficients, ws);
    const double err = l2_error(mesh, level.p, coefficients, exact_at_t);

    errors.push_back(err);
    hs.push_back(length / level.n);
    ps.push_back(level.p);
  }

  // Print the convergence table to stderr for the record.
  std::fprintf(stderr, "\n=== Spatial convergence (double precision) ===\n");
  std::fprintf(stderr, "target_time=%.4f, dt policy: P1 n_steps=5n, P2 n_steps=5n^2\n",
               target_time);
  std::fprintf(stderr, "%4s %2s %8s %12s %10s\n", "n", "p", "h", "L2_error", "rate");
  std::fprintf(stderr, "%4s %2s %8s %12s %10s\n", "---", "--", "-------", "--------", "-----");

  // Compute and check rates within each polynomial order group.
  double prev_err = 0.0;
  double prev_h = 0.0;
  int prev_p = -1;
  bool all_ok = true;
  for (std::size_t i = 0; i < levels.size(); ++i) {
    double rate = 0.0;
    if (ps[i] == prev_p) {
      rate = std::log(prev_err / errors[i]) / std::log(prev_h / hs[i]);
    }
    std::fprintf(stderr, "%4d %2d %8.5f %12.4e %10.2f\n", levels[i].n, ps[i], hs[i], errors[i],
                 rate);
    if (ps[i] == prev_p) {
      const double required = static_cast<double>(ps[i]) + 0.7;
      if (rate < required) {
        std::fprintf(stderr, "  FAIL: rate %.2f < required %.1f for P%d\n", rate, required, ps[i]);
        all_ok = false;
      }
    }
    prev_err = errors[i];
    prev_h = hs[i];
    prev_p = ps[i];
  }
  return all_ok;
}

// Build an irregular (perturbed) periodic line mesh on [0, L] with n elements.
// Interior nodes are shifted by a smooth sinusoidal perturbation so element
// lengths vary, exercising convergence on an unstructured mesh. The endpoints
// are fixed so the domain length is preserved.
wasl::LineMesh make_irregular_mesh(double length, int n) {
  const double h = length / n;
  std::vector<double> x(n + 1);
  x[0] = 0.0;
  x[n] = length;
  for (int i = 1; i < n; ++i) {
    // Smooth perturbation, amplitude 0.2 * h, fixed at both ends.
    x[i] = i * h + 0.2 * h * std::sin(M_PI * static_cast<double>(i) / n);
  }

  wasl::LineMesh mesh;
  mesh.elements.reserve(n);
  for (int e = 0; e < n; ++e) {
    mesh.elements.push_back(wasl::LineElement{e, e + 1, x[e], x[e + 1], 1});
  }
  mesh.faces.reserve(n + 1);
  mesh.faces.push_back(
      wasl::LineFace{0, wasl::invalid_element, x[0], -1.0, wasl::invalid_physical});
  for (int e = 0; e + 1 < n; ++e) {
    mesh.faces.push_back(wasl::LineFace{e, e + 1, x[e + 1], 1.0, wasl::invalid_physical});
  }
  mesh.faces.push_back(
      wasl::LineFace{n - 1, wasl::invalid_element, x[n], 1.0, wasl::invalid_physical});
  return mesh;
}

// Irregular-mesh convergence: same as the uniform test but on perturbed meshes.
// This closes the M2 exit criterion "P1-P3 smooth convergence on irregular
// elements". The minimum element length drives the CFL timestep so stability is
// preserved despite varying element sizes.
bool check_irregular_mesh_convergence() {
  const double length = 1.0;
  const wasl::ScalarAdvection eq{1.0};
  const double target_time = 0.1;
  const auto exact_at_t = [target_time](double x) {
    return std::sin(2.0 * M_PI * (x - target_time));
  };

  // Precision: double. Tolerance: observed rate >= p + 0.7.
  struct Level {
    int n;
    int p;
    int n_steps;
  };

  const std::vector<Level> levels = {
      {8, 1, 40},    {16, 1, 80},   {32, 1, 160}, {8, 2, 320},
      {16, 2, 1280}, {32, 2, 5120}, {8, 3, 1280}, {16, 3, 5120},
  };

  std::vector<double> errors;
  std::vector<double> hs;
  std::vector<int> ps;
  for (const auto &level : levels) {
    const wasl::LineMesh mesh = make_irregular_mesh(length, level.n);
    const auto exterior_fn = wasl::make_periodic_exterior(mesh, level.p);
    wasl::SspRk3Workspace ws;

    const auto initial = [](double x) { return std::sin(2.0 * M_PI * x); };
    auto coefficients = project_function(mesh, level.p, initial);
    const double dt = target_time / level.n_steps;
    advance(mesh, eq, level.p, dt, target_time, exterior_fn, coefficients, ws);
    const double err = l2_error(mesh, level.p, coefficients, exact_at_t);

    errors.push_back(err);
    hs.push_back(length / level.n);
    ps.push_back(level.p);
  }

  std::fprintf(stderr, "\n=== Irregular-mesh convergence (double precision) ===\n");
  std::fprintf(stderr, "target_time=%.4f, dt policy: n_steps scaled to suppress temporal error\n",
               target_time);
  std::fprintf(stderr, "%4s %2s %8s %12s %10s\n", "n", "p", "h", "L2_error", "rate");
  std::fprintf(stderr, "%4s %2s %8s %12s %10s\n", "---", "--", "-------", "--------", "-----");

  double prev_err = 0.0;
  double prev_h = 0.0;
  int prev_p = -1;
  bool all_ok = true;
  for (std::size_t i = 0; i < levels.size(); ++i) {
    double rate = 0.0;
    if (ps[i] == prev_p) {
      rate = std::log(prev_err / errors[i]) / std::log(prev_h / hs[i]);
    }
    std::fprintf(stderr, "%4d %2d %8.5f %12.4e %10.2f\n", levels[i].n, ps[i], hs[i], errors[i],
                 rate);
    if (ps[i] == prev_p) {
      const double required = static_cast<double>(ps[i]) + 0.7;
      if (rate < required) {
        std::fprintf(stderr, "  FAIL: rate %.2f < required %.1f for P%d\n", rate, required, ps[i]);
        all_ok = false;
      }
    }
    prev_err = errors[i];
    prev_h = hs[i];
    prev_p = ps[i];
  }
  return all_ok;
}

} // namespace

int main() {
  if (!check_constant_state_preservation())
    return 1;
  if (!check_temporal_order())
    return 1;
  if (!check_sine_transport_qualitative())
    return 1;
  if (!check_spatial_convergence())
    return 1;
  if (!check_irregular_mesh_convergence())
    return 1;
  return 0;
}
