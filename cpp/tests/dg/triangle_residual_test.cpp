// Conservation, constant-state, and smooth-convergence checks for the CPU
// scalar-advection residual on affine triangles.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <numbers>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/dg/triangle_residual.hpp"
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace {

struct StructuredMesh {
  wasl::TriangleMesh mesh;
  int cells_per_axis;
};

StructuredMesh make_unit_square_mesh(int cells_per_axis) {
  wasl::MshDocument document;
  document.surface_entities.push_back(
      wasl::MshSurfaceEntity{1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, {1}, {}});
  for (int j = 0; j <= cells_per_axis; ++j) {
    for (int i = 0; i <= cells_per_axis; ++i) {
      const auto tag = static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1);
      document.nodes.push_back(wasl::MshNode{tag, 2, 1, static_cast<double>(i) / cells_per_axis,
                                             static_cast<double>(j) / cells_per_axis, 0.0});
    }
  }
  auto node = [cells_per_axis](int i, int j) {
    return static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1);
  };
  std::int64_t element_tag = 1;
  for (int j = 0; j < cells_per_axis; ++j) {
    for (int i = 0; i < cells_per_axis; ++i) {
      const auto lower_left = node(i, j);
      const auto lower_right = node(i + 1, j);
      const auto upper_right = node(i + 1, j + 1);
      const auto upper_left = node(i, j + 1);
      document.elements.push_back(
          wasl::MshElement{element_tag++, 2, 1, 2, {lower_left, lower_right, upper_right}});
      document.elements.push_back(
          wasl::MshElement{element_tag++, 2, 1, 2, {lower_left, upper_right, upper_left}});
    }
  }
  return {wasl::build_triangle_mesh(document), cells_per_axis};
}

double evaluate(const double *modes, int polynomial_order, double xi, double eta) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  std::vector<double> basis(mode_count);
  wasl::dubiner_basis(xi, eta, polynomial_order, basis.data());
  double state = 0.0;
  for (std::size_t mode = 0; mode < mode_count; ++mode)
    state += modes[mode] * basis[mode];
  return state;
}

std::vector<double> boundary_constant(const wasl::TriangleMesh &mesh, int polynomial_order,
                                      double value) {
  return std::vector<double>(mesh.faces.size() * static_cast<std::size_t>(polynomial_order + 1),
                             value);
}

bool check_constant_state() {
  const auto structured = make_unit_square_mesh(2);
  const wasl::ScalarAdvection2D equation{0.7, -1.1};
  for (int p = 0; p <= 2; ++p) {
    const std::size_t mode_count = wasl::triangle_mode_count(p);
    const double state = 0.8;
    std::vector<double> coefficients(structured.mesh.num_elements() * mode_count, 0.0);
    for (std::size_t element = 0; element < structured.mesh.num_elements(); ++element)
      coefficients[element * mode_count] = state * std::sqrt(2.0);
    const auto exterior = boundary_constant(structured.mesh, p, state);
    std::vector<double> rhs(coefficients.size());
    wasl::triangle_scalar_advection_residual(structured.mesh, equation, p, coefficients.data(),
                                             exterior.data(), rhs.data());
    for (double value : rhs) {
      if (std::fabs(value) > 2e-11)
        return false;
    }
  }
  return true;
}

bool check_shared_face_conservation(int p) {
  const auto structured = make_unit_square_mesh(3);
  const wasl::ScalarAdvection2D equation{0.9, -0.4};
  const std::size_t mode_count = wasl::triangle_mode_count(p);
  std::vector<double> coefficients(structured.mesh.num_elements() * mode_count);
  for (std::size_t i = 0; i < coefficients.size(); ++i)
    coefficients[i] = 0.2 * std::sin(0.7 * static_cast<double>(i + 1));
  std::vector<double> exterior(structured.mesh.faces.size() * (p + 1), 0.0);
  for (std::size_t i = 0; i < exterior.size(); ++i)
    exterior[i] = 0.1 * std::cos(0.3 * static_cast<double>(i + 1));

  std::vector<double> rhs(coefficients.size());
  wasl::triangle_scalar_advection_residual(structured.mesh, equation, p, coefficients.data(),
                                           exterior.data(), rhs.data());
  double integral_rate = 0.0;
  for (std::size_t element = 0; element < structured.mesh.num_elements(); ++element) {
    integral_rate += std::sqrt(2.0) *
                     wasl::triangle_element_map(structured.mesh, element).jacobian_determinant() *
                     rhs[element * mode_count];
  }

  const auto edge_rules = wasl::triangle_edge_rules(p);
  double boundary_flux = 0.0;
  for (std::size_t face_index = 0; face_index < structured.mesh.faces.size(); ++face_index) {
    const auto &face = structured.mesh.faces[face_index];
    if (face.right_element != wasl::invalid_element)
      continue;
    const auto map = wasl::triangle_element_map(structured.mesh, face.left_element);
    const auto &rule = edge_rules[static_cast<std::size_t>(face.left_local_edge)];
    const double reference_length = face.left_local_edge == 1 ? 2.0 * std::sqrt(2.0) : 2.0;
    const double scale = map.edge_length(face.left_local_edge) / reference_length;
    const double *modes =
        coefficients.data() + static_cast<std::size_t>(face.left_element) * mode_count;
    for (int q = 0; q <= p; ++q) {
      const double interior = evaluate(modes, p, rule.xi[q], rule.eta[q]);
      const double flux =
          wasl::upwind_normal_flux(equation, face.outward_normal_x, face.outward_normal_y, interior,
                                   exterior[face_index * (p + 1) + q]);
      boundary_flux += rule.weights[q] * scale * flux;
    }
  }
  return std::fabs(integral_rate + boundary_flux) < 2e-12;
}

std::vector<double> project(const wasl::TriangleMesh &mesh, int polynomial_order,
                            const std::function<double(double, double)> &function) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  const auto rule = wasl::triangle_volume_rule(5);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = wasl::triangle_element_map(mesh, element);
    for (int q = 0; q < rule.size(); ++q) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      wasl::dubiner_basis(rule.xi[q], rule.eta[q], polynomial_order, basis.data());
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        coefficients[element * mode_count + mode] += rule.weights[q] * function(x, y) * basis[mode];
    }
  }
  return coefficients;
}

double periodic_trace(const StructuredMesh &structured, int polynomial_order,
                      const std::vector<double> &coefficients, double x, double y) {
  const int n = structured.cells_per_axis;
  const int i = std::min(static_cast<int>(x * n), n - 1);
  const int j = std::min(static_cast<int>(y * n), n - 1);
  const double local_x = x * n - i;
  const double local_y = y * n - j;
  const int triangle = local_y <= local_x ? 0 : 1;
  const std::size_t element = static_cast<std::size_t>(2 * (j * n + i) + triangle);
  const auto map = wasl::triangle_element_map(structured.mesh, element);
  double xi = 0.0, eta = 0.0;
  map.physical_to_reference(x, y, xi, eta);
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  return evaluate(coefficients.data() + element * mode_count, polynomial_order, xi, eta);
}

std::vector<double> periodic_exterior(const StructuredMesh &structured, int polynomial_order,
                                      const std::vector<double> &coefficients) {
  const auto rules = wasl::triangle_edge_rules(polynomial_order);
  const int edge_points = polynomial_order + 1;
  std::vector<double> exterior(structured.mesh.faces.size() * edge_points, 0.0);
  for (std::size_t face_index = 0; face_index < structured.mesh.faces.size(); ++face_index) {
    const auto &face = structured.mesh.faces[face_index];
    if (face.right_element != wasl::invalid_element)
      continue;
    const auto map = wasl::triangle_element_map(structured.mesh, face.left_element);
    const auto &rule = rules[static_cast<std::size_t>(face.left_local_edge)];
    for (int q = 0; q < edge_points; ++q) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      x -= face.outward_normal_x;
      y -= face.outward_normal_y;
      exterior[face_index * edge_points + q] =
          periodic_trace(structured, polynomial_order, coefficients, x, y);
    }
  }
  return exterior;
}

template <typename VelocityFunction>
void advance(const StructuredMesh &structured, VelocityFunction velocity, double maximum_speed,
             int polynomial_order, double target_time, std::vector<double> &coefficients) {
  const double h = 1.0 / structured.cells_per_axis;
  const double requested_dt = 0.08 * h / (maximum_speed * (2 * polynomial_order + 1));
  const int steps = static_cast<int>(std::ceil(target_time / requested_dt));
  const double dt = target_time / steps;
  std::vector<double> initial(coefficients.size()), stage(coefficients.size()),
      rhs(coefficients.size());
  initial = coefficients;
  for (int step = 0; step < steps; ++step) {
    initial = coefficients;
    auto exterior = periodic_exterior(structured, polynomial_order, coefficients);
    wasl::triangle_scalar_advection_residual(structured.mesh, velocity, polynomial_order,
                                             coefficients.data(), exterior.data(), rhs.data());
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      stage[i] = coefficients[i] + dt * rhs[i];

    exterior = periodic_exterior(structured, polynomial_order, stage);
    wasl::triangle_scalar_advection_residual(structured.mesh, velocity, polynomial_order,
                                             stage.data(), exterior.data(), rhs.data());
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      stage[i] = 0.75 * initial[i] + 0.25 * (stage[i] + dt * rhs[i]);

    exterior = periodic_exterior(structured, polynomial_order, stage);
    wasl::triangle_scalar_advection_residual(structured.mesh, velocity, polynomial_order,
                                             stage.data(), exterior.data(), rhs.data());
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      coefficients[i] = (1.0 / 3.0) * initial[i] + (2.0 / 3.0) * (stage[i] + dt * rhs[i]);
  }
}

double l2_error(const wasl::TriangleMesh &mesh, int polynomial_order,
                const std::vector<double> &coefficients,
                const std::function<double(double, double)> &exact) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  const auto rule = wasl::triangle_volume_rule(5);
  std::vector<double> basis(mode_count);
  double error_squared = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = wasl::triangle_element_map(mesh, element);
    for (int q = 0; q < rule.size(); ++q) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      const double difference = evaluate(coefficients.data() + element * mode_count,
                                         polynomial_order, rule.xi[q], rule.eta[q]) -
                                exact(x, y);
      error_squared += rule.weights[q] * map.jacobian_determinant() * difference * difference;
    }
  }
  return std::sqrt(error_squared);
}

bool check_smooth_convergence() {
  const wasl::ScalarAdvection2D equation{0.8, -0.35};
  constexpr double target_time = 0.05;
  const auto initial = [](double x, double y) {
    return std::sin(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y);
  };
  const auto exact = [&](double x, double y) {
    return initial(x - equation.velocity_x * target_time, y - equation.velocity_y * target_time);
  };
  for (int p : {1, 2}) {
    double previous_error = 0.0;
    for (int n : {8, 16, 32}) {
      const auto structured = make_unit_square_mesh(n);
      auto coefficients = project(structured.mesh, p, initial);
      advance(
          structured, [equation](double, double) { return equation; },
          std::hypot(equation.velocity_x, equation.velocity_y), p, target_time, coefficients);
      const double error = l2_error(structured.mesh, p, coefficients, exact);
      if (n == 32) {
        const double rate = std::log(previous_error / error) / std::log(2.0);
        std::printf("triangle P%d convergence: rate %.2f, L2 error %.3e\n", p, rate, error);
        if (rate < p + 0.7) {
          std::fprintf(stderr, "triangle P%d convergence rate %.2f, error %.3e\n", p, rate, error);
          return false;
        }
      }
      previous_error = error;
    }
  }
  return true;
}

double total_mass(const wasl::TriangleMesh &mesh, int polynomial_order,
                  const std::vector<double> &coefficients) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  double mass = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    mass += std::sqrt(2.0) * wasl::triangle_element_map(mesh, element).jacobian_determinant() *
            coefficients[element * mode_count];
  }
  return mass;
}

bool check_rotating_cone() {
  constexpr int p = 2;
  constexpr int cells_per_axis = 8;
  constexpr double radius = 0.15;
  constexpr double angular_speed = 2.0 * std::numbers::pi;
  const auto structured = make_unit_square_mesh(cells_per_axis);
  const auto cone = [](double x, double y) {
    const double distance = std::hypot(x - 0.25, y - 0.5);
    return std::max(0.0, 1.0 - distance / radius);
  };
  const auto rotation = [](double x, double y) {
    return wasl::ScalarAdvection2D{-angular_speed * (y - 0.5), angular_speed * (x - 0.5)};
  };
  auto coefficients = project(structured.mesh, p, cone);
  const double initial_mass = total_mass(structured.mesh, p, coefficients);
  advance(structured, rotation, angular_speed / std::sqrt(2.0), p, 1.0, coefficients);
  const double mass_drift = std::fabs(total_mass(structured.mesh, p, coefficients) - initial_mass);
  const double error = l2_error(structured.mesh, p, coefficients, cone);
  std::printf("triangle rotating cone P%d: mass drift %.3e, L2 error %.3e\n", p, mass_drift, error);
  return mass_drift < 2e-12 && error < 0.03;
}

} // namespace

int main() {
  if (!check_constant_state())
    return 1;
  for (int p = 0; p <= 2; ++p) {
    if (!check_shared_face_conservation(p))
      return 1;
  }
  if (!check_smooth_convergence())
    return 1;
  if (!check_rotating_cone())
    return 1;
  return 0;
}
