#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <numbers>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/dg/euler_triangle_residual.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace {

constexpr wasl::Euler kEuler{1.4};

struct StructuredMesh {
  wasl::TriangleMesh mesh;
  int cells_per_axis;
};

StructuredMesh make_unit_square_mesh(int cells_per_axis) {
  wasl::MshDocument document;
  document.surface_entities.push_back(
      wasl::MshSurfaceEntity{1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, {1}, {}});
  for (int j = 0; j <= cells_per_axis; ++j)
    for (int i = 0; i <= cells_per_axis; ++i)
      document.nodes.push_back(wasl::MshNode{
          static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1), 2, 1,
          static_cast<double>(i) / cells_per_axis, static_cast<double>(j) / cells_per_axis, 0.0});
  auto node = [cells_per_axis](int i, int j) {
    return static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1);
  };
  std::int64_t tag = 1;
  for (int j = 0; j < cells_per_axis; ++j) {
    for (int i = 0; i < cells_per_axis; ++i) {
      const auto lower_left = node(i, j);
      const auto lower_right = node(i + 1, j);
      const auto upper_right = node(i + 1, j + 1);
      const auto upper_left = node(i, j + 1);
      document.elements.push_back(
          wasl::MshElement{tag++, 2, 1, 2, {lower_left, lower_right, upper_right}});
      document.elements.push_back(
          wasl::MshElement{tag++, 2, 1, 2, {lower_left, upper_right, upper_left}});
    }
  }
  return {wasl::build_triangle_mesh(document), cells_per_axis};
}

double component(const wasl::EulerConserved &state, int variable) {
  if (variable == 0)
    return state.density;
  if (variable == 1)
    return state.momentum_x;
  if (variable == 2)
    return state.momentum_y;
  return state.total_energy;
}

wasl::EulerConserved vortex(double x, double y, double time) {
  x -= time;
  y -= time;
  const double dx = x - 0.5;
  const double dy = y - 0.5;
  const double radius_squared = dx * dx + dy * dy;
  const double beta = 1.0;
  const double swirl = beta / (2.0 * std::numbers::pi) * std::exp(0.5 * (1.0 - radius_squared));
  const double temperature =
      1.0 - (kEuler.heat_capacity_ratio - 1.0) * beta * beta /
                (8.0 * kEuler.heat_capacity_ratio * std::numbers::pi * std::numbers::pi) *
                std::exp(1.0 - radius_squared);
  const double density = std::pow(temperature, 1.0 / (kEuler.heat_capacity_ratio - 1.0));
  return wasl::conserved(kEuler, {density, 1.0 - swirl * dy, 1.0 + swirl * dx,
                                  std::pow(density, kEuler.heat_capacity_ratio)});
}

std::vector<double> project(const wasl::TriangleMesh &mesh, int polynomial_order, double time) {
  const std::size_t element_count = mesh.num_elements();
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  const auto rule = wasl::triangle_volume_rule(5);
  std::vector<double> coefficients(wasl::Euler::variable_count * mode_count * element_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto map = wasl::triangle_element_map(mesh, element);
    for (int q = 0; q < rule.size(); ++q) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      const auto state = vortex(x, y, time);
      wasl::dubiner_basis(rule.xi[q], rule.eta[q], polynomial_order, basis.data());
      for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                       element] += rule.weights[q] * component(state, variable) * basis[mode];
    }
  }
  return coefficients;
}

std::vector<double> exact_boundary(const wasl::TriangleMesh &mesh, int polynomial_order,
                                   double time) {
  const int edge_points = polynomial_order + 1;
  const auto rules = wasl::triangle_edge_rules(polynomial_order);
  std::vector<double> exterior(wasl::Euler::variable_count * mesh.faces.size() * edge_points, 0.0);
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    if (face.right_element != wasl::invalid_element)
      continue;
    const auto map = wasl::triangle_element_map(mesh, face.left_element);
    const auto &rule = rules[static_cast<std::size_t>(face.left_local_edge)];
    for (int q = 0; q < edge_points; ++q) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      const auto state = vortex(x, y, time);
      for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
        exterior[(static_cast<std::size_t>(variable) * mesh.faces.size() + face_index) *
                     edge_points +
                 q] = component(state, variable);
    }
  }
  return exterior;
}

void residual(const wasl::TriangleMesh &mesh, int polynomial_order,
              const std::vector<double> &coefficients, double time, std::vector<double> &rhs) {
  const auto exterior = exact_boundary(mesh, polynomial_order, time);
  wasl::euler_triangle_residual(mesh, kEuler, wasl::hllc_flux, polynomial_order,
                                coefficients.data(), exterior.data(), rhs.data());
}

void advance(const StructuredMesh &structured, int polynomial_order, double target_time,
             std::vector<double> &coefficients) {
  const double requested_dt = 0.08 / structured.cells_per_axis / (3.0 * (2 * polynomial_order + 1));
  const int step_count = static_cast<int>(std::ceil(target_time / requested_dt));
  const double dt = target_time / step_count;
  std::vector<double> initial(coefficients.size()), stage(coefficients.size()),
      rhs(coefficients.size());
  double time = 0.0;
  for (int step = 0; step < step_count; ++step) {
    initial = coefficients;
    residual(structured.mesh, polynomial_order, coefficients, time, rhs);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      stage[i] = coefficients[i] + dt * rhs[i];
    residual(structured.mesh, polynomial_order, stage, time + dt, rhs);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      stage[i] = 0.75 * initial[i] + 0.25 * (stage[i] + dt * rhs[i]);
    residual(structured.mesh, polynomial_order, stage, time + 0.5 * dt, rhs);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      coefficients[i] = (1.0 / 3.0) * initial[i] + (2.0 / 3.0) * (stage[i] + dt * rhs[i]);
    time += dt;
  }
}

double density_error(const wasl::TriangleMesh &mesh, int polynomial_order,
                     const std::vector<double> &coefficients, double time) {
  const std::size_t element_count = mesh.num_elements();
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  const auto rule = wasl::triangle_volume_rule(5);
  std::vector<double> basis(mode_count);
  double squared_error = 0.0;
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto map = wasl::triangle_element_map(mesh, element);
    for (int q = 0; q < rule.size(); ++q) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(rule.xi[q], rule.eta[q], x, y);
      wasl::dubiner_basis(rule.xi[q], rule.eta[q], polynomial_order, basis.data());
      double density = 0.0;
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        density += coefficients[mode * element_count + element] * basis[mode];
      const double difference = density - vortex(x, y, time).density;
      squared_error += rule.weights[q] * map.jacobian_determinant() * difference * difference;
    }
  }
  return std::sqrt(squared_error);
}

bool check_free_stream() {
  const auto structured = make_unit_square_mesh(2);
  for (int p = 0; p <= 2; ++p) {
    const std::size_t mode_count = wasl::triangle_mode_count(p);
    const std::size_t element_count = structured.mesh.num_elements();
    const auto state = wasl::conserved(kEuler, {1.1, 0.4, -0.2, 0.9});
    std::vector<double> coefficients(wasl::Euler::variable_count * mode_count * element_count, 0.0);
    for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
      for (std::size_t element = 0; element < element_count; ++element)
        coefficients[static_cast<std::size_t>(variable) * mode_count * element_count + element] =
            component(state, variable) * std::sqrt(2.0);
    const int edge_points = p + 1;
    std::vector<double> exterior(wasl::Euler::variable_count * structured.mesh.faces.size() *
                                 edge_points);
    for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
      for (std::size_t face = 0; face < structured.mesh.faces.size(); ++face)
        for (int q = 0; q < edge_points; ++q)
          exterior[(static_cast<std::size_t>(variable) * structured.mesh.faces.size() + face) *
                       edge_points +
                   q] = component(state, variable);
    std::vector<double> rhs(coefficients.size());
    wasl::euler_triangle_residual(structured.mesh, kEuler, wasl::hllc_flux, p, coefficients.data(),
                                  exterior.data(), rhs.data());
    if (std::any_of(rhs.begin(), rhs.end(), [](double value) { return std::fabs(value) > 3e-10; }))
      return false;
  }
  return true;
}

bool check_vortex_convergence() {
  constexpr double target_time = 0.02;
  for (int p : {1, 2}) {
    double previous_error = 0.0;
    for (int n : {4, 8, 16}) {
      const auto structured = make_unit_square_mesh(n);
      auto coefficients = project(structured.mesh, p, 0.0);
      advance(structured, p, target_time, coefficients);
      const double error = density_error(structured.mesh, p, coefficients, target_time);
      if (n == 16) {
        const double rate = std::log(previous_error / error) / std::log(2.0);
        std::printf("Euler vortex P%d: rate %.2f, density L2 error %.3e\n", p, rate, error);
        if (rate < p + 0.6)
          return false;
      }
      previous_error = error;
    }
  }
  return true;
}

wasl::EulerConserved add(const wasl::EulerConserved &state, double scale,
                         const wasl::EulerConserved &change) {
  return {state.density + scale * change.density, state.momentum_x + scale * change.momentum_x,
          state.momentum_y + scale * change.momentum_y,
          state.total_energy + scale * change.total_energy};
}

bool run_shock_tube(const wasl::EulerPrimitive &left, const wasl::EulerPrimitive &right,
                    double final_time) {
  constexpr int cell_count = 200;
  const double dx = 1.0 / cell_count;
  std::vector<wasl::EulerConserved> states(cell_count), next(cell_count), fluxes(cell_count + 1);
  for (int cell = 0; cell < cell_count; ++cell)
    states[cell] = wasl::conserved(kEuler, cell < cell_count / 2 ? left : right);
  double time = 0.0;
  while (time < final_time) {
    double maximum_speed = 0.0;
    for (const auto &state : states)
      maximum_speed = std::max(maximum_speed,
                               wasl::euler_maximum_characteristic_speed(kEuler, state, 1.0, 0.0));
    const double dt = std::min(0.45 * dx / maximum_speed, final_time - time);
    fluxes.front() = wasl::euler_flux_x(kEuler, states.front());
    for (int face = 1; face < cell_count; ++face)
      fluxes[face] = wasl::hllc_flux(kEuler, states[face - 1], states[face], 1.0, 0.0);
    fluxes.back() = wasl::euler_flux_x(kEuler, states.back());
    for (int cell = 0; cell < cell_count; ++cell) {
      const wasl::EulerConserved divergence{fluxes[cell + 1].density - fluxes[cell].density,
                                            fluxes[cell + 1].momentum_x - fluxes[cell].momentum_x,
                                            fluxes[cell + 1].momentum_y - fluxes[cell].momentum_y,
                                            fluxes[cell + 1].total_energy -
                                                fluxes[cell].total_energy};
      next[cell] = add(states[cell], -dt / dx, divergence);
      if (!wasl::is_physical(kEuler, next[cell]))
        return false;
    }
    states.swap(next);
    time += dt;
  }
  const auto [minimum, maximum] =
      std::minmax_element(states.begin(), states.end(),
                          [](const auto &a, const auto &b) { return a.density < b.density; });
  std::printf("Euler shock tube: density range [%.3f, %.3f]\n", minimum->density, maximum->density);
  return maximum->density - minimum->density > 0.2;
}

bool check_shock_tubes() {
  return run_shock_tube({1.0, 0.0, 0.0, 1.0}, {0.125, 0.0, 0.0, 0.1}, 0.2) &&
         run_shock_tube({0.445, 0.698, 0.0, 3.528}, {0.5, 0.0, 0.0, 0.571}, 0.14);
}

} // namespace

int main() {
  return check_free_stream() && check_shock_tubes() && check_vortex_convergence() ? 0 : 1;
}
