#include "wasl/case/periodic_boundary.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/case/coefficient_layout.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace wasl {
namespace {

double euler_maximum_speed(const Euler &equation, const std::vector<double> &coefficients,
                           std::size_t element_count, std::size_t mode_count, int dimensions) {
  double result = 0.0;
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto state = average_euler_state(coefficients, element_count, mode_count, element);
    if (!is_physical(equation, state))
      throw std::runtime_error("simulation: Euler cell average is not physical");
    const auto value = primitive(equation, state);
    const auto sound_speed =
        std::sqrt(equation.heat_capacity_ratio * value.pressure / value.density);
    const auto velocity = dimensions == 1 ? std::fabs(value.velocity_x)
                                          : std::hypot(value.velocity_x, value.velocity_y);
    result = std::max(result, velocity + sound_speed);
  }
  return result;
}

} // namespace

PeriodicTriangleTrace periodic_triangle_trace(const TriangleMesh &mesh, int polynomial_order) {
  const auto rules = triangle_edge_rules(polynomial_order);
  const auto mode_count = triangle_mode_count(polynomial_order);
  const int point_count = polynomial_order + 1;
  PeriodicTriangleTrace result;
  result.elements.assign(mesh.faces.size() * point_count, invalid_element);
  result.basis.assign(mesh.faces.size() * point_count * mode_count, 0.0);
  result.point_count = point_count;
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();
  for (const auto &element : mesh.elements) {
    min_x = std::min({min_x, element.x0, element.x1, element.x2});
    max_x = std::max({max_x, element.x0, element.x1, element.x2});
    min_y = std::min({min_y, element.y0, element.y1, element.y2});
    max_y = std::max({max_y, element.y0, element.y1, element.y2});
  }
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    if (face.right_element != invalid_element) {
      continue;
    }
    const auto map = triangle_element_map(mesh, face.left_element);
    const auto &rule = rules.at(static_cast<std::size_t>(face.left_local_edge));
    for (int point = 0; point < point_count; ++point) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(rule.xi[point], rule.eta[point], x, y);
      if (std::fabs(face.outward_normal_x) > 0.5 && std::fabs(face.outward_normal_y) < 1e-10) {
        x += face.outward_normal_x < 0.0 ? max_x - min_x : min_x - max_x;
      } else if (std::fabs(face.outward_normal_y) > 0.5 &&
                 std::fabs(face.outward_normal_x) < 1e-10) {
        y += face.outward_normal_y < 0.0 ? max_y - min_y : min_y - max_y;
      } else {
        throw std::runtime_error("simulation: 2D periodic boundaries require an axis-aligned box");
      }
      const auto index = face_index * point_count + static_cast<std::size_t>(point);
      for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
        double xi = 0.0;
        double eta = 0.0;
        triangle_element_map(mesh, element).physical_to_reference(x, y, xi, eta);
        if (xi >= -1.0 - 1e-10 && eta >= -1.0 - 1e-10 && xi + eta <= 1e-10) {
          result.elements[index] = static_cast<std::int32_t>(element);
          dubiner_basis(xi, eta, polynomial_order, result.basis.data() + index * mode_count);
          break;
        }
      }
      if (result.elements[index] == invalid_element) {
        throw std::runtime_error("simulation: periodic 2D trace has no donor element");
      }
    }
  }
  return result;
}

std::vector<double> periodic_triangle_exterior(const PeriodicTriangleTrace &trace,
                                               const std::vector<double> &coefficients,
                                               std::size_t mode_count) {
  std::vector<double> exterior(trace.elements.size(), 0.0);
  for (std::size_t point = 0; point < trace.elements.size(); ++point) {
    if (trace.elements[point] == invalid_element) {
      continue;
    }
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      exterior[point] +=
          coefficients[static_cast<std::size_t>(trace.elements[point]) * mode_count + mode] *
          trace.basis[point * mode_count + mode];
    }
  }
  return exterior;
}

std::vector<double> periodic_line_euler_exterior(const LineMesh &mesh, int polynomial_order,
                                                 const std::vector<double> &coefficients) {
  const auto mode_count = line_mode_count(polynomial_order);
  const auto element_count = mesh.num_elements();
  std::vector<double> left_basis(mode_count), right_basis(mode_count);
  legendre_basis(-1.0, polynomial_order, left_basis.data());
  legendre_basis(1.0, polynomial_order, right_basis.data());
  std::vector<double> exterior(Euler::variable_count * mesh.faces.size(), 0.0);
  for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
    if (mesh.faces[face].right_element != invalid_element)
      continue;
    const bool left_boundary = mesh.faces[face].outward_normal < 0.0;
    const auto donor = left_boundary ? element_count - 1 : 0;
    const auto &basis = left_boundary ? right_basis : left_basis;
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        exterior[static_cast<std::size_t>(variable) * mesh.faces.size() + face] +=
            coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                         donor] *
            basis[mode];
  }
  return exterior;
}

std::vector<double> periodic_triangle_euler_exterior(const PeriodicTriangleTrace &trace,
                                                     const std::vector<double> &coefficients,
                                                     std::size_t element_count,
                                                     std::size_t mode_count,
                                                     std::size_t face_count) {
  std::vector<double> exterior(Euler::variable_count * trace.elements.size(), 0.0);
  for (std::size_t point = 0; point < trace.elements.size(); ++point) {
    if (trace.elements[point] == invalid_element)
      continue;
    const auto face = point / static_cast<std::size_t>(trace.point_count);
    const auto quadrature_point = point % static_cast<std::size_t>(trace.point_count);
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        exterior[(static_cast<std::size_t>(variable) * face_count + face) * trace.point_count +
                 quadrature_point] +=
            coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                         static_cast<std::size_t>(trace.elements[point])] *
            trace.basis[point * mode_count + mode];
  }
  return exterior;
}

double triangle_scalar_cfl_timestep(const TriangleMesh &mesh, double velocity_x, double velocity_y,
                                    int polynomial_order, double cfl) {
  const auto speed = std::hypot(velocity_x, velocity_y);
  if (mesh.elements.empty() || speed == 0.0) {
    return 0.0;
  }
  double min_altitude = std::numeric_limits<double>::max();
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = triangle_element_map(mesh, element);
    const auto longest_edge =
        std::max({map.edge_length(0), map.edge_length(1), map.edge_length(2)});
    min_altitude = std::min(min_altitude, 4.0 * map.jacobian_determinant() / longest_edge);
  }
  return cfl * min_altitude / (speed * (2.0 * polynomial_order + 1.0));
}

double euler_cfl_timestep(const LineMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl) {
  double minimum_length = std::numeric_limits<double>::max();
  for (const auto &element : mesh.elements)
    minimum_length = std::min(minimum_length, element.x_right - element.x_left);
  const auto speed =
      euler_maximum_speed(equation, coefficients, mesh.num_elements(), mode_count, 1);
  return speed > 0.0 ? cfl * minimum_length / (speed * (2 * polynomial_order + 1)) : 0.0;
}

double euler_cfl_timestep(const TriangleMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl) {
  double minimum_altitude = std::numeric_limits<double>::max();
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = triangle_element_map(mesh, element);
    const auto longest = std::max({map.edge_length(0), map.edge_length(1), map.edge_length(2)});
    minimum_altitude = std::min(minimum_altitude, 4.0 * map.jacobian_determinant() / longest);
  }
  const auto speed =
      euler_maximum_speed(equation, coefficients, mesh.num_elements(), mode_count, 2);
  return speed > 0.0 ? cfl * minimum_altitude / (speed * (2 * polynomial_order + 1)) : 0.0;
}

} // namespace wasl
