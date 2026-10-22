#include "wasl/case/initial_state.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/quadrature.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace wasl {
namespace {

double initial_value(const ScalarAdvectionCase &config, double x) {
  if (config.initial_sine) {
    return config.initial_sine_amplitude * std::sin(2.0 * std::numbers::pi * x);
  }
  for (const auto &region : config.initial_intervals) {
    if (x >= region.lower && x <= region.upper) {
      return region.value;
    }
  }
  return config.initial_default;
}

double initial_value(const ScalarAdvectionCase &config, double x, double y) {
  for (const auto &region : config.initial_circles) {
    const auto dx = x - region.center_x;
    const auto dy = y - region.center_y;
    if (dx * dx + dy * dy <= region.radius * region.radius) {
      return region.value;
    }
  }
  return config.initial_default;
}

double euler_component(const EulerConserved &state, int variable) {
  if (variable == 0)
    return state.density;
  if (variable == 1)
    return state.momentum_x;
  if (variable == 2)
    return state.momentum_y;
  return state.total_energy;
}

EulerPrimitive initial_value(const EulerCase &config, double x, double y = 0.0) {
  const EulerRegion *fallback = nullptr;
  for (const auto &region : config.initial_regions) {
    if (region.shape == EulerRegion::Shape::default_region) {
      fallback = &region;
    } else if (region.shape == EulerRegion::Shape::interval && x >= region.lower &&
               x <= region.upper) {
      return region.value;
    } else if (region.shape == EulerRegion::Shape::circle) {
      const auto dx = x - region.center_x;
      const auto dy = y - region.center_y;
      if (dx * dx + dy * dy <= region.radius * region.radius)
        return region.value;
    }
  }
  if (fallback == nullptr)
    throw std::runtime_error("simulation: Euler initial state requires a default region");
  return fallback->value;
}

} // namespace

std::vector<double> initial_coefficients(const LineMesh &mesh, const EulerCase &config,
                                         const Euler &equation) {
  const auto mode_count = line_mode_count(config.polynomial_order);
  const auto element_count = mesh.num_elements();
  const auto quadrature = gauss_legendre(7);
  std::vector<double> coefficients(Euler::variable_count * mode_count * element_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto map = element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      const auto xi = quadrature.points[point];
      legendre_basis(xi, config.polynomial_order, basis.data());
      const auto state = conserved(equation, initial_value(config, map.reference_to_physical(xi)));
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                       element] +=
              quadrature.weights[point] * euler_component(state, variable) * basis[mode];
    }
  }
  return coefficients;
}

std::vector<double> initial_coefficients(const TriangleMesh &mesh, const EulerCase &config,
                                         const Euler &equation) {
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto element_count = mesh.num_elements();
  const auto quadrature = triangle_volume_rule(std::max(config.polynomial_order, 5));
  std::vector<double> coefficients(Euler::variable_count * mode_count * element_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto map = triangle_element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(quadrature.xi[point], quadrature.eta[point], x, y);
      dubiner_basis(quadrature.xi[point], quadrature.eta[point], config.polynomial_order,
                    basis.data());
      const auto state = conserved(equation, initial_value(config, x, y));
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                       element] +=
              quadrature.weights[point] * euler_component(state, variable) * basis[mode];
    }
  }
  return coefficients;
}

std::vector<double> initial_coefficients(const LineMesh &mesh, const ScalarAdvectionCase &config) {
  const auto mode_count = line_mode_count(config.polynomial_order);
  const auto quadrature = gauss_legendre(7);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      const auto xi = quadrature.points[point];
      legendre_basis(xi, config.polynomial_order, basis.data());
      const auto value = initial_value(config, map.reference_to_physical(xi));
      for (std::size_t mode = 0; mode < mode_count; ++mode) {
        coefficients[element * mode_count + mode] +=
            quadrature.weights[point] * value * basis[mode];
      }
    }
  }
  return coefficients;
}

std::vector<double> initial_coefficients(const TriangleMesh &mesh,
                                         const ScalarAdvectionCase &config) {
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto quadrature = triangle_volume_rule(std::max(config.polynomial_order, 5));
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  std::vector<double> basis(mode_count);
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = triangle_element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      double x = 0.0;
      double y = 0.0;
      map.reference_to_physical(quadrature.xi[point], quadrature.eta[point], x, y);
      dubiner_basis(quadrature.xi[point], quadrature.eta[point], config.polynomial_order,
                    basis.data());
      const auto value = initial_value(config, x, y);
      for (std::size_t mode = 0; mode < mode_count; ++mode) {
        coefficients[element * mode_count + mode] +=
            quadrature.weights[point] * value * basis[mode];
      }
    }
  }
  return coefficients;
}

} // namespace wasl
