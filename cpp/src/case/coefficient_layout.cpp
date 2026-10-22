#include "wasl/case/coefficient_layout.hpp"

#include <cmath>

#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {

std::vector<double> checkpoint_coefficients(const std::vector<double> &coefficients,
                                            std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(coefficients.size());
  for (std::size_t element = 0; element < element_count; ++element) {
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      result[mode * element_count + element] = coefficients[element * mode_count + mode];
    }
  }
  return result;
}

std::vector<double> solver_coefficients(const std::vector<double> &coefficients,
                                        std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(coefficients.size());
  for (std::size_t element = 0; element < element_count; ++element) {
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      result[element * mode_count + mode] = coefficients[mode * element_count + element];
    }
  }
  return result;
}

std::vector<double> solver_coefficients(const Checkpoint &checkpoint) {
  return solver_coefficients(checkpoint.coefficients, checkpoint.element_count,
                             checkpoint.mode_count);
}

std::vector<double> scalar_cell_averages(const std::vector<double> &coefficients,
                                         std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(element_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    result[element] = coefficients[element * mode_count] / std::sqrt(2.0);
  }
  return result;
}

double scalar_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    result +=
        std::sqrt(2.0) * element_map(mesh, element).jacobian() * coefficients[element * mode_count];
  }
  return result;
}

double scalar_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    result += std::sqrt(2.0) * triangle_element_map(mesh, element).jacobian_determinant() *
              coefficients[element * mode_count];
  }
  return result;
}

EulerConserved average_euler_state(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count,
                                   std::size_t element) {
  const auto scale = 1.0 / std::sqrt(2.0);
  return {scale * coefficients[element], scale * coefficients[mode_count * element_count + element],
          scale * coefficients[2 * mode_count * element_count + element],
          scale * coefficients[3 * mode_count * element_count + element]};
}

std::vector<double> euler_averages(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(Euler::variable_count * element_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto state = average_euler_state(coefficients, element_count, mode_count, element);
    result[element * Euler::variable_count + 0] = state.density;
    result[element * Euler::variable_count + 1] = state.momentum_x;
    result[element * Euler::variable_count + 2] = state.momentum_y;
    result[element * Euler::variable_count + 3] = state.total_energy;
  }
  return result;
}

double euler_density_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                          std::size_t) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    result += std::sqrt(2.0) * element_map(mesh, element).jacobian() * coefficients[element];
  }
  return result;
}

double euler_density_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                          std::size_t) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    result += std::sqrt(2.0) * triangle_element_map(mesh, element).jacobian_determinant() *
              coefficients[element];
  }
  return result;
}

} // namespace wasl
