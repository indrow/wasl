#include "wasl/case/euler_diagnostics.hpp"

#include <cmath>

#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/euler.hpp"

namespace wasl {

std::vector<double> euler_averages(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count) {
  const ModalLayout layout{mode_count, element_count};
  std::vector<double> result(Euler::variable_count * element_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto state = euler_cell_average(layout, coefficients.data(), element);
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      result[element * Euler::variable_count + static_cast<std::size_t>(variable)] =
          euler_component(state, variable);
  }
  return result;
}

double euler_density_mass(const LineMesh &mesh, const std::vector<double> &coefficients) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element)
    result += std::sqrt(2.0) * element_map(mesh, element).jacobian() * coefficients[element];
  return result;
}

double euler_density_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element)
    result +=
        std::sqrt(2.0) * element_map(mesh, element).jacobian_determinant() * coefficients[element];
  return result;
}

} // namespace wasl
