#include "wasl/case/scalar_diagnostics.hpp"

#include <cmath>

namespace wasl {

// Mode zero is phi_0 = 1/sqrt(2): the average is c_0 / sqrt(2) and the element
// integral is sqrt(2) det(J) c_0.

std::vector<double> scalar_cell_averages(const std::vector<double> &coefficients,
                                         std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(element_count);
  for (std::size_t element = 0; element < element_count; ++element)
    result[element] = coefficients[element * mode_count] / std::sqrt(2.0);
  return result;
}

double scalar_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element)
    result +=
        std::sqrt(2.0) * element_map(mesh, element).jacobian() * coefficients[element * mode_count];
  return result;
}

double scalar_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count) {
  double result = 0.0;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element)
    result += std::sqrt(2.0) * element_map(mesh, element).jacobian_determinant() *
              coefficients[element * mode_count];
  return result;
}

} // namespace wasl
