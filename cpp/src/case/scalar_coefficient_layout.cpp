#include "wasl/case/scalar_coefficient_layout.hpp"

namespace wasl {

std::vector<double> checkpoint_coefficients(const std::vector<double> &coefficients,
                                            std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(coefficients.size());
  for (std::size_t element = 0; element < element_count; ++element)
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      result[mode * element_count + element] = coefficients[element * mode_count + mode];
  return result;
}

std::vector<double> solver_coefficients(const std::vector<double> &coefficients,
                                        std::size_t element_count, std::size_t mode_count) {
  std::vector<double> result(coefficients.size());
  for (std::size_t element = 0; element < element_count; ++element)
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      result[element * mode_count + mode] = coefficients[mode * element_count + element];
  return result;
}

std::vector<double> solver_coefficients(const Checkpoint &checkpoint) {
  return solver_coefficients(checkpoint.coefficients, checkpoint.element_count,
                             checkpoint.mode_count);
}

} // namespace wasl
