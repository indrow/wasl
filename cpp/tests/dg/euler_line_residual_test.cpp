#include <algorithm>
#include <cmath>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/dg/euler/line_residual.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/mesh/line/line_mesh.hpp"

namespace {

double component(const wasl::EulerConserved &state, int variable) {
  if (variable == 0)
    return state.density;
  if (variable == 1)
    return state.momentum_x;
  if (variable == 2)
    return state.momentum_y;
  return state.total_energy;
}

} // namespace

int main() {
  wasl::LineMesh mesh;
  mesh.elements = {{1, 2, 0.0, 0.4, 1}, {2, 3, 0.4, 1.0, 1}};
  mesh.faces = {{0, wasl::invalid_element, 0.0, -1.0, 2},
                {0, 1, 0.4, 1.0, wasl::invalid_physical},
                {1, wasl::invalid_element, 1.0, 1.0, 3}};
  constexpr int polynomial_order = 2;
  const auto mode_count = wasl::line_mode_count(polynomial_order);
  const auto element_count = mesh.num_elements();
  const wasl::Euler equation{1.4};
  const auto state = wasl::conserved(equation, {1.1, 0.4, 0.0, 0.9});
  std::vector<double> coefficients(wasl::Euler::variable_count * mode_count * element_count, 0.0);
  for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
    for (std::size_t element = 0; element < element_count; ++element)
      coefficients[static_cast<std::size_t>(variable) * mode_count * element_count + element] =
          component(state, variable) * std::sqrt(2.0);
  std::vector<double> exterior(wasl::Euler::variable_count * mesh.faces.size());
  for (int variable = 0; variable < wasl::Euler::variable_count; ++variable)
    for (std::size_t face = 0; face < mesh.faces.size(); ++face)
      exterior[static_cast<std::size_t>(variable) * mesh.faces.size() + face] =
          component(state, variable);
  std::vector<double> rhs(coefficients.size());
  wasl::euler_line_residual(mesh, equation, wasl::hllc_flux, polynomial_order, coefficients.data(),
                            exterior.data(), rhs.data());
  return std::any_of(rhs.begin(), rhs.end(), [](double value) { return std::fabs(value) > 2e-12; });
}
