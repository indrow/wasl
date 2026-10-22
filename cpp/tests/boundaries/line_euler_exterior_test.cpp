// Euler exterior states at the ends of a line: periodic ends read the opposite
// element's far trace, outflow ends their owner's near trace, and inflow is
// rejected.

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/line_euler_exterior.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line/line_mesh.hpp"

namespace {

constexpr double kTolerance = 1e-14;
constexpr int kOrder = 2;

wasl::LineMesh three_elements() {
  wasl::LineMesh mesh;
  mesh.elements = {{1, 2, 0.0, 0.3, 1}, {2, 3, 0.3, 0.7, 1}, {3, 4, 0.7, 1.0, 1}};
  mesh.faces = {{0, wasl::invalid_element, 0.0, -1.0, 2},
                {0, 1, 0.3, 1.0, wasl::invalid_physical},
                {1, 2, 0.7, 1.0, wasl::invalid_physical},
                {2, wasl::invalid_element, 1.0, 1.0, 3}};
  return mesh;
}

// Trace of variable v of one element at xi.
double trace(const std::vector<double> &coefficients, const wasl::ModalLayout &layout,
             std::size_t variable, std::size_t element, double xi) {
  std::vector<double> basis(layout.mode_count);
  wasl::legendre_basis(xi, kOrder, basis.data());
  double value = 0.0;
  for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
    value += coefficients[layout.index(variable, mode, element)] * basis[mode];
  return value;
}

bool ends_match(const wasl::LineBoundaryConditions &conditions, std::size_t left_donor,
                double left_xi, std::size_t right_donor, double right_xi) {
  const auto mesh = three_elements();
  const wasl::ModalLayout layout{wasl::line_mode_count(kOrder), mesh.num_elements()};
  std::vector<double> coefficients(wasl::Euler::variable_count * layout.mode_count *
                                   layout.element_count);
  for (std::size_t index = 0; index < coefficients.size(); ++index)
    coefficients[index] = 0.1 * static_cast<double>(index % 7) - 0.05 * static_cast<double>(index);
  const auto exterior = wasl::line_euler_exterior(mesh, kOrder, conditions, coefficients);
  const auto faces = mesh.faces.size();
  for (std::size_t variable = 0; variable < wasl::Euler::variable_count; ++variable) {
    if (std::fabs(exterior[variable * faces] -
                  trace(coefficients, layout, variable, left_donor, left_xi)) > kTolerance ||
        std::fabs(exterior[variable * faces + faces - 1] -
                  trace(coefficients, layout, variable, right_donor, right_xi)) > kTolerance ||
        exterior[variable * faces + 1] != 0.0 || exterior[variable * faces + 2] != 0.0)
      return false;
  }
  return true;
}

bool inflow_is_rejected() {
  const wasl::BoundaryCondition inflow{wasl::BoundaryKind::inflow, 1.0};
  try {
    wasl::line_euler_exterior(three_elements(), kOrder, {inflow, inflow},
                              std::vector<double>(36, 0.0));
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

} // namespace

int main() {
  const wasl::BoundaryCondition outflow{wasl::BoundaryKind::outflow, 0.0};
  const bool ok = ends_match({}, 2, 1.0, 0, -1.0) &&
                  ends_match({outflow, outflow}, 0, -1.0, 2, 1.0) && inflow_is_rejected();
  return ok ? 0 : 1;
}
