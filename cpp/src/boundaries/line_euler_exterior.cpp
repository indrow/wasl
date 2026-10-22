#include "wasl/boundaries/line_euler_exterior.hpp"

#include <cstddef>
#include <stdexcept>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/euler.hpp"

namespace wasl {

std::vector<double> line_euler_exterior(const LineMesh &mesh, int polynomial_order,
                                        const LineBoundaryConditions &conditions,
                                        const std::vector<double> &coefficients) {
  const ModalLayout layout{line_mode_count(polynomial_order), mesh.num_elements()};
  std::vector<double> left_basis(layout.mode_count), right_basis(layout.mode_count);
  legendre_basis(-1.0, polynomial_order, left_basis.data());
  legendre_basis(1.0, polynomial_order, right_basis.data());

  std::vector<double> exterior(Euler::variable_count * mesh.faces.size(), 0.0);
  for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
    const auto &value = mesh.faces[face];
    if (!value.is_boundary())
      continue;
    const bool left_end = value.outward_normal < 0.0;
    const auto kind = (left_end ? conditions.left : conditions.right).kind;
    if (kind == BoundaryKind::inflow)
      throw std::invalid_argument("line Euler exterior: inflow is not implemented");
    // Periodic ends read the opposite element's far endpoint; outflow ends
    // read their own owner's near endpoint.
    const bool periodic = kind == BoundaryKind::periodic;
    const auto donor = periodic ? (left_end ? layout.element_count - 1 : 0)
                                : static_cast<std::size_t>(value.left_element);
    const bool right_endpoint = periodic ? left_end : !left_end;
    const auto &basis = right_endpoint ? right_basis : left_basis;
    for (std::size_t variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
        exterior[variable * mesh.faces.size() + face] +=
            coefficients[layout.index(variable, mode, donor)] * basis[mode];
  }
  return exterior;
}

} // namespace wasl
