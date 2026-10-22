#include "wasl/time/cfl_timestep.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/primitive.hpp"

namespace wasl {
namespace {

double minimum_length(const LineMesh &mesh) {
  double result = std::numeric_limits<double>::max();
  for (const auto &element : mesh.elements)
    result = std::min(result, element.x_right - element.x_left);
  return result;
}

// Smallest altitude 2 * area / longest edge, with area = 2 det(J).
double minimum_altitude(const TriangleMesh &mesh) {
  double result = std::numeric_limits<double>::max();
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const auto map = element_map(mesh, element);
    const auto longest = std::max({map.edge_length(0), map.edge_length(1), map.edge_length(2)});
    result = std::min(result, 4.0 * map.jacobian_determinant() / longest);
  }
  return result;
}

// Largest |u| + c over cell averages; 1D uses |u_x|, 2D the speed |u|.
double euler_maximum_speed(const Euler &equation, const std::vector<double> &coefficients,
                           std::size_t element_count, std::size_t mode_count, int dimensions) {
  const ModalLayout layout{mode_count, element_count};
  double result = 0.0;
  for (std::size_t element = 0; element < element_count; ++element) {
    const auto state = euler_cell_average(layout, coefficients.data(), element);
    if (!is_physical(equation, state))
      throw std::runtime_error("simulation: Euler cell average is not physical");
    const auto value = primitive(equation, state);
    const auto velocity = dimensions == 1 ? std::fabs(value.velocity_x)
                                          : std::hypot(value.velocity_x, value.velocity_y);
    result = std::max(result, velocity + sound_speed(equation, value));
  }
  return result;
}

double euler_limit(double length, double speed, int polynomial_order, double cfl) {
  return speed > 0.0 ? cfl * length / (speed * (2 * polynomial_order + 1)) : 0.0;
}

} // namespace

double advection_cfl_timestep(const LineMesh &mesh, const ScalarAdvection &eq, int polynomial_order,
                              double cfl) {
  if (mesh.elements.empty() || eq.wavespeed == 0.0)
    return 0.0;
  return cfl * minimum_length(mesh) / (std::fabs(eq.wavespeed) * (2.0 * polynomial_order + 1.0));
}

double triangle_scalar_cfl_timestep(const TriangleMesh &mesh, double velocity_x, double velocity_y,
                                    int polynomial_order, double cfl) {
  const auto speed = std::hypot(velocity_x, velocity_y);
  if (mesh.elements.empty() || speed == 0.0)
    return 0.0;
  return cfl * minimum_altitude(mesh) / (speed * (2.0 * polynomial_order + 1.0));
}

double euler_cfl_timestep(const LineMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl) {
  const auto speed =
      euler_maximum_speed(equation, coefficients, mesh.num_elements(), mode_count, 1);
  return euler_limit(minimum_length(mesh), speed, polynomial_order, cfl);
}

double euler_cfl_timestep(const TriangleMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl) {
  const auto speed =
      euler_maximum_speed(equation, coefficients, mesh.num_elements(), mode_count, 2);
  return euler_limit(minimum_altitude(mesh), speed, polynomial_order, cfl);
}

} // namespace wasl
