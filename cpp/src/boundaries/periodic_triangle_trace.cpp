#include "wasl/boundaries/periodic_triangle_trace.hpp"

#include <cmath>
#include <stdexcept>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/triangle/bounding_box.hpp"
#include "wasl/quadrature/triangle_edge_rules.hpp"

namespace wasl {
namespace {

// Translate a point on a box side to the opposite side.
void periodic_image(const TriangleFace &face, const BoundingBox &box, double &x, double &y) {
  if (std::fabs(face.outward_normal_x) > 0.5 && std::fabs(face.outward_normal_y) < 1e-10)
    x += face.outward_normal_x < 0.0 ? box.max_x - box.min_x : box.min_x - box.max_x;
  else if (std::fabs(face.outward_normal_y) > 0.5 && std::fabs(face.outward_normal_x) < 1e-10)
    y += face.outward_normal_y < 0.0 ? box.max_y - box.min_y : box.min_y - box.max_y;
  else
    throw std::runtime_error("simulation: 2D periodic boundaries require an axis-aligned box");
}

// ponytail: linear donor search is quadratic in boundary size; add a spatial
// index when periodic meshes grow beyond a few thousand boundary edges.
std::int32_t find_donor(const TriangleMesh &mesh, double x, double y, double &xi, double &eta) {
  constexpr double tolerance = 1e-10;
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    element_map(mesh, element).physical_to_reference(x, y, xi, eta);
    if (xi >= -1.0 - tolerance && eta >= -1.0 - tolerance && xi + eta <= tolerance)
      return static_cast<std::int32_t>(element);
  }
  throw std::runtime_error("simulation: periodic 2D trace has no donor element");
}

} // namespace

PeriodicTriangleTrace periodic_triangle_trace(const TriangleMesh &mesh, int polynomial_order) {
  const auto rules = triangle_edge_rules(polynomial_order);
  const auto mode_count = triangle_mode_count(polynomial_order);
  const int point_count = polynomial_order + 1;
  const BoundingBox box = bounding_box(mesh);

  PeriodicTriangleTrace result;
  result.elements.assign(mesh.faces.size() * point_count, invalid_element);
  result.basis.assign(mesh.faces.size() * point_count * mode_count, 0.0);
  result.point_count = point_count;
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    if (!face.is_boundary())
      continue;
    const auto map = element_map(mesh, static_cast<std::size_t>(face.left_element));
    const auto &rule = rules.at(static_cast<std::size_t>(face.left_local_edge));
    for (int point = 0; point < point_count; ++point) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(rule.xi[point], rule.eta[point], x, y);
      periodic_image(face, box, x, y);
      const auto index = face_index * point_count + static_cast<std::size_t>(point);
      double xi = 0.0, eta = 0.0;
      result.elements[index] = find_donor(mesh, x, y, xi, eta);
      dubiner_basis(xi, eta, polynomial_order, result.basis.data() + index * mode_count);
    }
  }
  return result;
}

std::vector<double> periodic_triangle_exterior(const PeriodicTriangleTrace &trace,
                                               const std::vector<double> &coefficients,
                                               std::size_t mode_count) {
  std::vector<double> exterior(trace.elements.size(), 0.0);
  for (std::size_t point = 0; point < trace.elements.size(); ++point) {
    if (trace.elements[point] == invalid_element)
      continue;
    const auto donor = static_cast<std::size_t>(trace.elements[point]);
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      exterior[point] +=
          coefficients[donor * mode_count + mode] * trace.basis[point * mode_count + mode];
  }
  return exterior;
}

std::vector<double> periodic_triangle_euler_exterior(const PeriodicTriangleTrace &trace,
                                                     const std::vector<double> &coefficients,
                                                     std::size_t element_count,
                                                     std::size_t mode_count,
                                                     std::size_t face_count) {
  const ModalLayout layout{mode_count, element_count};
  const auto point_count = static_cast<std::size_t>(trace.point_count);
  std::vector<double> exterior(Euler::variable_count * trace.elements.size(), 0.0);
  for (std::size_t point = 0; point < trace.elements.size(); ++point) {
    if (trace.elements[point] == invalid_element)
      continue;
    const auto donor = static_cast<std::size_t>(trace.elements[point]);
    const auto face = point / point_count;
    const auto quadrature_point = point % point_count;
    for (std::size_t variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        exterior[(variable * face_count + face) * point_count + quadrature_point] +=
            coefficients[layout.index(variable, mode, donor)] *
            trace.basis[point * mode_count + mode];
  }
  return exterior;
}

} // namespace wasl
