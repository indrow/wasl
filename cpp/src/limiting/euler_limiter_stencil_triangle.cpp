#include "wasl/limiting/euler_limiter_stencil.hpp"

#include <vector>

#include "limiting/limiter_stencil_assembly.hpp"
#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/limiting/weno/candidate_transform.hpp"
#include "wasl/mesh/triangle/bounding_box.hpp"
#include "wasl/quadrature/triangle_edge_rules.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {
namespace {

// Each element takes its characteristic normal from the first face that gives
// it a neighbor, oriented outward from that element.
void connect_neighbors(EulerLimiterStencil &stencil, const TriangleMesh &mesh,
                       const std::vector<std::int32_t> &periodic_donors, int edge_point_count) {
  std::vector<unsigned char> has_normal(stencil.element_count, 0);
  const auto set_normal = [&](std::size_t element, double nx, double ny) {
    if (has_normal[element])
      return;
    stencil.normal_x[element] = nx;
    stencil.normal_y[element] = ny;
    has_normal[element] = 1;
  };

  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    const auto left = static_cast<std::size_t>(face.left_element);
    auto right = face.right_element;
    if (face.is_boundary() && !periodic_donors.empty())
      right = periodic_donors[face_index * static_cast<std::size_t>(edge_point_count)];
    if (right == invalid_element)
      continue;
    append_unique_neighbor(stencil, left, right);
    set_normal(left, face.outward_normal_x, face.outward_normal_y);
    if (!face.is_boundary()) {
      append_unique_neighbor(stencil, static_cast<std::size_t>(right), face.left_element);
      set_normal(static_cast<std::size_t>(right), -face.outward_normal_x, -face.outward_normal_y);
    }
  }
}

double centroid_gap(double t0, double t1, double t2, double s0, double s1, double s2) {
  return (t0 + t1 + t2 - s0 - s1 - s2) / 3.0;
}

// Column s of each transform is the projection of the neighbor's mode s.
void build_transforms(EulerLimiterStencil &stencil, const TriangleMesh &mesh) {
  allocate_transforms(stencil);
  const BoundingBox box = bounding_box(mesh);
  std::vector<double> source(stencil.mode_count, 0.0), projected(stencil.mode_count);
  for (std::size_t element = 0; element < stencil.element_count; ++element) {
    const auto target = element_map(mesh, element);
    for (int candidate = 0; candidate < stencil.maximum_neighbors; ++candidate) {
      const auto neighbor = stencil.neighbor(element, candidate);
      if (neighbor == invalid_element)
        continue;
      auto source_map = element_map(mesh, static_cast<std::size_t>(neighbor));
      source_map.shift(periodic_offset(centroid_gap(target.x0, target.x1, target.x2, source_map.x0,
                                                    source_map.x1, source_map.x2),
                                       box.width()),
                       periodic_offset(centroid_gap(target.y0, target.y1, target.y2, source_map.y0,
                                                    source_map.y1, source_map.y2),
                                       box.height()));
      auto *matrix = stencil.transform(element, candidate);
      for (std::size_t source_mode = 0; source_mode < stencil.mode_count; ++source_mode) {
        source[source_mode] = 1.0;
        project_neighbor_polynomial(source_map, target, stencil.polynomial_order, source.data(),
                                    projected.data());
        source[source_mode] = 0.0;
        for (std::size_t target_mode = 0; target_mode < stencil.mode_count; ++target_mode)
          matrix[target_mode * stencil.mode_count + source_mode] = projected[target_mode];
      }
    }
  }
}

void build_check_basis(EulerLimiterStencil &stencil) {
  const auto append = [&](double xi, double eta) {
    const auto offset = stencil.check_basis.size();
    stencil.check_basis.resize(offset + stencil.mode_count);
    dubiner_basis(xi, eta, stencil.polynomial_order, stencil.check_basis.data() + offset);
  };
  const auto volume = triangle_volume_rule(stencil.polynomial_order);
  for (int point = 0; point < volume.size(); ++point)
    append(volume.xi[point], volume.eta[point]);
  for (const auto &edge : triangle_edge_rules(stencil.polynomial_order))
    for (std::size_t point = 0; point < edge.xi.size(); ++point)
      append(edge.xi[point], edge.eta[point]);
}

} // namespace

EulerLimiterStencil make_euler_limiter_stencil(const TriangleMesh &mesh, int polynomial_order,
                                               const std::vector<std::int32_t> &periodic_donors,
                                               int edge_point_count) {
  EulerLimiterStencil stencil;
  stencil.element_count = mesh.num_elements();
  stencil.mode_count = triangle_mode_count(polynomial_order);
  stencil.polynomial_order = polynomial_order;
  stencil.maximum_neighbors = 3;
  stencil.highest_degree_begin =
      polynomial_order > 0 ? triangle_mode_count(polynomial_order - 1) : stencil.mode_count;
  allocate_neighbors(stencil);
  connect_neighbors(stencil, mesh, periodic_donors, edge_point_count);
  build_transforms(stencil, mesh);
  build_check_basis(stencil);
  return stencil;
}

} // namespace wasl
