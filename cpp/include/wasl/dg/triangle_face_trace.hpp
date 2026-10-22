// Map a quadrature point on a triangle face from the left element's reference
// coordinates into the right element's, and scale reference edge weights to
// physical length.

#pragma once

#include "wasl/mesh/triangle/triangle_geometry.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// The left point is first reduced to the edge parameter in [-1, 1] that runs
// along the left element's traversal direction. Opposite traversals (the usual
// interior case) reverse the parameter; the right element's local edge then
// places it on its own reference edge.
inline void neighbor_edge_point(const TriangleFace &face, double left_xi, double left_eta,
                                double &right_xi, double &right_eta) noexcept {
  double parameter = face.left_local_edge == 0   ? left_xi
                     : face.left_local_edge == 1 ? left_eta
                                                 : -left_eta;
  if (face.left_vertex0 != face.right_vertex0)
    parameter = -parameter;
  right_xi = -1.0;
  right_eta = -1.0;
  if (face.right_local_edge == 0) {
    right_xi = parameter;
  } else if (face.right_local_edge == 1) {
    right_xi = -parameter;
    right_eta = parameter;
  } else {
    right_eta = -parameter;
  }
}

// Physical over reference length of the left element's face edge.
inline double face_edge_scale(const TriangleMap &left_map, int left_local_edge) noexcept {
  return left_map.edge_length(left_local_edge) / reference_edge_length(left_local_edge);
}

} // namespace wasl
