// Setup data for immediate-neighbor Euler limiting, computed once per mesh.
//
// For each element the stencil stores its face neighbors (periodic images
// included), the matrices mapping a neighbor's modes into the element's basis,
// the normal used for the characteristic decomposition, and the basis values
// at the positivity check points.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

struct EulerLimiterStencil {
  std::size_t element_count = 0;
  std::size_t mode_count = 0;
  int polynomial_order = 0;
  int maximum_neighbors = 0;
  std::size_t highest_degree_begin = 0;   // first mode of the highest degree
  std::vector<std::int32_t> neighbors;    // [element][candidate], invalid_element if absent
  std::vector<double> transforms;         // [element][candidate][target mode][source mode]
  std::vector<double> normal_x, normal_y; // [element] characteristic direction
  std::vector<double> check_basis;        // [check point][mode]

  std::int32_t neighbor(std::size_t element, int candidate) const {
    return neighbors[element * static_cast<std::size_t>(maximum_neighbors) +
                     static_cast<std::size_t>(candidate)];
  }

  const double *transform(std::size_t element, int candidate) const {
    return transforms.data() + transform_offset(element, candidate);
  }
  double *transform(std::size_t element, int candidate) {
    return transforms.data() + transform_offset(element, candidate);
  }

  std::size_t check_point_count() const { return check_basis.size() / mode_count; }

private:
  std::size_t transform_offset(std::size_t element, int candidate) const {
    return (element * static_cast<std::size_t>(maximum_neighbors) +
            static_cast<std::size_t>(candidate)) *
           mode_count * mode_count;
  }
};

// Whether the two boundary elements of a line are neighbors across the ends.
enum class LineEnds { periodic, open };

// Line elements use their endpoint neighbors. Periodic ends join the two
// boundary elements; open ends leave them with one candidate each, so a state
// is never blended with the far side of the domain. Check points are the Gauss
// points and endpoints.
EulerLimiterStencil make_euler_limiter_stencil(const LineMesh &mesh, int polynomial_order,
                                               LineEnds ends);

// Triangle elements use up to three edge neighbors. periodic_donors, packed
// [face][edge point] as in PeriodicTriangleTrace, supplies boundary neighbors
// when non-empty. Check points are the volume and edge quadrature points.
EulerLimiterStencil make_euler_limiter_stencil(const TriangleMesh &mesh, int polynomial_order,
                                               const std::vector<std::int32_t> &periodic_donors,
                                               int edge_point_count);

} // namespace wasl
