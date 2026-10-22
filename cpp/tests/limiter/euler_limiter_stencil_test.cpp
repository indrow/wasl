// Line limiter stencils: periodic ends join the boundary elements, open ends
// leave each boundary element with its single interior neighbor.

#include <cstdint>
#include <vector>

#include "wasl/limiting/euler_limiter_stencil.hpp"
#include "wasl/mesh/line/line_mesh.hpp"

namespace {

wasl::LineMesh four_elements() {
  wasl::LineMesh mesh;
  for (int element = 0; element < 4; ++element)
    mesh.elements.push_back({element + 1, element + 2, 0.25 * element, 0.25 * (element + 1), 1});
  mesh.faces.push_back({0, wasl::invalid_element, 0.0, -1.0, wasl::invalid_physical});
  for (int element = 0; element < 3; ++element)
    mesh.faces.push_back({element, element + 1, 0.25 * (element + 1), 1.0, wasl::invalid_physical});
  mesh.faces.push_back({3, wasl::invalid_element, 1.0, 1.0, wasl::invalid_physical});
  return mesh;
}

std::vector<std::int32_t> neighbors(const wasl::EulerLimiterStencil &stencil, int element) {
  std::vector<std::int32_t> result;
  for (int candidate = 0; candidate < stencil.maximum_neighbors; ++candidate)
    if (stencil.neighbor(element, candidate) != wasl::invalid_element)
      result.push_back(stencil.neighbor(element, candidate));
  return result;
}

} // namespace

int main() {
  const auto mesh = four_elements();
  const auto periodic = wasl::make_euler_limiter_stencil(mesh, 2, wasl::LineEnds::periodic);
  const auto open = wasl::make_euler_limiter_stencil(mesh, 2, wasl::LineEnds::open);
  using Neighbors = std::vector<std::int32_t>;
  const bool ok = neighbors(periodic, 0) == Neighbors{1, 3} &&
                  neighbors(periodic, 3) == Neighbors{2, 0} && neighbors(open, 0) == Neighbors{1} &&
                  neighbors(open, 3) == Neighbors{2} && neighbors(open, 1) == Neighbors{0, 2} &&
                  neighbors(periodic, 2) == Neighbors{1, 3};
  return ok ? 0 : 1;
}
