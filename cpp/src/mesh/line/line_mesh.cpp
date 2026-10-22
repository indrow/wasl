#include "wasl/mesh/line/line_mesh.hpp"

#include <algorithm>

namespace wasl {

LineMap element_map(const LineMesh &mesh, std::size_t element) {
  return mesh.elements.at(element).map();
}

std::size_t count_interior_faces(const LineMesh &mesh) noexcept {
  return static_cast<std::size_t>(
      std::count_if(mesh.faces.begin(), mesh.faces.end(),
                    [](const LineFace &face) { return !face.is_boundary(); }));
}

std::size_t count_boundary_faces(const LineMesh &mesh) noexcept {
  return mesh.faces.size() - count_interior_faces(mesh);
}

} // namespace wasl
