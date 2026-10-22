#include "wasl/mesh/triangle/triangle_mesh.hpp"

#include <algorithm>

namespace wasl {

TriangleMap element_map(const TriangleMesh &mesh, std::size_t element) {
  return mesh.elements.at(element).map();
}

std::size_t count_interior_faces(const TriangleMesh &mesh) noexcept {
  return static_cast<std::size_t>(
      std::count_if(mesh.faces.begin(), mesh.faces.end(),
                    [](const TriangleFace &face) { return !face.is_boundary(); }));
}

std::size_t count_boundary_faces(const TriangleMesh &mesh) noexcept {
  return mesh.faces.size() - count_interior_faces(mesh);
}

} // namespace wasl
