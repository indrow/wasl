// Tests for MSH 4.1 triangle reading and canonical triangle mesh topology.
//
// The fixture is a hand-authored MSH 4.1 file for the unit square split into two
// triangles by the diagonal (0,0)-(1,1). It matches the exact format Gmsh emits.
// Element order is scrambled (T1 listed before T0) and one triangle is written
// with reversed vertex order (clockwise), so the tests exercise orientation
// normalization and deterministic edge ownership independent of record order.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include "wasl/mesh/line_mesh.hpp" // invalid_element, invalid_physical
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/triangle_geometry.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace {

bool approx(double a, double b, double tol = 1e-12) { return std::fabs(a - b) <= tol; }

// Unit square [0,1]x[0,1] split by the diagonal node1-node3 into two triangles:
//   T0: nodes 1,2,3 (lower right, CCW)
//   T1: nodes 1,3,4 (upper left, CCW)
// The file lists T1 first and writes T1 reversed as 1,4,3 (clockwise) to test
// orientation normalization. Boundary physical tags (curves):
//   bottom (nodes 1-2): physical 10 "bottom"
//   right  (nodes 2-3): physical 11 "right"
//   top    (nodes 4-3): physical 12 "top"
//   left   (nodes 1-4): physical 13 "left"
// Surface (both triangles): physical 1 "fluid".
constexpr const char *kSquareMsh = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
1 10 "bottom"
1 11 "right"
1 12 "top"
1 13 "left"
2 1 "fluid"
$EndPhysicalNames
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 10 2 1 -2
2 1 0 0 1 1 0 1 11 2 2 -3
3 0 1 0 1 1 0 1 12 2 4 -3
4 0 0 0 0 1 0 1 13 2 1 -4
1 0 0 0 1 1 0 1 1 4 1 2 3 4
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
2 1 4 3
1 1 2 3
$EndElements
)MSH";

std::string write_fixture(const char *content) {
  char tmpl[] = "/tmp/wasl_fixture_XXXXXX";
  const int fd = mkstemp(tmpl);
  if (fd < 0)
    return {};
  close(fd);
  const std::string path = tmpl;
  std::ofstream out(path, std::ios::trunc);
  if (!out) {
    std::remove(path.c_str());
    return {};
  }
  out << content;
  out.close();
  return path;
}

bool check_parse_counts() {
  const std::string path = write_fixture(kSquareMsh);
  const auto doc = wasl::read_msh(path);
  std::remove(path.c_str());

  if (doc.nodes.size() != 4)
    return false;
  if (doc.elements.size() != 2)
    return false;
  if (doc.surface_entities.size() != 1)
    return false;
  if (doc.curve_entities.size() != 4)
    return false;
  if (doc.physical_names.size() != 5)
    return false;
  return true;
}

bool check_orientation_and_counts() {
  const std::string path = write_fixture(kSquareMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_triangle_mesh(doc);
  std::remove(path.c_str());

  if (mesh.num_elements() != 2)
    return false;
  // Every element must have positive (CCW) orientation regardless of input.
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (wasl::triangle_element_map(mesh, e).orientation() != 1)
      return false;
  }
  // Each triangle has area 0.5; the diagonal edge is shared.
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!approx(wasl::triangle_element_map(mesh, e).area(), 0.5))
      return false;
  }
  return true;
}

// Two triangles share one interior edge (the diagonal). Four boundary edges and
// one interior edge give five faces total.
bool check_edge_topology() {
  const std::string path = write_fixture(kSquareMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_triangle_mesh(doc);
  std::remove(path.c_str());

  if (mesh.faces.size() != 5)
    return false;
  if (wasl::count_interior_faces(mesh) != 1)
    return false;
  if (wasl::count_boundary_faces(mesh) != 4)
    return false;

  // Find the interior face and verify deterministic left<right ownership and
  // that both sides record their traversal.
  const wasl::TriangleFace *interior = nullptr;
  for (const auto &face : mesh.faces) {
    if (face.right_element != wasl::invalid_element) {
      interior = &face;
      break;
    }
  }
  if (interior == nullptr)
    return false;
  if (interior->left_element >= interior->right_element)
    return false;
  // The diagonal is between nodes 1 and 3; the sorted key is (1, 3).
  if (interior->vertex_a != 1 || interior->vertex_b != 3)
    return false;
  // Per-side vertex order must be recorded on both sides.
  if (interior->left_vertex0 == 0 || interior->right_vertex0 == 0)
    return false;

  // Boundary faces have exactly one owner and an invalid right element.
  for (const auto &face : mesh.faces) {
    if (face.right_element == wasl::invalid_element) {
      if (face.left_element == wasl::invalid_element)
        return false;
    }
  }
  return true;
}

// Boundary outward normals point out of the square. The bottom edge normal is
// (0, -1) and the top edge normal is (0, 1) (up to numerical round-off).
bool check_boundary_normals() {
  const std::string path = write_fixture(kSquareMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_triangle_mesh(doc);
  std::remove(path.c_str());

  for (const auto &face : mesh.faces) {
    if (face.right_element != wasl::invalid_element)
      continue;
    // Outward normals are unit length.
    if (!approx(face.outward_normal_x * face.outward_normal_x +
                    face.outward_normal_y * face.outward_normal_y,
                1.0))
      return false;
    // Identify the bottom (nodes 1,2) and top (nodes 3,4) edges.
    const bool is_bottom = (face.vertex_a == 1 && face.vertex_b == 2);
    const bool is_top = (face.vertex_a == 3 && face.vertex_b == 4);
    if (is_bottom && !(approx(face.outward_normal_x, 0.0) && approx(face.outward_normal_y, -1.0)))
      return false;
    if (is_top && !(approx(face.outward_normal_x, 0.0) && approx(face.outward_normal_y, 1.0)))
      return false;
  }
  return true;
}

bool check_physical_groups() {
  const std::string path = write_fixture(kSquareMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_triangle_mesh(doc);
  std::remove(path.c_str());

  // Both elements belong to the "fluid" surface physical group (tag 1).
  if (mesh.physical_to_elements.at(1).size() != 2)
    return false;
  if (mesh.physical_names.at({2, 1}) != "fluid")
    return false;
  if (mesh.physical_names.at({1, 10}) != "bottom")
    return false;
  // Each of the four boundary curves owns one face.
  for (int tag : {10, 11, 12, 13}) {
    if (mesh.physical_to_faces.at(tag).size() != 1)
      return false;
  }
  return true;
}

bool check_unsupported_element_type() {
  // A 4-node quad (element type 3) must be rejected explicitly.
  constexpr const char *quad = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 10 2 1 -2
2 1 0 0 1 1 0 1 11 2 2 -3
3 0 1 0 1 1 0 1 12 2 4 -3
4 0 0 0 0 1 0 1 13 2 1 -4
1 0 0 0 1 1 0 1 1 4 1 2 3 4
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 1 1 1
2 1 3 1
1 1 2 3 4
$EndElements
)MSH";
  const std::string path = write_fixture(quad);
  bool threw = false;
  try {
    static_cast<void>(wasl::read_msh(path));
  } catch (const std::runtime_error &error) {
    threw = std::string{error.what()}.find("unsupported element type") != std::string::npos;
  }
  std::remove(path.c_str());
  return threw;
}

bool check_missing_node_error() {
  wasl::MshDocument document;
  document.nodes.push_back(wasl::MshNode{1, 2, 1, 0.0, 0.0, 0.0});
  document.nodes.push_back(wasl::MshNode{2, 2, 1, 1.0, 0.0, 0.0});
  // References node 3, which is absent.
  document.elements.push_back(wasl::MshElement{1, 2, 1, 2, {1, 2, 3}});
  try {
    static_cast<void>(wasl::build_triangle_mesh(document));
  } catch (const std::runtime_error &error) {
    return std::string{error.what()}.find("missing node 3") != std::string::npos;
  }
  return false;
}

bool check_no_triangles_error() {
  wasl::MshDocument document;
  // Only a line element: no triangles.
  document.elements.push_back(wasl::MshElement{1, 1, 1, 1, {1, 2}});
  try {
    static_cast<void>(wasl::build_triangle_mesh(document));
  } catch (const std::runtime_error &error) {
    return std::string{error.what()}.find("no triangle elements") != std::string::npos;
  }
  return false;
}

} // namespace

int main() {
  if (!check_parse_counts())
    return 1;
  if (!check_orientation_and_counts())
    return 1;
  if (!check_edge_topology())
    return 1;
  if (!check_boundary_normals())
    return 1;
  if (!check_physical_groups())
    return 1;
  if (!check_unsupported_element_type())
    return 1;
  if (!check_missing_node_error())
    return 1;
  if (!check_no_triangles_error())
    return 1;
  return 0;
}
