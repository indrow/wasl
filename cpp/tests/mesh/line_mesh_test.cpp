// Tests for MSH 4.1 reading and canonical line mesh construction.
//
// The fixture is a hand-authored MSH 4.1 file with scrambled element order and
// one reversed element, so the tests exercise orientation normalization and
// deterministic left/right connectivity independent of record order.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "wasl/mesh/line/build_line_mesh.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/msh/msh_reader.hpp"

namespace {

// A four-element line mesh on [0, 1] with deliberately scrambled element order
// in the MSH file and one reversed element. Gmsh assigns element tags in write
// order; Wasl re-sorts by left x, so element 0 is always the leftmost cell.
// Curve 1 is reversed: its MSH record lists node 4 before node 3.
constexpr const char *kScrambledMsh = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
0 2 "left"
0 3 "right"
1 1 "fluid"
$EndPhysicalNames
$Entities
5 4 0 0
1 0 0 0 1 2
2 0.25 0 0 0
3 0.5 0 0 0
4 0.75 0 0 0
5 1 0 0 1 3
1 0 0 0 0.25 0 0 1 1 2 1 -2
2 0.25 0 0 0.5 0 0 1 1 2 2 -3
3 0.5 0 0 0.75 0 0 1 1 2 4 -3
4 0.75 0 0 1 0 0 1 1 2 4 -5
$EndEntities
$Nodes
5 5 1 5
0 1 0 1
1
0 0 0
0 2 0 1
2
0.25 0 0
0 3 0 1
3
0.5 0 0
0 4 0 1
4
0.75 0 0
0 5 0 1
5
1 0 0
$EndNodes
$Elements
6 6 1 6
0 1 15 1
1 1
0 5 15 1
2 5
1 3 1 1
3 4 3
1 2 1 1
4 2 3
1 1 1 1
5 1 2
1 4 1 1
6 4 5
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
  const std::string path = write_fixture(kScrambledMsh);
  const auto doc = wasl::read_msh(path);
  std::remove(path.c_str());

  if (doc.nodes.size() != 5)
    return false;
  if (doc.elements.size() != 6)
    return false; // 4 lines + 2 points
  if (doc.physical_names.size() != 3)
    return false;
  if (doc.curve_entities.size() != 4)
    return false;
  if (doc.point_entities.size() != 5)
    return false;
  return true;
}

bool check_mesh_orientation_and_order() {
  const std::string path = write_fixture(kScrambledMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_line_mesh(doc);
  std::remove(path.c_str());

  if (mesh.num_elements() != 4)
    return false;

  // Elements must be sorted by increasing left endpoint x regardless of the
  // scrambled MSH record order. Curve 3 was written reversed (4 before 3); it
  // must still appear as x_left=0.5, x_right=0.75.
  const double expected_left[] = {0.0, 0.25, 0.5, 0.75};
  const double expected_right[] = {0.25, 0.5, 0.75, 1.0};
  for (std::size_t i = 0; i < mesh.num_elements(); ++i) {
    if (mesh.elements[i].x_left != expected_left[i])
      return false;
    if (mesh.elements[i].x_right != expected_right[i])
      return false;
    if (mesh.elements[i].left_node > mesh.elements[i].right_node)
      return false;
  }
  return true;
}

bool check_connectivity() {
  const std::string path = write_fixture(kScrambledMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_line_mesh(doc);
  std::remove(path.c_str());

  // Four elements share three interior endpoints and two boundaries: 5 faces.
  if (mesh.faces.size() != 5)
    return false;
  if (wasl::count_interior_faces(mesh) != 3)
    return false;
  if (wasl::count_boundary_faces(mesh) != 2)
    return false;

  // Left boundary: outward normal -1, no left neighbor (boundary has no right
  // element in Wasl's convention; here left_element is the single owner).
  const auto &left = mesh.faces.front();
  if (left.right_element != wasl::invalid_element)
    return false;
  if (left.outward_normal != -1.0)
    return false;
  if (left.physical_tag != 2)
    return false;

  // Right boundary: outward normal +1.
  const auto &right = mesh.faces.back();
  if (right.right_element != wasl::invalid_element)
    return false;
  if (right.outward_normal != 1.0)
    return false;
  if (right.physical_tag != 3)
    return false;

  if (mesh.physical_to_faces.at(2) != std::vector<std::int32_t>{0})
    return false;
  if (mesh.physical_to_faces.at(3) != std::vector<std::int32_t>{4})
    return false;

  // Interior faces pair consecutive elements.
  for (std::size_t i = 1; i + 1 < mesh.faces.size(); ++i) {
    const auto &f = mesh.faces[i];
    if (f.left_element + 1 != f.right_element)
      return false;
    if (f.outward_normal != 1.0)
      return false; // canonical interior direction
  }
  return true;
}

bool check_constant_geometry() {
  // Equal-length elements produce a uniform Jacobian. This is the constant-
  // geometry invariant: a structured sub-mesh maps to a single metric.
  const std::string path = write_fixture(kScrambledMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_line_mesh(doc);
  std::remove(path.c_str());

  for (std::size_t i = 0; i < mesh.num_elements(); ++i) {
    const auto map = wasl::element_map(mesh, i);
    if (map.length() != 0.25)
      return false;
    if (map.jacobian() != 0.125)
      return false;
  }
  return true;
}

bool check_physical_groups() {
  const std::string path = write_fixture(kScrambledMsh);
  const auto doc = wasl::read_msh(path);
  const auto mesh = wasl::build_line_mesh(doc);
  std::remove(path.c_str());

  // All four line elements belong to the "fluid" physical group (tag 1).
  if (mesh.physical_to_elements.at(1).size() != 4)
    return false;

  // Physical names resolve by tag.
  if (mesh.physical_names.at({1, 1}) != "fluid")
    return false;
  if (mesh.physical_names.at({0, 2}) != "left")
    return false;
  if (mesh.physical_names.at({0, 3}) != "right")
    return false;
  return true;
}

bool check_invalid_node_tag_error() {
  constexpr const char *invalid = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 1 1 1
0 1 0 1
not-a-node-tag
0 0 0
$EndNodes
)MSH";
  const std::string path = write_fixture(invalid);
  try {
    static_cast<void>(wasl::read_msh(path));
  } catch (const std::runtime_error &error) {
    std::remove(path.c_str());
    return std::string{error.what()}.find("node tag") != std::string::npos;
  }
  std::remove(path.c_str());
  return false;
}

bool check_missing_node_error() {
  wasl::MshDocument document;
  document.nodes.push_back(wasl::MshNode{1, 0, 1, 0.0, 0.0, 0.0});
  document.elements.push_back(wasl::MshElement{1, 1, 1, 1, {1, 2}});
  try {
    static_cast<void>(wasl::build_line_mesh(document));
  } catch (const std::runtime_error &error) {
    return std::string{error.what()}.find("missing node 2") != std::string::npos;
  }
  return false;
}

} // namespace

int main() {
  if (!check_parse_counts())
    return 1;
  if (!check_mesh_orientation_and_order())
    return 1;
  if (!check_connectivity())
    return 1;
  if (!check_constant_geometry())
    return 1;
  if (!check_physical_groups())
    return 1;
  if (!check_invalid_node_tag_error())
    return 1;
  if (!check_missing_node_error())
    return 1;
  return 0;
}
