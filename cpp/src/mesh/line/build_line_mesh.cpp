#include "wasl/mesh/line/build_line_mesh.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "mesh/msh/msh_node_lookup.hpp"

namespace wasl {
namespace {

// Line elements oriented left to right and sorted by position.
std::vector<LineElement> oriented_elements(const MshDocument &document,
                                           const MshNodeLookup &nodes) {
  std::vector<LineElement> elements;
  for (const auto &element : document.elements) {
    if (element.element_type != msh_line_element || element.nodes.size() != 2)
      continue;
    const auto &first = nodes.at(element.nodes[0]);
    const auto &second = nodes.at(element.nodes[1]);
    if (first.x == second.x)
      throw std::runtime_error("line mesh: zero-length element " + std::to_string(element.tag));

    const bool forward = first.x < second.x;
    const auto &left = forward ? first : second;
    const auto &right = forward ? second : first;
    elements.push_back(
        LineElement{left.tag, right.tag, left.x, right.x,
                    document.first_physical_tag(element.entity_dim, element.entity_tag)});
  }
  if (elements.empty())
    throw std::runtime_error("line mesh: document contains no line elements");

  std::stable_sort(elements.begin(), elements.end(),
                   [](const LineElement &a, const LineElement &b) { return a.x_left < b.x_left; });
  return elements;
}

int node_physical_tag(const MshDocument &document, const MshNode &node) {
  return document.first_physical_tag(node.entity_dim, node.entity_tag);
}

// Left boundary, one face per adjacent pair, right boundary.
std::vector<LineFace> endpoint_faces(const std::vector<LineElement> &elements,
                                     const MshDocument &document, const MshNodeLookup &nodes) {
  std::vector<LineFace> faces;
  faces.reserve(elements.size() + 1);
  faces.push_back(LineFace{0, invalid_element, elements.front().x_left, LineMap::left_face_normal,
                           node_physical_tag(document, nodes.at(elements.front().left_node))});

  for (std::size_t i = 0; i + 1 < elements.size(); ++i) {
    if (elements[i].right_node != elements[i + 1].left_node)
      throw std::runtime_error("line mesh: gap or discontinuity between elements " +
                               std::to_string(i) + " and " + std::to_string(i + 1));
    faces.push_back(LineFace{static_cast<std::int32_t>(i), static_cast<std::int32_t>(i + 1),
                             elements[i].x_right, LineMap::right_face_normal, invalid_physical});
  }

  faces.push_back(LineFace{static_cast<std::int32_t>(elements.size() - 1), invalid_element,
                           elements.back().x_right, LineMap::right_face_normal,
                           node_physical_tag(document, nodes.at(elements.back().right_node))});
  return faces;
}

} // namespace

LineMesh build_line_mesh(const MshDocument &document) {
  const MshNodeLookup nodes(document, "line mesh");
  LineMesh mesh;
  mesh.physical_names = physical_names_of(document);
  mesh.elements = oriented_elements(document, nodes);
  mesh.faces = endpoint_faces(mesh.elements, document, nodes);

  for (std::size_t i = 0; i < mesh.elements.size(); ++i)
    mesh.physical_to_elements[mesh.elements[i].physical_tag].push_back(
        static_cast<std::int32_t>(i));
  for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
    const int tag = mesh.faces[i].physical_tag;
    if (tag != invalid_physical)
      mesh.physical_to_faces[tag].push_back(static_cast<std::int32_t>(i));
  }
  return mesh;
}

} // namespace wasl
