#include "wasl/mesh/line_mesh.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "wasl/mesh/msh_reader.hpp"

namespace wasl {
namespace {

using NodeLookup = std::unordered_map<std::int64_t, const MshNode *>;

NodeLookup index_nodes(const MshDocument &document) {
  NodeLookup nodes;
  nodes.reserve(document.nodes.size());
  for (const auto &node : document.nodes) {
    if (!nodes.emplace(node.tag, &node).second) {
      throw std::runtime_error("line mesh: duplicate node tag " + std::to_string(node.tag));
    }
  }
  return nodes;
}

const MshNode &require_node(const NodeLookup &nodes, std::int64_t tag) {
  const auto found = nodes.find(tag);
  if (found == nodes.end()) {
    throw std::runtime_error("line mesh: element references missing node " + std::to_string(tag));
  }
  return *found->second;
}

int node_physical_tag(const MshDocument &document, const MshNode &node) {
  return document.first_physical_tag(node.entity_dim, node.entity_tag);
}

} // namespace

LineMesh build_line_mesh(const MshDocument &document) {
  LineMesh mesh;
  for (const auto &physical : document.physical_names) {
    mesh.physical_names[{physical.dim, physical.tag}] = physical.name;
  }

  const auto nodes = index_nodes(document);
  for (const auto &element : document.elements) {
    if (element.element_type != 1 || element.nodes.size() != 2) {
      continue;
    }
    const auto &first = require_node(nodes, element.nodes[0]);
    const auto &second = require_node(nodes, element.nodes[1]);
    if (first.x == second.x) {
      throw std::runtime_error("line mesh: zero-length element " + std::to_string(element.tag));
    }

    const bool forward = first.x < second.x;
    mesh.elements.push_back(LineElement{
        forward ? first.tag : second.tag,
        forward ? second.tag : first.tag,
        forward ? first.x : second.x,
        forward ? second.x : first.x,
        document.first_physical_tag(element.entity_dim, element.entity_tag),
    });
  }
  if (mesh.elements.empty()) {
    throw std::runtime_error("line mesh: document contains no line elements");
  }

  std::stable_sort(
      mesh.elements.begin(), mesh.elements.end(),
      [](const LineElement &left, const LineElement &right) { return left.x_left < right.x_left; });

  for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
    mesh.physical_to_elements[mesh.elements[i].physical_tag].push_back(
        static_cast<std::int32_t>(i));
  }

  mesh.faces.reserve(mesh.elements.size() + 1);
  const auto &left_node = require_node(nodes, mesh.elements.front().left_node);
  mesh.faces.push_back(LineFace{0, invalid_element, mesh.elements.front().x_left, -1.0,
                                node_physical_tag(document, left_node)});

  for (std::size_t i = 0; i + 1 < mesh.elements.size(); ++i) {
    const auto &left = mesh.elements[i];
    const auto &right = mesh.elements[i + 1];
    if (left.right_node != right.left_node) {
      throw std::runtime_error("line mesh: gap or discontinuity between elements " +
                               std::to_string(i) + " and " + std::to_string(i + 1));
    }
    mesh.faces.push_back(LineFace{static_cast<std::int32_t>(i), static_cast<std::int32_t>(i + 1),
                                  left.x_right, 1.0, invalid_physical});
  }

  const auto last = static_cast<std::int32_t>(mesh.elements.size() - 1);
  const auto &right_node = require_node(nodes, mesh.elements.back().right_node);
  mesh.faces.push_back(LineFace{last, invalid_element, mesh.elements.back().x_right, 1.0,
                                node_physical_tag(document, right_node)});

  for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
    const int tag = mesh.faces[i].physical_tag;
    if (tag != invalid_physical) {
      mesh.physical_to_faces[tag].push_back(static_cast<std::int32_t>(i));
    }
  }
  return mesh;
}

LineMap element_map(const LineMesh &mesh, std::size_t element) {
  const auto &value = mesh.elements.at(element);
  return LineMap{value.x_left, value.x_right};
}

std::size_t count_interior_faces(const LineMesh &mesh) noexcept {
  return static_cast<std::size_t>(
      std::count_if(mesh.faces.begin(), mesh.faces.end(),
                    [](const LineFace &face) { return face.right_element != invalid_element; }));
}

std::size_t count_boundary_faces(const LineMesh &mesh) noexcept {
  return mesh.faces.size() - count_interior_faces(mesh);
}

} // namespace wasl
