// An irregular unit-interval line mesh whose ends are the physical points
// "left" and "right", shared by the case-level boundary tests.

#pragma once

#include <cmath>

#include "wasl/mesh/line/build_line_mesh.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/msh/msh_document.hpp"

namespace wasl_test {

inline wasl::LineMesh named_line(int count) {
  constexpr int left_tag = 2;
  constexpr int right_tag = 3;
  wasl::MshDocument document;
  document.physical_names = {{0, left_tag, "left"}, {0, right_tag, "right"}, {1, 1, "fluid"}};
  document.point_entities = {{1, 0.0, 0.0, 0.0, {left_tag}}, {2, 1.0, 0.0, 0.0, {right_tag}}};
  document.curve_entities = {{1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, {1}, {1, -2}}};
  for (int node = 0; node <= count; ++node) {
    double x = static_cast<double>(node) / count;
    int entity_dim = 1, entity_tag = 1;
    if (node == 0) {
      entity_dim = 0;
    } else if (node == count) {
      entity_dim = 0;
      entity_tag = 2;
    } else {
      x += 0.3 / count * std::sin(5.0 * node);
    }
    document.nodes.push_back({node + 1, entity_dim, entity_tag, x, 0.0, 0.0});
  }
  for (int element = 0; element < count; ++element)
    document.elements.push_back(
        {element + 1, 1, 1, wasl::msh_line_element, {element + 1, element + 2}});
  return wasl::build_line_mesh(document);
}

} // namespace wasl_test
