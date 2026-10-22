#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "wasl/mesh/line_geometry.hpp"

namespace wasl {

struct MshDocument;

constexpr std::int32_t invalid_element = -1;
constexpr int invalid_physical = -1;

struct LineElement {
  std::int64_t left_node;
  std::int64_t right_node;
  double x_left;
  double x_right;
  int physical_tag;
};

struct LineFace {
  std::int32_t left_element;
  std::int32_t right_element;
  double x;
  double outward_normal;
  int physical_tag;
};

struct LineMesh {
  std::vector<LineElement> elements;
  std::vector<LineFace> faces;
  std::map<std::pair<int, int>, std::string> physical_names;
  std::unordered_map<int, std::vector<std::int32_t>> physical_to_elements;
  std::unordered_map<int, std::vector<std::int32_t>> physical_to_faces;

  std::size_t num_elements() const noexcept { return elements.size(); }
};

LineMesh build_line_mesh(const MshDocument &document);
LineMap element_map(const LineMesh &mesh, std::size_t element);
std::size_t count_interior_faces(const LineMesh &mesh) noexcept;
std::size_t count_boundary_faces(const LineMesh &mesh) noexcept;

} // namespace wasl
