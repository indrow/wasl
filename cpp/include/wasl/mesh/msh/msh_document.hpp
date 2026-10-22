// In-memory image of the Gmsh MSH 4.1 sections Wasl consumes.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace wasl {

struct MshPhysicalName {
  int dim;
  int tag;
  std::string name;
};

struct MshPointEntity {
  int tag;
  double x;
  double y;
  double z;
  std::vector<int> physical_tags;
};

struct MshCurveEntity {
  int tag;
  double min_x;
  double min_y;
  double min_z;
  double max_x;
  double max_y;
  double max_z;
  std::vector<int> physical_tags;
  std::vector<int> bounding_points;
};

struct MshSurfaceEntity {
  int tag;
  double min_x;
  double min_y;
  double min_z;
  double max_x;
  double max_y;
  double max_z;
  std::vector<int> physical_tags;
  std::vector<int> bounding_curves;
};

struct MshNode {
  std::int64_t tag;
  int entity_dim;
  int entity_tag;
  double x;
  double y;
  double z;
};

// Gmsh element type codes Wasl reads.
constexpr int msh_line_element = 1;
constexpr int msh_triangle_element = 2;
constexpr int msh_point_element = 15;

struct MshElement {
  std::int64_t tag;
  int entity_dim;
  int entity_tag;
  int element_type;
  std::vector<std::int64_t> nodes;
};

struct MshDocument {
  std::vector<MshPhysicalName> physical_names;
  std::vector<MshPointEntity> point_entities;
  std::vector<MshCurveEntity> curve_entities;
  std::vector<MshSurfaceEntity> surface_entities;
  std::vector<MshNode> nodes;
  std::vector<MshElement> elements;

  // First physical tag of an entity, or invalid_physical when it has none.
  int first_physical_tag(int entity_dim, int entity_tag) const noexcept;
};

} // namespace wasl
