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

  int first_physical_tag(int entity_dim, int entity_tag) const noexcept;
};

MshDocument read_msh(const std::string &path);

} // namespace wasl
