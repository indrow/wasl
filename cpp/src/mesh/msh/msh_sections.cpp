#include "mesh/msh/msh_sections.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace wasl {
namespace {

// Entity bounding boxes share one record shape across curves and surfaces.
template <typename Entity>
void read_bounding_box(MshTokenReader &reader, Entity &entity, const std::string &kind) {
  entity.min_x = reader.read<double>((kind + " minimum x").c_str());
  entity.min_y = reader.read<double>((kind + " minimum y").c_str());
  entity.min_z = reader.read<double>((kind + " minimum z").c_str());
  entity.max_x = reader.read<double>((kind + " maximum x").c_str());
  entity.max_y = reader.read<double>((kind + " maximum y").c_str());
  entity.max_z = reader.read<double>((kind + " maximum z").c_str());
}

MshPointEntity read_point_entity(MshTokenReader &reader) {
  MshPointEntity entity{};
  entity.tag = reader.read<int>("point-entity tag");
  entity.x = reader.read<double>("point-entity x");
  entity.y = reader.read<double>("point-entity y");
  entity.z = reader.read<double>("point-entity z");
  entity.physical_tags =
      reader.read_ints(reader.count("point physical-tag count"), "point physical tag");
  return entity;
}

MshCurveEntity read_curve_entity(MshTokenReader &reader) {
  MshCurveEntity entity{};
  entity.tag = reader.read<int>("curve-entity tag");
  read_bounding_box(reader, entity, "curve");
  entity.physical_tags =
      reader.read_ints(reader.count("curve physical-tag count"), "curve physical tag");
  entity.bounding_points =
      reader.read_ints(reader.count("curve boundary-point count"), "curve boundary point");
  return entity;
}

MshSurfaceEntity read_surface_entity(MshTokenReader &reader) {
  MshSurfaceEntity entity{};
  entity.tag = reader.read<int>("surface-entity tag");
  read_bounding_box(reader, entity, "surface");
  entity.physical_tags =
      reader.read_ints(reader.count("surface physical-tag count"), "surface physical tag");
  entity.bounding_curves =
      reader.read_ints(reader.count("surface bounding-curve count"), "surface bounding curve");
  return entity;
}

int nodes_per_element(int element_type) {
  switch (element_type) {
  case msh_line_element:
    return 2;
  case msh_triangle_element:
    return 3;
  case msh_point_element:
    return 1;
  default:
    throw std::runtime_error("MSH reader: unsupported element type " +
                             std::to_string(element_type));
  }
}

} // namespace

void parse_mesh_format(MshTokenReader &reader) {
  if (reader.read<std::string>("mesh format version") != "4.1")
    throw std::runtime_error("MSH reader: only version 4.1 is supported");
  if (reader.read<int>("mesh file type") != 0)
    throw std::runtime_error("MSH reader: only ASCII files are supported");
  if (reader.read<int>("mesh data size") != 8)
    throw std::runtime_error("MSH reader: expected data size 8");
  reader.expect("$EndMeshFormat");
}

void parse_physical_names(MshTokenReader &reader, MshDocument &document) {
  const int count = reader.count("physical-name count");
  document.physical_names.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const int dimension = reader.read<int>("physical-name dimension");
    const int tag = reader.read<int>("physical-name tag");
    document.physical_names.push_back(
        MshPhysicalName{dimension, tag, reader.read_quoted("physical name")});
  }
  reader.expect("$EndPhysicalNames");
}

void parse_entities(MshTokenReader &reader, MshDocument &document) {
  const int point_count = reader.count("point-entity count");
  const int curve_count = reader.count("curve-entity count");
  const int surface_count = reader.count("surface-entity count");
  if (reader.count("volume-entity count") != 0)
    throw std::runtime_error("MSH reader: volume entities are not supported");

  document.point_entities.reserve(static_cast<std::size_t>(point_count));
  for (int i = 0; i < point_count; ++i)
    document.point_entities.push_back(read_point_entity(reader));
  document.curve_entities.reserve(static_cast<std::size_t>(curve_count));
  for (int i = 0; i < curve_count; ++i)
    document.curve_entities.push_back(read_curve_entity(reader));
  document.surface_entities.reserve(static_cast<std::size_t>(surface_count));
  for (int i = 0; i < surface_count; ++i)
    document.surface_entities.push_back(read_surface_entity(reader));
  reader.expect("$EndEntities");
}

void parse_nodes(MshTokenReader &reader, MshDocument &document) {
  const int block_count = reader.count("node-block count");
  const int node_count = reader.count("node count");
  reader.read<std::int64_t>("minimum node tag");
  reader.read<std::int64_t>("maximum node tag");

  document.nodes.reserve(static_cast<std::size_t>(node_count));
  std::unordered_set<std::int64_t> seen_tags;
  for (int block = 0; block < block_count; ++block) {
    const int entity_dim = reader.read<int>("node-block entity dimension");
    const int entity_tag = reader.read<int>("node-block entity tag");
    if (reader.read<int>("node-block parametric flag") != 0)
      throw std::runtime_error("MSH reader: parametric nodes are not supported");
    const int count = reader.count("nodes in block");

    // A block lists all node tags first, then all coordinates.
    std::vector<std::int64_t> tags;
    tags.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      const auto tag = reader.read<std::int64_t>("node tag");
      if (!seen_tags.insert(tag).second)
        throw std::runtime_error("MSH reader: duplicate node tag " + std::to_string(tag));
      tags.push_back(tag);
    }
    for (const auto tag : tags) {
      const double x = reader.read<double>("node x");
      const double y = reader.read<double>("node y");
      const double z = reader.read<double>("node z");
      document.nodes.push_back(MshNode{tag, entity_dim, entity_tag, x, y, z});
    }
  }
  if (document.nodes.size() != static_cast<std::size_t>(node_count))
    throw std::runtime_error("MSH reader: node count does not match block contents");
  reader.expect("$EndNodes");
}

void parse_elements(MshTokenReader &reader, MshDocument &document) {
  const int block_count = reader.count("element-block count");
  const int element_count = reader.count("element count");
  reader.read<std::int64_t>("minimum element tag");
  reader.read<std::int64_t>("maximum element tag");

  document.elements.reserve(static_cast<std::size_t>(element_count));
  std::unordered_set<std::int64_t> seen_tags;
  for (int block = 0; block < block_count; ++block) {
    const int entity_dim = reader.read<int>("element-block entity dimension");
    const int entity_tag = reader.read<int>("element-block entity tag");
    const int element_type = reader.read<int>("element type");
    const int count = reader.count("elements in block");
    const int node_count = nodes_per_element(element_type);
    for (int i = 0; i < count; ++i) {
      MshElement element{};
      element.tag = reader.read<std::int64_t>("element tag");
      if (!seen_tags.insert(element.tag).second)
        throw std::runtime_error("MSH reader: duplicate element tag " +
                                 std::to_string(element.tag));
      element.entity_dim = entity_dim;
      element.entity_tag = entity_tag;
      element.element_type = element_type;
      element.nodes.reserve(static_cast<std::size_t>(node_count));
      for (int node = 0; node < node_count; ++node)
        element.nodes.push_back(reader.read<std::int64_t>("element node tag"));
      document.elements.push_back(std::move(element));
    }
  }
  if (document.elements.size() != static_cast<std::size_t>(element_count))
    throw std::runtime_error("MSH reader: element count does not match block contents");
  reader.expect("$EndElements");
}

void skip_section(MshTokenReader &reader, const std::string &section) {
  if (section.empty() || section.front() != '$')
    throw std::runtime_error("MSH reader: invalid section name " + section);
  const std::string end = "$End" + section.substr(1);
  std::string token;
  while (reader.next(token)) {
    if (token == end)
      return;
  }
  throw std::runtime_error("MSH reader: missing " + end);
}

} // namespace wasl
