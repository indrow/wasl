#include "wasl/mesh/msh_reader.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace wasl {
namespace {

class TokenReader {
public:
  explicit TokenReader(const std::string &path) : input_{path} {
    if (!input_) {
      throw std::runtime_error("MSH reader: cannot open " + path);
    }
  }

  template <class T> T read(const char *context) {
    T value{};
    if (!(input_ >> value)) {
      throw std::runtime_error("MSH reader: invalid or missing " + std::string{context});
    }
    return value;
  }

  std::string read_name() {
    std::string value;
    if (!(input_ >> std::quoted(value))) {
      throw std::runtime_error("MSH reader: invalid or missing physical name");
    }
    return value;
  }

  bool next(std::string &value) { return static_cast<bool>(input_ >> value); }

  void expect(const char *expected) {
    const auto token = read<std::string>(expected);
    if (token != expected) {
      throw std::runtime_error("MSH reader: expected " + std::string{expected} + ", got " + token);
    }
  }

  int count(const char *context) {
    const auto value = read<long long>(context);
    if (value < 0 || value > std::numeric_limits<int>::max()) {
      throw std::runtime_error("MSH reader: invalid " + std::string{context});
    }
    return static_cast<int>(value);
  }

private:
  std::ifstream input_;
};

std::vector<int> read_ints(TokenReader &reader, int count, const char *context) {
  std::vector<int> values;
  values.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    values.push_back(reader.read<int>(context));
  }
  return values;
}

void parse_mesh_format(TokenReader &reader) {
  if (reader.read<std::string>("mesh format version") != "4.1") {
    throw std::runtime_error("MSH reader: only version 4.1 is supported");
  }
  if (reader.read<int>("mesh file type") != 0) {
    throw std::runtime_error("MSH reader: only ASCII files are supported");
  }
  if (reader.read<int>("mesh data size") != 8) {
    throw std::runtime_error("MSH reader: expected data size 8");
  }
  reader.expect("$EndMeshFormat");
}

void parse_physical_names(TokenReader &reader, MshDocument &document) {
  const int count = reader.count("physical-name count");
  document.physical_names.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    document.physical_names.push_back(MshPhysicalName{reader.read<int>("physical-name dimension"),
                                                      reader.read<int>("physical-name tag"),
                                                      reader.read_name()});
  }
  reader.expect("$EndPhysicalNames");
}

void parse_entities(TokenReader &reader, MshDocument &document) {
  const int point_count = reader.count("point-entity count");
  const int curve_count = reader.count("curve-entity count");
  const int surface_count = reader.count("surface-entity count");
  const int volume_count = reader.count("volume-entity count");
  if (volume_count != 0) {
    throw std::runtime_error("MSH reader: volume entities are not supported");
  }

  document.point_entities.reserve(static_cast<std::size_t>(point_count));
  for (int i = 0; i < point_count; ++i) {
    MshPointEntity entity{};
    entity.tag = reader.read<int>("point-entity tag");
    entity.x = reader.read<double>("point-entity x");
    entity.y = reader.read<double>("point-entity y");
    entity.z = reader.read<double>("point-entity z");
    entity.physical_tags =
        read_ints(reader, reader.count("point physical-tag count"), "point physical tag");
    document.point_entities.push_back(std::move(entity));
  }

  document.curve_entities.reserve(static_cast<std::size_t>(curve_count));
  for (int i = 0; i < curve_count; ++i) {
    MshCurveEntity entity{};
    entity.tag = reader.read<int>("curve-entity tag");
    entity.min_x = reader.read<double>("curve minimum x");
    entity.min_y = reader.read<double>("curve minimum y");
    entity.min_z = reader.read<double>("curve minimum z");
    entity.max_x = reader.read<double>("curve maximum x");
    entity.max_y = reader.read<double>("curve maximum y");
    entity.max_z = reader.read<double>("curve maximum z");
    entity.physical_tags =
        read_ints(reader, reader.count("curve physical-tag count"), "curve physical tag");
    entity.bounding_points =
        read_ints(reader, reader.count("curve boundary-point count"), "curve boundary point");
    document.curve_entities.push_back(std::move(entity));
  }

  document.surface_entities.reserve(static_cast<std::size_t>(surface_count));
  for (int i = 0; i < surface_count; ++i) {
    MshSurfaceEntity entity{};
    entity.tag = reader.read<int>("surface-entity tag");
    entity.min_x = reader.read<double>("surface minimum x");
    entity.min_y = reader.read<double>("surface minimum y");
    entity.min_z = reader.read<double>("surface minimum z");
    entity.max_x = reader.read<double>("surface maximum x");
    entity.max_y = reader.read<double>("surface maximum y");
    entity.max_z = reader.read<double>("surface maximum z");
    entity.physical_tags =
        read_ints(reader, reader.count("surface physical-tag count"), "surface physical tag");
    entity.bounding_curves =
        read_ints(reader, reader.count("surface bounding-curve count"), "surface bounding curve");
    document.surface_entities.push_back(std::move(entity));
  }
  reader.expect("$EndEntities");
}

void parse_nodes(TokenReader &reader, MshDocument &document) {
  const int block_count = reader.count("node-block count");
  const int node_count = reader.count("node count");
  reader.read<std::int64_t>("minimum node tag");
  reader.read<std::int64_t>("maximum node tag");

  document.nodes.reserve(static_cast<std::size_t>(node_count));
  std::unordered_set<std::int64_t> seen_tags;
  for (int block = 0; block < block_count; ++block) {
    const int entity_dim = reader.read<int>("node-block entity dimension");
    const int entity_tag = reader.read<int>("node-block entity tag");
    if (reader.read<int>("node-block parametric flag") != 0) {
      throw std::runtime_error("MSH reader: parametric nodes are not supported");
    }
    const int count = reader.count("nodes in block");
    std::vector<std::int64_t> tags;
    tags.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      const auto tag = reader.read<std::int64_t>("node tag");
      if (!seen_tags.insert(tag).second) {
        throw std::runtime_error("MSH reader: duplicate node tag " + std::to_string(tag));
      }
      tags.push_back(tag);
    }
    for (const auto tag : tags) {
      document.nodes.push_back(MshNode{tag, entity_dim, entity_tag, reader.read<double>("node x"),
                                       reader.read<double>("node y"),
                                       reader.read<double>("node z")});
    }
  }
  if (document.nodes.size() != static_cast<std::size_t>(node_count)) {
    throw std::runtime_error("MSH reader: node count does not match block contents");
  }
  reader.expect("$EndNodes");
}

int nodes_per_element(int element_type) {
  if (element_type == 1) {
    return 2; // 2-node line
  }
  if (element_type == 2) {
    return 3; // 3-node triangle
  }
  if (element_type == 15) {
    return 1; // node
  }
  throw std::runtime_error("MSH reader: unsupported element type " + std::to_string(element_type));
}

void parse_elements(TokenReader &reader, MshDocument &document) {
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
      if (!seen_tags.insert(element.tag).second) {
        throw std::runtime_error("MSH reader: duplicate element tag " +
                                 std::to_string(element.tag));
      }
      element.entity_dim = entity_dim;
      element.entity_tag = entity_tag;
      element.element_type = element_type;
      element.nodes.reserve(static_cast<std::size_t>(node_count));
      for (int node = 0; node < node_count; ++node) {
        element.nodes.push_back(reader.read<std::int64_t>("element node tag"));
      }
      document.elements.push_back(std::move(element));
    }
  }
  if (document.elements.size() != static_cast<std::size_t>(element_count)) {
    throw std::runtime_error("MSH reader: element count does not match block contents");
  }
  reader.expect("$EndElements");
}

void skip_section(TokenReader &reader, const std::string &section) {
  if (section.empty() || section.front() != '$') {
    throw std::runtime_error("MSH reader: invalid section name " + section);
  }
  const std::string end = "$End" + section.substr(1);
  std::string token;
  while (reader.next(token)) {
    if (token == end) {
      return;
    }
  }
  throw std::runtime_error("MSH reader: missing " + end);
}

} // namespace

int MshDocument::first_physical_tag(int entity_dim, int entity_tag) const noexcept {
  if (entity_dim == 0) {
    for (const auto &entity : point_entities) {
      if (entity.tag == entity_tag) {
        return entity.physical_tags.empty() ? -1 : entity.physical_tags.front();
      }
    }
  } else if (entity_dim == 1) {
    for (const auto &entity : curve_entities) {
      if (entity.tag == entity_tag) {
        return entity.physical_tags.empty() ? -1 : entity.physical_tags.front();
      }
    }
  } else if (entity_dim == 2) {
    for (const auto &entity : surface_entities) {
      if (entity.tag == entity_tag) {
        return entity.physical_tags.empty() ? -1 : entity.physical_tags.front();
      }
    }
  }
  return -1;
}

MshDocument read_msh(const std::string &path) {
  TokenReader reader{path};
  MshDocument document;
  std::string section;
  while (reader.next(section)) {
    if (section == "$MeshFormat") {
      parse_mesh_format(reader);
    } else if (section == "$PhysicalNames") {
      parse_physical_names(reader, document);
    } else if (section == "$Entities") {
      parse_entities(reader, document);
    } else if (section == "$Nodes") {
      parse_nodes(reader, document);
    } else if (section == "$Elements") {
      parse_elements(reader, document);
    } else {
      skip_section(reader, section);
    }
  }
  return document;
}

} // namespace wasl
