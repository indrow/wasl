#include "wasl/mesh/msh/msh_document.hpp"

#include <vector>

#include "wasl/mesh/mesh_index.hpp"

namespace wasl {
namespace {

template <typename Entity>
int first_physical_tag_of(const std::vector<Entity> &entities, int entity_tag) noexcept {
  for (const auto &entity : entities) {
    if (entity.tag == entity_tag)
      return entity.physical_tags.empty() ? invalid_physical : entity.physical_tags.front();
  }
  return invalid_physical;
}

} // namespace

int MshDocument::first_physical_tag(int entity_dim, int entity_tag) const noexcept {
  switch (entity_dim) {
  case 0:
    return first_physical_tag_of(point_entities, entity_tag);
  case 1:
    return first_physical_tag_of(curve_entities, entity_tag);
  case 2:
    return first_physical_tag_of(surface_entities, entity_tag);
  default:
    return invalid_physical;
  }
}

} // namespace wasl
