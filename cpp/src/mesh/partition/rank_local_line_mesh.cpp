#include "wasl/mesh/partition/rank_local_line_mesh.hpp"

namespace wasl {

RankLocalLineMesh rank_local_line_mesh(const LineMesh &global,
                                       const std::vector<std::int32_t> &owned,
                                       const std::vector<std::int32_t> &ghosts) {
  RankLocalLineMesh result;
  std::int32_t local = 0;
  for (const auto *elements : {&owned, &ghosts}) {
    for (const auto global_element : *elements) {
      result.mesh.elements.push_back(global.elements.at(global_element));
      result.global_to_local[global_element] = local++;
    }
  }
  for (const auto &face : global.faces) {
    const bool left_local = result.global_to_local.contains(face.left_element);
    const bool right_local =
        face.is_boundary() || result.global_to_local.contains(face.right_element);
    if (!left_local || !right_local)
      continue;
    auto local_face = face;
    local_face.left_element = result.global_to_local.at(face.left_element);
    if (!face.is_boundary())
      local_face.right_element = result.global_to_local.at(face.right_element);
    result.mesh.faces.push_back(local_face);
  }
  return result;
}

} // namespace wasl
