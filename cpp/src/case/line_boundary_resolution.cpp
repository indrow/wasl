#include "case/line_boundary_resolution.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>

#include "wasl/mesh/prepared/prepared_physical_names.hpp"

namespace wasl {
namespace {

// Physical points tag the ends of a line mesh.
constexpr int boundary_dimension = 0;

std::string boundary_face_name(const PhysicalNames &names, const LineFace &face) {
  const auto found = names.find({boundary_dimension, face.physical_tag});
  if (face.physical_tag == invalid_physical || found == names.end())
    throw std::runtime_error("boundaries: the line boundary face at x = " + std::to_string(face.x) +
                             " has no physical name");
  return found->second;
}

const NamedBoundary &configured(const NamedBoundaries &boundaries, const std::string &name) {
  const auto found = std::find_if(boundaries.begin(), boundaries.end(),
                                  [&](const NamedBoundary &value) { return value.name == name; });
  if (found == boundaries.end())
    throw std::runtime_error("boundaries." + name + ": the mesh boundary has no condition");
  return *found;
}

} // namespace

bool all_periodic(const NamedBoundaries &boundaries) {
  return std::all_of(boundaries.begin(), boundaries.end(), [](const NamedBoundary &value) {
    return value.condition.kind == BoundaryKind::periodic;
  });
}

LineBoundaryConditions resolve_line_boundaries(const PreparedMesh &prepared, const LineMesh &mesh,
                                               const NamedBoundaries &boundaries) {
  if (boundaries.empty())
    return {};
  const auto names = prepared_physical_names(prepared);
  LineBoundaryConditions conditions;
  std::set<std::string> used;
  for (const auto &face : mesh.faces) {
    if (!face.is_boundary())
      continue;
    const auto name = boundary_face_name(names, face);
    used.insert(name);
    (face.outward_normal < 0.0 ? conditions.left : conditions.right) =
        configured(boundaries, name).condition;
  }
  for (const auto &boundary : boundaries) {
    if (!used.contains(boundary.name))
      throw std::runtime_error("boundaries." + boundary.name +
                               ": no line boundary face has this physical name");
  }
  if ((conditions.left.kind == BoundaryKind::periodic) !=
      (conditions.right.kind == BoundaryKind::periodic))
    throw std::runtime_error("boundaries: periodic conditions must be applied to both line ends");
  return conditions;
}

} // namespace wasl
