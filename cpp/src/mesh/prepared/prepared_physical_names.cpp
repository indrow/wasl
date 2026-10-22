#include "wasl/mesh/prepared/prepared_physical_names.hpp"

#include <cstddef>
#include <string>

namespace wasl {

PhysicalNames prepared_physical_names(const PreparedMesh &mesh) {
  PhysicalNames names;
  for (std::size_t index = 0; index < mesh.physical_name_tags.size(); ++index) {
    const auto begin = mesh.physical_name_bytes.begin() +
                       static_cast<std::ptrdiff_t>(mesh.physical_name_offsets[index]);
    const auto end = mesh.physical_name_bytes.begin() +
                     static_cast<std::ptrdiff_t>(mesh.physical_name_offsets[index + 1]);
    names[{mesh.physical_name_dimensions[index], mesh.physical_name_tags[index]}] =
        std::string(begin, end);
  }
  return names;
}

} // namespace wasl
