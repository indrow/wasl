#include "wasl/mesh/msh/msh_reader.hpp"

#include <string>

#include "mesh/msh/msh_sections.hpp"
#include "mesh/msh/msh_token_reader.hpp"

namespace wasl {

MshDocument read_msh(const std::string &path) {
  MshTokenReader reader{path};
  MshDocument document;
  std::string section;
  while (reader.next(section)) {
    if (section == "$MeshFormat")
      parse_mesh_format(reader);
    else if (section == "$PhysicalNames")
      parse_physical_names(reader, document);
    else if (section == "$Entities")
      parse_entities(reader, document);
    else if (section == "$Nodes")
      parse_nodes(reader, document);
    else if (section == "$Elements")
      parse_elements(reader, document);
    else
      skip_section(reader, section);
  }
  return document;
}

} // namespace wasl
