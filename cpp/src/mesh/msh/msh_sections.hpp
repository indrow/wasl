// Parsers for individual MSH 4.1 sections. Each consumes the body after the
// section header through its closing $End marker.

#pragma once

#include <string>

#include "mesh/msh/msh_token_reader.hpp"
#include "wasl/mesh/msh/msh_document.hpp"

namespace wasl {

void parse_mesh_format(MshTokenReader &reader);
void parse_physical_names(MshTokenReader &reader, MshDocument &document);
void parse_entities(MshTokenReader &reader, MshDocument &document);
void parse_nodes(MshTokenReader &reader, MshDocument &document);
void parse_elements(MshTokenReader &reader, MshDocument &document);
void skip_section(MshTokenReader &reader, const std::string &section);

} // namespace wasl
