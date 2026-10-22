// Text-file helpers shared by the CSV and VTK writers.

#pragma once

#include <filesystem>
#include <fstream>
#include <string>

namespace wasl {

void create_parent_directories(const std::filesystem::path &path);

// Open for writing with round-trip double precision; throws on failure.
std::ofstream open_text_output(const std::filesystem::path &path);

// Escape a string for an XML attribute value.
std::string xml_escape(const std::string &value);

} // namespace wasl
