#include "wasl/io/pvtu_writer.hpp"

#include <stdexcept>

#include "io/text_output.hpp"

namespace wasl {

void write_pvtu(const std::filesystem::path &path, const std::vector<std::string> &pieces,
                std::size_t component_count) {
  if (pieces.empty() || component_count == 0)
    throw std::runtime_error("PVTU output: pieces and components must be non-empty");
  auto output = open_text_output(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n<PUnstructuredGrid GhostLevel=\"0\">\n"
            "<PCellData Scalars=\"state\"><PDataArray type=\"Float64\" Name=\"state\" "
            "NumberOfComponents=\""
         << component_count
         << "\"/></PCellData>\n"
            "<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n";
  for (const auto &piece : pieces)
    output << "<Piece Source=\"" << xml_escape(piece) << "\"/>\n";
  output << "</PUnstructuredGrid></VTKFile>\n";
}

} // namespace wasl
