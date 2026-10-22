#include "wasl/io/pvd_writer.hpp"

#include "io/text_output.hpp"

namespace wasl {

void write_pvd(const std::filesystem::path &path,
               const std::vector<std::pair<double, std::string>> &snapshots) {
  auto output = open_text_output(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n<Collection>\n";
  for (const auto &[time, file] : snapshots)
    output << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\""
           << xml_escape(file) << "\"/>\n";
  output << "</Collection></VTKFile>\n";
}

} // namespace wasl
