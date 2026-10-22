#include "wasl/io/vtu_writer.hpp"

#include <ostream>
#include <stdexcept>

#include "io/text_output.hpp"
#include "wasl/mesh/prepared/prepared_mesh_validation.hpp"

namespace wasl {
namespace {

// VTK cell types: 3 is VTK_LINE, 5 is VTK_TRIANGLE.
int vtk_cell_type(int dimension) { return dimension == 1 ? 3 : 5; }

void write_points(std::ostream &output, const PreparedMesh &mesh) {
  output << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (double value : mesh.points)
    output << value << ' ';
  output << "\n</DataArray></Points>\n";
}

void write_cells(std::ostream &output, const PreparedMesh &mesh) {
  const auto cell_count = mesh.element_count();
  output << "<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto index : mesh.element_nodes) {
    if (index < 0 || static_cast<std::size_t>(index) >= mesh.node_tags.size())
      throw std::runtime_error("VTU output: element references invalid local node index");
    output << index << ' ';
  }
  output << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t cell = 1; cell <= cell_count; ++cell)
    output << cell * static_cast<std::size_t>(mesh.nodes_per_element) << ' ';
  output << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t cell = 0; cell < cell_count; ++cell)
    output << vtk_cell_type(mesh.dimension) << ' ';
  output << "\n</DataArray>\n</Cells>\n";
}

void write_cell_data(std::ostream &output, const std::vector<double> &cell_values,
                     std::size_t component_count) {
  output << "<CellData Scalars=\"state\"><DataArray type=\"Float64\" Name=\"state\" "
            "NumberOfComponents=\""
         << component_count << "\" format=\"ascii\">\n";
  for (double value : cell_values)
    output << value << ' ';
  output << "\n</DataArray></CellData>\n";
}

} // namespace

void write_vtu(const std::filesystem::path &path, const PreparedMesh &mesh,
               const std::vector<double> &cell_values, std::size_t component_count) {
  validate_prepared_mesh(mesh);
  if (component_count == 0 || cell_values.size() != mesh.element_count() * component_count)
    throw std::runtime_error("VTU output: cell-value size mismatch");
  auto output = open_text_output(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<UnstructuredGrid><Piece NumberOfPoints=\""
         << mesh.node_tags.size() << "\" NumberOfCells=\"" << mesh.element_count() << "\">\n";
  write_points(output, mesh);
  write_cells(output, mesh);
  write_cell_data(output, cell_values, component_count);
  output << "</Piece></UnstructuredGrid></VTKFile>\n";
}

} // namespace wasl
