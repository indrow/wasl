#include "wasl/case/case_mesh.hpp"

#include <stdexcept>

#include "wasl/io/prepared_mesh_file.hpp"

namespace wasl {

CaseMesh load_case_mesh(const std::filesystem::path &prepared_mesh, int dimensions) {
  CaseMesh mesh;
  mesh.path = std::filesystem::weakly_canonical(prepared_mesh);
  mesh.prepared = read_prepared_mesh(mesh.path);
  if (mesh.prepared.dimension != dimensions)
    throw std::runtime_error("simulation: case dimension does not match prepared mesh");
  return mesh;
}

} // namespace wasl
