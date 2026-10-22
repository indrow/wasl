#include "case/scalar_advection_drivers.hpp"

#include "wasl/case/scalar_coefficient_layout.hpp"
#include "wasl/case/scalar_diagnostics.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// The scalar solver keeps U[element][mode]; checkpoints persist
// U[variable][mode][element].
template <typename Mesh>
void add_scalar_state_operations(SerialCaseModel &model, const Mesh &mesh) {
  const auto element_count = mesh.num_elements();
  const auto mode_count = model.mode_count;
  model.from_checkpoint = [=](const SerialCaseModel::State &stored) {
    return solver_coefficients(stored, element_count, mode_count);
  };
  model.to_checkpoint = [=](const SerialCaseModel::State &state) {
    return checkpoint_coefficients(state, element_count, mode_count);
  };
  model.mass = [&mesh, mode_count](const SerialCaseModel::State &state) {
    return scalar_mass(mesh, state, mode_count);
  };
  model.cell_averages = [=](const SerialCaseModel::State &state) {
    return scalar_cell_averages(state, element_count, mode_count);
  };
}

template void add_scalar_state_operations(SerialCaseModel &, const LineMesh &);
template void add_scalar_state_operations(SerialCaseModel &, const TriangleMesh &);

} // namespace wasl
