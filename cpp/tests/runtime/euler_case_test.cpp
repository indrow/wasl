#include <cmath>
#include <filesystem>

#include "wasl/case/euler_case.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/io/checkpoint.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/prepared_mesh_file.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/prepared/prepare_mesh.hpp"
#include "wasl/mesh/triangle/build_triangle_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace {

wasl::EulerCase config(const std::filesystem::path &mesh, const std::filesystem::path &output,
                       int dimensions, bool limiting = false) {
  wasl::EulerCase result;
  result.files.prepared_mesh = mesh;
  result.files.output_directory = output;
  result.dimensions = dimensions;
  result.limiting = limiting;
  result.polynomial_order = limiting ? 2 : 1;
  result.cfl = 0.1;
  result.run.end_time = 0.01;
  result.run.max_steps = 100;
  result.run.history_interval = 100;
  result.run.visualization_interval = 100;
  result.run.checkpoint_interval = 100;
  result.run.keep_checkpoints = 1;
  result.initial_regions.push_back(
      {wasl::EulerRegion::Shape::default_region, 0.0, 0.0, 0.0, 0.0, 0.0,
       limiting ? wasl::EulerPrimitive{0.8, -0.1, dimensions == 2 ? 0.1 : 0.0, 0.7}
                : wasl::EulerPrimitive{1.1, 0.4, dimensions == 2 ? -0.2 : 0.0, 0.9}});
  if (limiting) {
    result.run.end_time = 0.001;
    result.initial_regions.push_back(
        {dimensions == 1 ? wasl::EulerRegion::Shape::interval : wasl::EulerRegion::Shape::circle,
         0.0, 0.5, 0.25, 0.5, 0.3, wasl::EulerPrimitive{1.2, 0.2, 0.0, 1.0}});
  }
  return result;
}

bool constant(const wasl::Checkpoint &checkpoint, const wasl::EulerConserved &state) {
  const double expected[] = {state.density, state.momentum_x, state.momentum_y, state.total_energy};
  for (std::size_t variable = 0; variable < checkpoint.variable_count; ++variable) {
    for (std::size_t element = 0; element < checkpoint.element_count; ++element) {
      const auto mode_zero =
          (variable * checkpoint.mode_count) * checkpoint.element_count + element;
      if (std::fabs(checkpoint.coefficients[mode_zero] - expected[variable] * std::sqrt(2.0)) >
          3e-11)
        return false;
      for (std::size_t mode = 1; mode < checkpoint.mode_count; ++mode)
        if (std::fabs(checkpoint.coefficients[mode_zero + mode * checkpoint.element_count]) > 3e-11)
          return false;
    }
  }
  return true;
}

} // namespace

int main() {
  const auto root = std::filesystem::temp_directory_path() / "wasl-euler-case-test";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  wasl::LineMesh line;
  line.elements = {{1, 2, 0.0, 0.4, 1}, {2, 3, 0.4, 1.0, 1}};
  line.faces = {{0, wasl::invalid_element, 0.0, -1.0, 2},
                {0, 1, 0.4, 1.0, wasl::invalid_physical},
                {1, wasl::invalid_element, 1.0, 1.0, 3}};
  const auto line_path = root / "line.h5";
  wasl::write_prepared_mesh(line_path, wasl::prepare_mesh(line));

  wasl::MshDocument document;
  document.surface_entities.push_back(
      wasl::MshSurfaceEntity{1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, {1}, {}});
  document.nodes = {{1, 2, 1, 0.0, 0.0, 0.0},
                    {2, 2, 1, 1.0, 0.0, 0.0},
                    {3, 2, 1, 1.0, 1.0, 0.0},
                    {4, 2, 1, 0.0, 1.0, 0.0}};
  document.elements = {{1, 2, 1, 2, {1, 2, 3}}, {2, 2, 1, 2, {1, 3, 4}}};
  const auto triangle_path = root / "triangle.h5";
  wasl::write_prepared_mesh(triangle_path, wasl::prepare_mesh(wasl::build_triangle_mesh(document)));

  const wasl::Euler equation{1.4};
  for (const int dimensions : {1, 2}) {
    const auto mesh = dimensions == 1 ? line_path : triangle_path;
    const auto run =
        wasl::run_euler_case(config(mesh, root / std::to_string(dimensions), dimensions));
    const auto checkpoint = wasl::read_checkpoint(run.checkpoint);
    const auto state = wasl::conserved(equation, {1.1, 0.4, dimensions == 2 ? -0.2 : 0.0, 0.9});
    if (!constant(checkpoint, state))
      return 1;
    const auto limited = wasl::run_euler_case(
        config(mesh, root / (std::to_string(dimensions) + "-limited"), dimensions, true));
    const auto limited_checkpoint = wasl::read_checkpoint(limited.checkpoint);
    for (std::size_t element = 0; element < limited_checkpoint.element_count; ++element) {
      const auto scale = 1.0 / std::sqrt(2.0);
      const wasl::EulerConserved average{
          scale * limited_checkpoint.coefficients[element],
          scale * limited_checkpoint.coefficients[limited_checkpoint.mode_count *
                                                      limited_checkpoint.element_count +
                                                  element],
          scale * limited_checkpoint.coefficients[2 * limited_checkpoint.mode_count *
                                                      limited_checkpoint.element_count +
                                                  element],
          scale * limited_checkpoint.coefficients[3 * limited_checkpoint.mode_count *
                                                      limited_checkpoint.element_count +
                                                  element]};
      if (!wasl::is_physical(equation, average))
        return 1;
    }
  }
  std::filesystem::remove_all(root);
  return 0;
}
