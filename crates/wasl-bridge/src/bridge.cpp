#include "wasl-bridge/include/bridge.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/core/simulation.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

#ifdef WASL_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::string string(rust::Str value) { return std::string(value.data(), value.size()); }

} // namespace

rust::String prepare_mesh(rust::Str source_path, rust::Str output_path) {
  const auto source = string(source_path);
  const auto output = string(output_path);
  const auto document = wasl::read_msh(source);
  const bool has_triangles =
      std::any_of(document.elements.begin(), document.elements.end(),
                  [](const auto &element) { return element.element_type == 2; });
  const auto prepared = has_triangles ? wasl::prepare_mesh(wasl::build_triangle_mesh(document))
                                      : wasl::prepare_mesh(wasl::build_line_mesh(document));
  wasl::write_prepared_mesh(output, prepared);
  return "prepared " + std::to_string(prepared.element_physical_tags.size()) + " elements";
}

rust::String inspect_checkpoint(rust::Str checkpoint_path) {
  const auto checkpoint = wasl::read_checkpoint(string(checkpoint_path));
  std::ostringstream summary;
  summary.precision(17);
  summary << "step=" << checkpoint.step << " time=" << checkpoint.time
          << " variables=" << checkpoint.variable_count << " modes=" << checkpoint.mode_count
          << " elements=" << checkpoint.element_count;
  return summary.str();
}

rust::String inspect_prepared_mesh(rust::Str mesh_path) {
  const auto mesh = wasl::read_prepared_mesh(string(mesh_path));
  return std::to_string(mesh.dimension) + "D mesh with " +
         std::to_string(mesh.element_physical_tags.size()) + " elements";
}

rust::String export_checkpoint(rust::Str checkpoint_path, rust::Str mesh_path,
                               rust::Str output_path) {
  const auto checkpoint = wasl::read_checkpoint(string(checkpoint_path));
  const auto mesh = wasl::read_prepared_mesh(string(mesh_path));
  if (checkpoint.element_count != mesh.element_physical_tags.size()) {
    throw std::runtime_error("checkpoint export: mesh element count mismatch");
  }
  // Mode zero is phi_0 = 1/sqrt(2), so the cell average is c_0 / sqrt(2).
  const double inverse_sqrt2 = 1.0 / std::sqrt(2.0);
  std::vector<double> averages(checkpoint.variable_count * checkpoint.element_count);
  for (std::size_t variable = 0; variable < checkpoint.variable_count; ++variable) {
    for (std::size_t element = 0; element < checkpoint.element_count; ++element) {
      const auto source = (variable * checkpoint.mode_count) * checkpoint.element_count + element;
      averages[element * checkpoint.variable_count + variable] =
          checkpoint.coefficients[source] * inverse_sqrt2;
    }
  }
  wasl::write_vtu(string(output_path), mesh, averages, checkpoint.variable_count);
  return "exported " + std::to_string(checkpoint.element_count) + " cells";
}

rust::String prepare_mpi_run_directory(rust::Str root_path, rust::Str effective_config) {
  return wasl::prepare_mpi_run_directory(string(root_path), string(effective_config)).string();
}

rust::String
run_scalar_advection_case(rust::Str mesh_path, rust::Str output_path, rust::Str restart_path,
                          bool mpi, int dimensions, int polynomial_order, double velocity_x,
                          double velocity_y, double cfl, double start_time, double end_time,
                          std::uint64_t max_steps, std::uint64_t history_interval,
                          std::uint64_t visualization_interval, std::uint64_t checkpoint_interval,
                          std::uint64_t keep_checkpoints, double initial_default, bool initial_sine,
                          double initial_sine_amplitude, rust::Slice<const double> intervals,
                          rust::Slice<const double> circles) {
  if (intervals.size() % 3 != 0) {
    throw std::runtime_error("initial.regions: invalid interval data");
  }
  if (circles.size() % 4 != 0) {
    throw std::runtime_error("initial.regions: invalid circle data");
  }
  wasl::ScalarAdvectionCase config;
  config.prepared_mesh = string(mesh_path);
  config.output_directory = string(output_path);
  config.restart_checkpoint = string(restart_path);
  config.mpi = mpi;
  config.dimensions = dimensions;
  config.polynomial_order = polynomial_order;
  config.wavespeed = velocity_x;
  config.velocity_y = velocity_y;
  config.cfl = cfl;
  config.start_time = start_time;
  config.end_time = end_time;
  config.max_steps = max_steps;
  config.history_interval = history_interval;
  config.visualization_interval = visualization_interval;
  config.checkpoint_interval = checkpoint_interval;
  config.keep_checkpoints = keep_checkpoints;
  config.initial_default = initial_default;
  config.initial_sine = initial_sine;
  config.initial_sine_amplitude = initial_sine_amplitude;
  for (std::size_t index = 0; index < intervals.size(); index += 3) {
    config.initial_intervals.push_back(
        {intervals[index], intervals[index + 1], intervals[index + 2]});
  }
  for (std::size_t index = 0; index < circles.size(); index += 4) {
    config.initial_circles.push_back(
        {circles[index], circles[index + 1], circles[index + 2], circles[index + 3]});
  }
  const auto summary = wasl::run_scalar_advection_case(config);
#ifdef WASL_ENABLE_MPI
  if (mpi) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
      return "";
    }
  }
#endif
  std::ostringstream message;
  message.precision(17);
  message << "step=" << summary.step << " time=" << summary.time
          << " checkpoint=" << summary.checkpoint.string();
  return message.str();
}

rust::String run_euler_case(rust::Str mesh_path, rust::Str output_path, rust::Str restart_path,
                            bool limiting, double sensor_threshold, double density_floor,
                            double pressure_floor, int dimensions, int polynomial_order,
                            double heat_capacity_ratio, double cfl, double start_time,
                            double end_time, std::uint64_t max_steps,
                            std::uint64_t history_interval, std::uint64_t visualization_interval,
                            std::uint64_t checkpoint_interval, std::uint64_t keep_checkpoints,
                            rust::Slice<const double> regions) {
  constexpr std::size_t region_width = 10;
  if (regions.size() % region_width != 0)
    throw std::runtime_error("initial.regions: invalid Euler region data");
  wasl::EulerCase config;
  config.prepared_mesh = string(mesh_path);
  config.output_directory = string(output_path);
  config.restart_checkpoint = string(restart_path);
  config.limiting = limiting;
  config.sensor_threshold = sensor_threshold;
  config.density_floor = density_floor;
  config.pressure_floor = pressure_floor;
  config.dimensions = dimensions;
  config.polynomial_order = polynomial_order;
  config.heat_capacity_ratio = heat_capacity_ratio;
  config.cfl = cfl;
  config.start_time = start_time;
  config.end_time = end_time;
  config.max_steps = max_steps;
  config.history_interval = history_interval;
  config.visualization_interval = visualization_interval;
  config.checkpoint_interval = checkpoint_interval;
  config.keep_checkpoints = keep_checkpoints;
  for (std::size_t index = 0; index < regions.size(); index += region_width) {
    wasl::EulerRegion region;
    const auto shape = static_cast<int>(regions[index]);
    if (shape < 0 || shape > 2)
      throw std::runtime_error("initial.regions: invalid Euler region shape");
    region.shape = static_cast<wasl::EulerRegion::Shape>(shape);
    region.lower = regions[index + 1];
    region.upper = regions[index + 2];
    region.center_x = regions[index + 3];
    region.center_y = regions[index + 4];
    region.radius = regions[index + 5];
    region.value = {regions[index + 6], regions[index + 7], regions[index + 8], regions[index + 9]};
    config.initial_regions.push_back(region);
  }
  const auto summary = wasl::run_euler_case(config);
  std::ostringstream message;
  message.precision(17);
  message << "step=" << summary.step << " time=" << summary.time
          << " checkpoint=" << summary.checkpoint.string();
  return message.str();
}
