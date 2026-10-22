#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "wasl/equations/euler/euler.hpp"

namespace wasl {

struct ScalarInterval {
  double lower;
  double upper;
  double value;
};

struct ScalarCircle {
  double center_x;
  double center_y;
  double radius;
  double value;
};

struct ScalarAdvectionCase {
  std::filesystem::path prepared_mesh;
  std::filesystem::path output_directory;
  std::filesystem::path restart_checkpoint;
  bool mpi = false;
  int dimensions = 1;
  int polynomial_order = 0;
  double wavespeed = 0.0;
  double velocity_y = 0.0;
  double cfl = 0.0;
  double start_time = 0.0;
  double end_time = 0.0;
  std::uint64_t max_steps = 0;
  std::uint64_t history_interval = 0;
  std::uint64_t visualization_interval = 0;
  std::uint64_t checkpoint_interval = 0;
  std::uint64_t keep_checkpoints = 0;
  double initial_default = 0.0;
  bool initial_sine = false;
  double initial_sine_amplitude = 0.0;
  std::vector<ScalarInterval> initial_intervals;
  std::vector<ScalarCircle> initial_circles;
};

struct EulerRegion {
  enum class Shape { default_region, interval, circle };
  Shape shape = Shape::default_region;
  double lower = 0.0;
  double upper = 0.0;
  double center_x = 0.0;
  double center_y = 0.0;
  double radius = 0.0;
  EulerPrimitive value{};
};

struct EulerCase {
  std::filesystem::path prepared_mesh;
  std::filesystem::path output_directory;
  std::filesystem::path restart_checkpoint;
  bool limiting = false;
  double sensor_threshold = 1e-3;
  double density_floor = 1e-10;
  double pressure_floor = 1e-10;
  int dimensions = 1;
  int polynomial_order = 0;
  double heat_capacity_ratio = 1.4;
  double cfl = 0.0;
  double start_time = 0.0;
  double end_time = 0.0;
  std::uint64_t max_steps = 0;
  std::uint64_t history_interval = 0;
  std::uint64_t visualization_interval = 0;
  std::uint64_t checkpoint_interval = 0;
  std::uint64_t keep_checkpoints = 0;
  std::vector<EulerRegion> initial_regions;
};

struct SimulationSummary {
  std::uint64_t step;
  double time;
  std::filesystem::path checkpoint;
};

SimulationSummary run_scalar_advection_case(const ScalarAdvectionCase &config);
SimulationSummary run_euler_case(const EulerCase &config);
std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &root,
                                                const std::string &effective_config);

} // namespace wasl
