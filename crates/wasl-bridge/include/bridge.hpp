#pragma once

#include <cstdint>

#include "rust/cxx.h"

rust::String prepare_mesh(rust::Str source_path, rust::Str output_path);
rust::String inspect_prepared_mesh(rust::Str mesh_path);
rust::String inspect_checkpoint(rust::Str checkpoint_path);
rust::String export_checkpoint(rust::Str checkpoint_path, rust::Str mesh_path,
                               rust::Str output_path);
rust::String prepare_mpi_run_directory(rust::Str root_path, rust::Str effective_config);
rust::String run_scalar_advection_case(
    rust::Str mesh_path, rust::Str output_path, rust::Str restart_path, bool mpi, int dimensions,
    int polynomial_order, double velocity_x, double velocity_y, double cfl, double start_time,
    double end_time, std::uint64_t max_steps, std::uint64_t history_interval,
    std::uint64_t visualization_interval, std::uint64_t checkpoint_interval,
    std::uint64_t keep_checkpoints, double initial_default, bool initial_sine,
    double initial_sine_amplitude, rust::Slice<const double> intervals,
    rust::Slice<const double> circles, rust::Slice<const rust::String> boundary_names,
    rust::Slice<const std::uint8_t> boundary_kinds, rust::Slice<const double> boundary_values);
rust::String run_euler_case(rust::Str mesh_path, rust::Str output_path, rust::Str restart_path,
                            bool limiting, double sensor_threshold, double density_floor,
                            double pressure_floor, int dimensions, int polynomial_order,
                            double heat_capacity_ratio, double cfl, double start_time,
                            double end_time, std::uint64_t max_steps,
                            std::uint64_t history_interval, std::uint64_t visualization_interval,
                            std::uint64_t checkpoint_interval, std::uint64_t keep_checkpoints,
                            rust::Slice<const double> regions,
                            rust::Slice<const rust::String> boundary_names,
                            rust::Slice<const std::uint8_t> boundary_kinds,
                            rust::Slice<const double> boundary_values);
