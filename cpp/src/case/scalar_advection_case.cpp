#include "wasl/core/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "wasl/basis/reference_element.hpp"
#include "wasl/boundaries/periodic.hpp"
#include "wasl/case/coefficient_layout.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/case/periodic_boundary.hpp"
#include "wasl/case/prepared_mesh.hpp"
#include "wasl/case/run_schedule.hpp"
#include "wasl/case/validation.hpp"
#include "wasl/dg/residual_operator.hpp"
#include "wasl/dg/triangle_residual.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

#ifdef WASL_ENABLE_MPI
#include <cstdlib>

#include <mpi.h>

#include "wasl/communication/halo_exchange.hpp"
#include "wasl/communication/host_staged_mpi.hpp"
#include "wasl/mesh/partition.hpp"
#endif

namespace wasl {
namespace {

SimulationSummary run_line_scalar_advection_case(const ScalarAdvectionCase &config,
                                                 const std::filesystem::path &prepared_path,
                                                 const PreparedMesh &prepared,
                                                 const LineMesh &mesh) {
  const auto mode_count = line_mode_count(config.polynomial_order);
  const auto element_count = mesh.num_elements();
  auto coefficients = config.restart_checkpoint.empty() ? initial_coefficients(mesh, config)
                                                        : std::vector<double>{};

  const auto mass = [&](const std::vector<double> &state) {
    return scalar_mass(mesh, state, mode_count);
  };
  const auto averages = [&](const std::vector<double> &state) {
    return scalar_cell_averages(state, element_count, mode_count);
  };
  const auto adapt_state = [&](const std::vector<double> &stored) {
    return solver_coefficients(stored, element_count, mode_count);
  };
  const auto checkpoint_state = [&](const std::vector<double> &state) {
    return checkpoint_coefficients(state, element_count, mode_count);
  };

  const ScalarAdvection equation{config.wavespeed};
  const auto exterior = make_periodic_exterior(mesh, config.polynomial_order);
  const auto cfl_step = advection_cfl_timestep(mesh, equation, config.polynomial_order, config.cfl);

  const auto initialize = [](const std::vector<double> &) {};
  const auto advance = [&](double time_step, std::vector<double> &stage,
                           SspRk3Workspace &workspace) {
    ssp_rk3_step(mesh, equation, config.polynomial_order, time_step, exterior, stage, workspace);
  };
  const auto timestep = [&](const std::vector<double> &) {
    return cfl_step > 0.0 ? cfl_step : std::numeric_limits<double>::max();
  };
  const auto synchronize = [](std::vector<double> &) {};

  return run_serial_case_loop(
      config.output_directory, prepared_path, prepared, mesh, mode_count, std::size_t{1},
      config.start_time, config.end_time, config.max_steps, config.history_interval,
      config.visualization_interval, config.checkpoint_interval, config.keep_checkpoints,
      std::move(coefficients), config.restart_checkpoint, mass, averages, adapt_state,
      checkpoint_state, initialize, advance, timestep, synchronize);
}

SimulationSummary run_triangle_scalar_advection_case(const ScalarAdvectionCase &config,
                                                     const std::filesystem::path &prepared_path,
                                                     const PreparedMesh &prepared) {
  if (config.mpi) {
    throw std::runtime_error("simulation: 2D case-level MPI is not implemented");
  }
  const auto mesh = triangle_mesh_from_prepared(prepared);
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto element_count = mesh.num_elements();
  auto coefficients = config.restart_checkpoint.empty() ? initial_coefficients(mesh, config)
                                                        : std::vector<double>{};

  const auto mass = [&](const std::vector<double> &state) {
    return scalar_mass(mesh, state, mode_count);
  };
  const auto averages = [&](const std::vector<double> &state) {
    return scalar_cell_averages(state, element_count, mode_count);
  };
  const auto adapt_state = [&](const std::vector<double> &stored) {
    return solver_coefficients(stored, element_count, mode_count);
  };
  const auto checkpoint_state = [&](const std::vector<double> &state) {
    return checkpoint_coefficients(state, element_count, mode_count);
  };

  const ScalarAdvection2D equation{config.wavespeed, config.velocity_y};
  const auto trace = periodic_triangle_trace(mesh, config.polynomial_order);
  const auto cfl_step = triangle_scalar_cfl_timestep(mesh, config.wavespeed, config.velocity_y,
                                                     config.polynomial_order, config.cfl);
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    auto exterior = periodic_triangle_exterior(trace, stage, mode_count);
    triangle_scalar_advection_residual(mesh, equation, config.polynomial_order, stage.data(),
                                       exterior.data(), rhs.data());
  };

  const auto initialize = [](const std::vector<double> &) {};
  const auto advance = [&](double time_step, std::vector<double> &stage,
                           SspRk3Workspace &workspace) {
    ssp_rk3_step(time_step, residual, stage, workspace);
  };
  const auto timestep = [&](const std::vector<double> &) {
    return cfl_step > 0.0 ? cfl_step : std::numeric_limits<double>::max();
  };
  const auto synchronize = [](std::vector<double> &) {};

  return run_serial_case_loop(
      config.output_directory, prepared_path, prepared, mesh, mode_count, std::size_t{1},
      config.start_time, config.end_time, config.max_steps, config.history_interval,
      config.visualization_interval, config.checkpoint_interval, config.keep_checkpoints,
      std::move(coefficients), config.restart_checkpoint, mass, averages, adapt_state,
      checkpoint_state, initialize, advance, timestep, synchronize);
}

#ifdef WASL_ENABLE_MPI

struct RankLineMesh {
  LineMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

void ensure_mpi() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    MPI_Init(nullptr, nullptr);
    std::atexit([] {
      int finalized = 0;
      MPI_Finalized(&finalized);
      if (finalized == 0) {
        MPI_Finalize();
      }
    });
  }
}

std::vector<std::int32_t> owned_elements(const Partition &partition, int rank) {
  std::vector<std::int32_t> owned;
  for (std::size_t element = 0; element < partition.element_owners.size(); ++element) {
    if (partition.element_owners[element] == rank) {
      owned.push_back(static_cast<std::int32_t>(element));
    }
  }
  return owned;
}

RankLineMesh rank_line_mesh(const LineMesh &global, const std::vector<std::int32_t> &owned,
                            const std::vector<std::int32_t> &ghosts) {
  RankLineMesh result;
  std::int32_t local = 0;
  for (const auto global_element : owned) {
    result.mesh.elements.push_back(global.elements.at(global_element));
    result.global_to_local[global_element] = local++;
  }
  for (const auto global_element : ghosts) {
    result.mesh.elements.push_back(global.elements.at(global_element));
    result.global_to_local[global_element] = local++;
  }
  for (const auto &face : global.faces) {
    const bool left_local = result.global_to_local.contains(face.left_element);
    const bool right_local = face.right_element == invalid_element ||
                             result.global_to_local.contains(face.right_element);
    if (!left_local || !right_local) {
      continue;
    }
    auto local_face = face;
    local_face.left_element = result.global_to_local.at(face.left_element);
    local_face.right_element = face.right_element == invalid_element
                                   ? invalid_element
                                   : result.global_to_local.at(face.right_element);
    result.mesh.faces.push_back(local_face);
  }
  return result;
}

void add_send(HaloSchedule &schedule, int peer, std::int32_t element) {
  const auto block = std::find_if(schedule.sends.begin(), schedule.sends.end(),
                                  [peer](const auto &send) { return send.to_rank == peer; });
  if (block == schedule.sends.end()) {
    schedule.sends.push_back({peer, {element}});
  } else if (std::find(block->elements.begin(), block->elements.end(), element) ==
             block->elements.end()) {
    block->elements.push_back(element);
    std::sort(block->elements.begin(), block->elements.end());
  }
}

void add_recv(HaloSchedule &schedule, int peer, std::int32_t element) {
  const auto block = std::find_if(schedule.recvs.begin(), schedule.recvs.end(),
                                  [peer](const auto &recv) { return recv.from_rank == peer; });
  if (block == schedule.recvs.end()) {
    schedule.recvs.push_back({peer, {element}});
  } else if (std::find(block->elements.begin(), block->elements.end(), element) ==
             block->elements.end()) {
    block->elements.push_back(element);
    std::sort(block->elements.begin(), block->elements.end());
  }
}

std::vector<double> local_exterior(const RankLineMesh &local, int polynomial_order,
                                   const std::vector<double> &coefficients, double left_boundary,
                                   double right_boundary) {
  const auto mode_count = line_mode_count(polynomial_order);
  std::vector<double> exterior(local.mesh.faces.size());
  for (std::size_t face = 0; face < local.mesh.faces.size(); ++face) {
    const auto right = local.mesh.faces[face].right_element;
    if (right == invalid_element) {
      exterior[face] = local.mesh.faces[face].outward_normal < 0.0 ? left_boundary : right_boundary;
    } else {
      exterior[face] =
          evaluate_state(coefficients.data() + right * mode_count, polynomial_order, -1.0);
    }
  }
  return exterior;
}

std::filesystem::path rank_checkpoint_path(const std::filesystem::path &path, int rank) {
  auto name = path.filename().string();
  const auto marker = name.rfind("-rank-");
  const auto extension = name.rfind('.');
  if (marker == std::string::npos || extension == std::string::npos || marker >= extension) {
    if (rank == 0) {
      return path;
    }
    throw std::runtime_error("simulation restart: MPI checkpoint must have a -rank-000000 suffix");
  }
  std::ostringstream rank_name;
  rank_name << std::setfill('0') << std::setw(6) << rank;
  name.replace(marker + 6, extension - marker - 6, rank_name.str());
  return path.parent_path() / name;
}

SimulationSummary run_distributed_scalar_advection_case(const ScalarAdvectionCase &config,
                                                        const std::filesystem::path &prepared_path,
                                                        PreparedMesh prepared,
                                                        const LineMesh &global_mesh) {
  ensure_mpi();
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
  if (rank_count != 2) {
    throw std::runtime_error("simulation: MPI scalar-advection driver requires exactly two ranks");
  }

  const auto partition = partition_mesh(prepared, rank_count);
  prepared.partition_rank_count = partition.rank_count;
  prepared.element_owners = partition.element_owners;
  prepared.face_owners = partition.face_owners;
  const auto owned = owned_elements(partition, rank);
  auto ghosts = rank_ghost_elements(prepared, partition, rank);
  const auto last = static_cast<std::int32_t>(global_mesh.num_elements() - 1);
  const auto periodic_ghost = rank == 0 ? last : 0;
  if (std::find(ghosts.begin(), ghosts.end(), periodic_ghost) == ghosts.end()) {
    ghosts.push_back(periodic_ghost);
    std::sort(ghosts.begin(), ghosts.end());
  }
  const auto local = rank_line_mesh(global_mesh, owned, ghosts);
  const auto mode_count = line_mode_count(config.polynomial_order);
  const auto global_element_count = global_mesh.num_elements();
  std::vector<double> coefficients(local.mesh.num_elements() * mode_count, 0.0);

  std::uint64_t step = 0;
  double time = config.start_time;
  if (config.restart_checkpoint.empty()) {
    const auto initial = initial_coefficients(global_mesh, config);
    for (const auto global : owned) {
      const auto local_element = local.global_to_local.at(global);
      std::copy_n(initial.data() + global * mode_count, mode_count,
                  coefficients.data() + local_element * mode_count);
    }
  } else {
    const auto checkpoint = read_checkpoint(rank_checkpoint_path(config.restart_checkpoint, rank));
    const std::vector<std::int64_t> expected(owned.begin(), owned.end());
    if (checkpoint.variable_count != 1 || checkpoint.mode_count != mode_count ||
        checkpoint.partition_rank != rank || checkpoint.partition_rank_count != rank_count ||
        checkpoint.global_element_count != global_element_count ||
        checkpoint.element_indices != expected ||
        std::filesystem::weakly_canonical(checkpoint.prepared_mesh) != prepared_path) {
      throw std::runtime_error("simulation restart: checkpoint layout or mesh provenance mismatch");
    }
    step = checkpoint.step;
    time = checkpoint.time;
    const auto restored = solver_coefficients(checkpoint);
    for (std::size_t local_element = 0; local_element < owned.size(); ++local_element) {
      std::copy_n(restored.data() + local_element * mode_count, mode_count,
                  coefficients.data() + local_element * mode_count);
    }
  }
  if (!std::isfinite(time) || time < config.start_time || time > config.end_time ||
      step > config.max_steps) {
    throw std::runtime_error("simulation restart: checkpoint is outside the configured run range");
  }

  std::filesystem::create_directories(config.output_directory);
  const auto owned_mesh = owned_mesh_slice(prepared, rank);
  std::vector<HistoryRow> history;
  std::vector<std::pair<double, std::string>> snapshots;
  std::deque<std::filesystem::path> checkpoints;
  SimulationSummary summary{step, time, {}};

  auto global_mass = [&] {
    double local_mass = 0.0;
    for (const auto global : owned) {
      const auto local_element = local.global_to_local.at(global);
      local_mass += std::sqrt(2.0) * element_map(global_mesh, global).jacobian() *
                    coefficients[local_element * mode_count];
    }
    double result = 0.0;
    MPI_Allreduce(&local_mass, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
  };
  auto write_visualization = [&] {
    const auto base = numbered_step_name("solution", step, "");
    std::ostringstream rank_name;
    rank_name << base << "-rank-" << std::setfill('0') << std::setw(6) << rank << ".vtu";
    std::vector<double> values;
    values.reserve(owned.size());
    for (const auto global : owned) {
      values.push_back(coefficients[local.global_to_local.at(global) * mode_count] /
                       std::sqrt(2.0));
    }
    write_vtu(config.output_directory / rank_name.str(), owned_mesh, values);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      std::vector<std::string> pieces;
      for (int source = 0; source < rank_count; ++source) {
        std::ostringstream piece;
        piece << base << "-rank-" << std::setfill('0') << std::setw(6) << source << ".vtu";
        pieces.push_back(piece.str());
      }
      const auto manifest = base + ".pvtu";
      write_pvtu(config.output_directory / manifest, pieces);
      snapshots.emplace_back(time, manifest);
    }
  };
  auto write_restart = [&] {
    const auto base = numbered_step_name("checkpoint", step, "");
    std::ostringstream name;
    name << base << "-rank-" << std::setfill('0') << std::setw(6) << rank << ".h5";
    const auto path = config.output_directory / name.str();
    Checkpoint checkpoint;
    checkpoint.time = time;
    checkpoint.step = step;
    checkpoint.variable_count = 1;
    checkpoint.mode_count = mode_count;
    checkpoint.element_count = owned.size();
    checkpoint.prepared_mesh = prepared_path.string();
    checkpoint.partition_rank = rank;
    checkpoint.partition_rank_count = rank_count;
    checkpoint.global_element_count = global_element_count;
    checkpoint.element_indices.assign(owned.begin(), owned.end());
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      for (const auto global : owned) {
        checkpoint.coefficients.push_back(
            coefficients[local.global_to_local.at(global) * mode_count + mode]);
      }
    }
    write_checkpoint(path, checkpoint);
    checkpoints.push_back(path);
    if (checkpoints.size() > config.keep_checkpoints) {
      std::filesystem::remove(checkpoints.front());
      checkpoints.pop_front();
    }
    summary.checkpoint = path;
  };

  history.push_back({step, time, 0.0, global_mass()});
  write_visualization();
  auto schedule = build_halo_schedule(prepared, partition, rank);
  const int peer = rank == 0 ? rank_count - 1 : 0;
  add_send(schedule, peer, rank == 0 ? 0 : last);
  add_recv(schedule, peer, periodic_ghost);
  HostStagedHaloExchange exchange(MPI_COMM_WORLD, schedule, mode_count, 0, mode_count,
                                  local.global_to_local);
  const ScalarAdvection equation{config.wavespeed};
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    exchange.start(stage.data());
    exchange.finish(stage.data());
    const auto left = evaluate_state(stage.data() + local.global_to_local.at(last) * mode_count,
                                     config.polynomial_order, 1.0);
    const auto right = evaluate_state(stage.data() + local.global_to_local.at(0) * mode_count,
                                      config.polynomial_order, -1.0);
    const auto exterior = local_exterior(local, config.polynomial_order, stage, left, right);
    scalar_advection_residual(local.mesh, equation, config.polynomial_order, stage.data(),
                              exterior.data(), rhs.data());
  };
  const auto cfl_step =
      advection_cfl_timestep(global_mesh, equation, config.polynomial_order, config.cfl);
  SspRk3Workspace workspace;
  while (time < config.end_time && step < config.max_steps) {
    const auto time_step =
        std::min(cfl_step > 0.0 ? cfl_step : config.end_time - time, config.end_time - time);
    if (time_step <= 0.0 || time + time_step == time) {
      throw std::runtime_error("simulation: timestep cannot advance time");
    }
    ssp_rk3_step(time_step, residual, coefficients, workspace);
    time += time_step;
    ++step;
    const bool final = time >= config.end_time || step == config.max_steps;
    if (step % config.history_interval == 0 || final) {
      const auto current_mass = global_mass();
      if (rank == 0) {
        history.push_back({step, time, time_step, current_mass});
      }
    }
    if (step % config.visualization_interval == 0 || final) {
      write_visualization();
    }
    if (step % config.checkpoint_interval == 0 || final) {
      write_restart();
    }
  }
  if (time < config.end_time) {
    throw std::runtime_error("simulation: max_steps reached before end_time");
  }
  if (summary.checkpoint.empty()) {
    write_restart();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    write_csv_history(config.output_directory / "history.csv", history);
    write_pvd(config.output_directory / "solution.pvd", snapshots);
  }
  summary.step = step;
  summary.time = time;
  return summary;
}

std::string broadcast_string(std::string value, int rank) {
  auto size = static_cast<std::uint64_t>(value.size());
  MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("MPI broadcast: string is too large");
  }
  if (rank != 0) {
    value.resize(size);
  }
  MPI_Bcast(value.data(), static_cast<int>(size), MPI_CHAR, 0, MPI_COMM_WORLD);
  return value;
}

} // namespace

std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &root,
                                                const std::string &effective_config) {
  ensure_mpi();
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
  if (rank_count != 2) {
    throw std::runtime_error("MPI run directory: exactly two ranks are required");
  }
  std::string path;
  std::string error;
  if (rank == 0) {
    try {
      std::filesystem::create_directories(root);
      for (int number = 1; number <= 999999; ++number) {
        std::ostringstream name;
        name << "run-" << std::setfill('0') << std::setw(6) << number;
        const auto candidate = root / name.str();
        if (std::filesystem::create_directory(candidate)) {
          std::ofstream output(candidate / "effective.toml");
          output.exceptions(std::ios::failbit | std::ios::badbit);
          output << effective_config;
          path = candidate.string();
          break;
        }
      }
      if (path.empty()) {
        throw std::runtime_error("no available run directory");
      }
    } catch (const std::exception &exception) {
      error = exception.what();
    }
  }
  error = broadcast_string(std::move(error), rank);
  path = broadcast_string(std::move(path), rank);
  if (!error.empty()) {
    throw std::runtime_error("MPI run directory: " + error);
  }
  return path;
}

#else

} // namespace

std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &root,
                                                const std::string &effective_config) {
  (void)root;
  (void)effective_config;
  throw std::runtime_error("MPI run directory: bridge was built without MPI support");
}

#endif

SimulationSummary run_scalar_advection_case(const ScalarAdvectionCase &config) {
  validate_scalar_advection_case(config);
  const auto prepared_path = std::filesystem::weakly_canonical(config.prepared_mesh);
  const auto prepared = read_prepared_mesh(prepared_path);
  if (config.dimensions != prepared.dimension) {
    throw std::runtime_error("simulation: case dimension does not match prepared mesh");
  }
  if (config.dimensions == 2) {
    return run_triangle_scalar_advection_case(config, prepared_path, prepared);
  }
  const auto mesh = line_mesh_from_prepared(prepared);
#ifdef WASL_ENABLE_MPI
  if (config.mpi) {
    return run_distributed_scalar_advection_case(config, prepared_path, prepared, mesh);
  }
#else
  if (config.mpi) {
    throw std::runtime_error("simulation: bridge was built without MPI support");
  }
#endif
  return run_line_scalar_advection_case(config, prepared_path, prepared, mesh);
}

} // namespace wasl
