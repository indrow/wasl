#include "wasl/case/mpi_run_directory.hpp"

#include <stdexcept>

#ifdef WASL_ENABLE_MPI

#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

#include "communication/mpi_session.hpp"

namespace wasl {
namespace {

// The first free run-NNNNNN directory, holding effective.toml.
std::filesystem::path create_run_directory(const std::filesystem::path &root,
                                           const std::string &effective_config) {
  std::filesystem::create_directories(root);
  for (int number = 1; number <= 999999; ++number) {
    std::ostringstream name;
    name << "run-" << std::setfill('0') << std::setw(6) << number;
    const auto candidate = root / name.str();
    if (std::filesystem::create_directory(candidate)) {
      std::ofstream output(candidate / "effective.toml");
      output.exceptions(std::ios::failbit | std::ios::badbit);
      output << effective_config;
      return candidate;
    }
  }
  throw std::runtime_error("no available run directory");
}

} // namespace

std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &root,
                                                const std::string &effective_config) {
  ensure_mpi_initialized();
  if (world_size() != 2)
    throw std::runtime_error("MPI run directory: exactly two ranks are required");
  // Rank 0 creates the directory; any failure is broadcast so every rank
  // throws together instead of deadlocking.
  std::string path;
  std::string error;
  if (world_rank() == 0) {
    try {
      path = create_run_directory(root, effective_config).string();
    } catch (const std::exception &exception) {
      error = exception.what();
    }
  }
  error = broadcast_from_root(std::move(error));
  path = broadcast_from_root(std::move(path));
  if (!error.empty())
    throw std::runtime_error("MPI run directory: " + error);
  return path;
}

} // namespace wasl

#else

namespace wasl {

std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &,
                                                const std::string &) {
  throw std::runtime_error("MPI run directory: bridge was built without MPI support");
}

} // namespace wasl

#endif
