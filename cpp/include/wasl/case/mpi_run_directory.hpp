// Run-directory creation agreed across MPI ranks.

#pragma once

#include <filesystem>
#include <string>

namespace wasl {

// Rank 0 creates the next free run-NNNNNN directory under root and writes
// effective.toml; every rank receives the same path. Throws without MPI
// support or unless exactly two ranks run.
std::filesystem::path prepare_mpi_run_directory(const std::filesystem::path &root,
                                                const std::string &effective_config);

} // namespace wasl
