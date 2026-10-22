// MPI_COMM_WORLD helpers for case-level drivers. Only compiled with
// WASL_ENABLE_MPI; MPI types stay out of the declarations.

#pragma once

#ifdef WASL_ENABLE_MPI

#include <string>

namespace wasl {

// Initialize MPI on first use and finalize it at process exit.
void ensure_mpi_initialized();

int world_rank();
int world_size();
void world_barrier();
double world_sum(double value);

// Rank 0's value on every rank.
std::string broadcast_from_root(std::string value);

} // namespace wasl

#endif
