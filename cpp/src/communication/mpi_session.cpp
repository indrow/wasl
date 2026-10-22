#ifdef WASL_ENABLE_MPI

#include "communication/mpi_session.hpp"

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>

#include <mpi.h>

namespace wasl {

void ensure_mpi_initialized() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0)
    return;
  MPI_Init(nullptr, nullptr);
  std::atexit([] {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized == 0)
      MPI_Finalize();
  });
}

int world_rank() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int world_size() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

void world_barrier() { MPI_Barrier(MPI_COMM_WORLD); }

double world_sum(double value) {
  double result = 0.0;
  MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return result;
}

std::string broadcast_from_root(std::string value) {
  auto size = static_cast<std::uint64_t>(value.size());
  MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("MPI broadcast: string is too large");
  if (world_rank() != 0)
    value.resize(size);
  MPI_Bcast(value.data(), static_cast<int>(size), MPI_CHAR, 0, MPI_COMM_WORLD);
  return value;
}

} // namespace wasl

#endif
