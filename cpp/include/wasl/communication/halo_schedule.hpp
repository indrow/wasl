// Which ghost coefficients each rank sends and receives.
//
// One MPI rank owns one partition (ARCHITECTURE.md) and imports the one
// face-adjacent ghost layer from partition.hpp. The schedule is transport
// independent: the in-process oracle and host-staged MPI both drive it.

#pragma once

#include <cstdint>
#include <vector>

#include "wasl/mesh/partition/partition.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// Global indices of this rank's owned elements that to_rank imports.
struct HaloSend {
  int to_rank = 0;
  std::vector<std::int32_t> elements;
};

// Global indices of ghost elements this rank imports from from_rank.
struct HaloRecv {
  int from_rank = 0;
  std::vector<std::int32_t> elements;
};

// Entries are sorted by neighbor rank, then element; the schedule set over all
// ranks is symmetric and deterministic.
struct HaloSchedule {
  std::vector<HaloSend> sends;
  std::vector<HaloRecv> recvs;
};

HaloSchedule build_halo_schedule(const PreparedMesh &mesh, const Partition &partition, int rank);

// Add one element to the send or receive block of a peer, keeping elements
// sorted and unique. Used for exchanges beyond the face-adjacent layer, such as
// periodic endpoint ghosts.
void add_halo_send(HaloSchedule &schedule, int peer, std::int32_t element);
void add_halo_recv(HaloSchedule &schedule, int peer, std::int32_t element);

} // namespace wasl
