// Host-staged halo exchange schedule and coefficient marshalling.
//
// One MPI rank owns one partition of elements (ARCHITECTURE.md). Each rank imports the
// one face-adjacent ghost-element layer computed by partition.hpp. This header
// describes what to send and receive and how to pack and unpack modal
// coefficients, independent of the transport. The transport is host-staged
// nonblocking MPI in production (ARCHITECTURE.md); the in-process exchange below is the
// deterministic oracle the regression drives without an MPI toolchain.
//
// Coefficient layout is U[variable][mode][element] (NUMERICAL_CONVENTIONS.md).
// On a rank, the element axis is owned elements followed by ghost elements, so
// the existing residual can read a neighbor's modes directly for a face whose
// other side is an imported ghost: only the coefficient array and the mesh
// indices change, not the residual.

#pragma once

#include <cstdint>
#include <map>
#include <vector>

namespace wasl {

struct PreparedMesh;
struct Partition;

// What this rank must send to and receive from each neighbor rank.
//
// sends[r].elements are global owned element indices that rank r needs as a
// ghost. recvs[r].elements are global ghost element indices this rank imports
// from rank r. Both are sorted; the schedule is symmetric and deterministic.
struct HaloSend {
  int to_rank = 0;
  std::vector<std::int32_t> elements; // global owned element indices to send
};
struct HaloRecv {
  int from_rank = 0;
  std::vector<std::int32_t> elements; // global ghost element indices received
};
struct HaloSchedule {
  std::vector<HaloSend> sends;
  std::vector<HaloRecv> recvs;
};

// Derive the send/recv schedule for one rank from ownership and the ghost layer.
// A rank sends its owned elements that some other rank imports, and receives the
// ghost elements it imports. Entries are sorted by neighbor rank then element.
HaloSchedule build_halo_schedule(const PreparedMesh &mesh, const Partition &partition, int rank);

// In-process oracle exchange for one variable: copy every rank's owned
// coefficients into the matching ghost slots of the receiving ranks. This is the
// semantics the MPI transport must reproduce; it runs without MPI so the
// regression proves halo correctness in the default CPU CI.
//
// Each rank's coefficient array is U[variable][mode][element] with the element
// axis laid out as owned elements then ghost elements.
// rank_global_to_local[r] maps a global element index to its local index in
// rank r's array, for every element the rank touches (owned or ghost). The
// schedule set must be consistent: every send element has a matching entry in
// the destination rank's map.
void exchange_halos_in_process(
    const std::vector<HaloSchedule> &schedules, std::size_t rank_count, std::size_t mode_count,
    std::size_t variable_offset, std::size_t variable_stride,
    const std::vector<const double *> &rank_coefficients,
    const std::vector<std::map<std::int32_t, std::int32_t>> &rank_global_to_local);

} // namespace wasl
