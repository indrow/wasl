#include "wasl/communication/halo_exchange.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp" // invalid_element
#include "wasl/mesh/partition.hpp"

namespace wasl {

namespace {

// Collect, for one rank, the ghost elements grouped by their owning rank. The
// result is sorted by owning rank and by element within each group (rank_ghost
// _elements returns a sorted list, so each group preserves that order).
std::map<int, std::vector<std::int32_t>> ghosts_by_owner(const PreparedMesh &mesh,
                                                         const Partition &partition, int rank) {
  const auto ghosts = rank_ghost_elements(mesh, partition, rank);
  std::map<int, std::vector<std::int32_t>> grouped;
  for (const auto global : ghosts) {
    const auto owner = partition.element_owners.at(static_cast<std::size_t>(global));
    grouped[owner].push_back(global);
  }
  return grouped;
}

} // namespace

HaloSchedule build_halo_schedule(const PreparedMesh &mesh, const Partition &partition, int rank) {
  if (rank < 0 || rank >= partition.rank_count) {
    throw std::runtime_error("halo schedule: rank is out of range");
  }

  HaloSchedule schedule;

  // Recvs: this rank's ghosts, grouped by the rank that owns them.
  const auto recv_groups = ghosts_by_owner(mesh, partition, rank);
  for (const auto &[from_rank, elements] : recv_groups) {
    schedule.recvs.push_back(HaloRecv{from_rank, elements});
  }

  // Sends: invert the recv relationship. For every other rank q, find which of
  // this rank's owned elements appear in q's ghost set, i.e. are face-adjacent
  // to an element q owns. Compute q's ghosts and keep those owned by this rank.
  for (int q = 0; q < partition.rank_count; ++q) {
    if (q == rank) {
      continue;
    }
    const auto q_ghosts = rank_ghost_elements(mesh, partition, q);
    std::vector<std::int32_t> to_send;
    for (const auto global : q_ghosts) {
      if (partition.element_owners.at(static_cast<std::size_t>(global)) == rank) {
        to_send.push_back(global);
      }
    }
    if (!to_send.empty()) {
      std::sort(to_send.begin(), to_send.end());
      schedule.sends.push_back(HaloSend{q, std::move(to_send)});
    }
  }
  return schedule;
}

void exchange_halos_in_process(
    const std::vector<HaloSchedule> &schedules, std::size_t rank_count, std::size_t mode_count,
    std::size_t variable_offset, std::size_t variable_stride,
    const std::vector<const double *> &rank_coefficients,
    const std::vector<std::map<std::int32_t, std::int32_t>> &rank_global_to_local) {
  if (schedules.size() != rank_count || rank_coefficients.size() != rank_count ||
      rank_global_to_local.size() != rank_count) {
    throw std::runtime_error("halo exchange: per-rank array size mismatch");
  }
  // For each sender rank, push its owned modes into the matching ghost slot of
  // the destination rank. The global element index pairs the send and recv
  // sides through each rank's global-to-local map.
  for (std::size_t sender = 0; sender < rank_count; ++sender) {
    for (const auto &send : schedules[sender].sends) {
      const auto to = static_cast<std::size_t>(send.to_rank);
      const auto &dst_map = rank_global_to_local[to];
      for (const auto global : send.elements) {
        const auto src_local = rank_global_to_local[sender].at(global);
        const auto found = dst_map.find(global);
        if (found == dst_map.end()) {
          throw std::runtime_error("halo exchange: send has no matching recv");
        }
        const auto *src = rank_coefficients[sender] + variable_offset +
                          static_cast<std::size_t>(src_local) * variable_stride;
        auto *dst = const_cast<double *>(rank_coefficients[to] + variable_offset +
                                         static_cast<std::size_t>(found->second) * variable_stride);
        std::copy_n(src, mode_count, dst);
      }
    }
  }
}

} // namespace wasl
