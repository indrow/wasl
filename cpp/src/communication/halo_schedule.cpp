#include "wasl/communication/halo_schedule.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>

namespace wasl {
namespace {

// A rank's ghosts grouped by owning rank, each group ascending.
std::map<int, std::vector<std::int32_t>> ghosts_by_owner(const PreparedMesh &mesh,
                                                         const Partition &partition, int rank) {
  std::map<int, std::vector<std::int32_t>> grouped;
  for (const auto global : rank_ghost_elements(mesh, partition, rank))
    grouped[partition.element_owners.at(static_cast<std::size_t>(global))].push_back(global);
  return grouped;
}

// Insert into the block for peer, creating it on first use.
template <typename Block, typename Peer>
void add_to_block(std::vector<Block> &blocks, int peer, std::int32_t element, Peer peer_of) {
  const auto block = std::find_if(blocks.begin(), blocks.end(),
                                  [&](const Block &value) { return peer_of(value) == peer; });
  if (block == blocks.end()) {
    blocks.push_back({peer, {element}});
  } else if (std::find(block->elements.begin(), block->elements.end(), element) ==
             block->elements.end()) {
    block->elements.push_back(element);
    std::sort(block->elements.begin(), block->elements.end());
  }
}

} // namespace

HaloSchedule build_halo_schedule(const PreparedMesh &mesh, const Partition &partition, int rank) {
  if (rank < 0 || rank >= partition.rank_count)
    throw std::runtime_error("halo schedule: rank is out of range");

  HaloSchedule schedule;
  for (auto &[from_rank, elements] : ghosts_by_owner(mesh, partition, rank))
    schedule.recvs.push_back(HaloRecv{from_rank, std::move(elements)});

  // Sends invert the receive relation: this rank sends each owned element that
  // appears in another rank's ghost set.
  for (int peer = 0; peer < partition.rank_count; ++peer) {
    if (peer == rank)
      continue;
    std::vector<std::int32_t> to_send;
    for (const auto global : rank_ghost_elements(mesh, partition, peer)) {
      if (partition.element_owners.at(static_cast<std::size_t>(global)) == rank)
        to_send.push_back(global);
    }
    if (!to_send.empty()) {
      std::sort(to_send.begin(), to_send.end());
      schedule.sends.push_back(HaloSend{peer, std::move(to_send)});
    }
  }
  return schedule;
}

void add_halo_send(HaloSchedule &schedule, int peer, std::int32_t element) {
  add_to_block(schedule.sends, peer, element, [](const HaloSend &send) { return send.to_rank; });
}

void add_halo_recv(HaloSchedule &schedule, int peer, std::int32_t element) {
  add_to_block(schedule.recvs, peer, element, [](const HaloRecv &recv) { return recv.from_rank; });
}

} // namespace wasl
