#include "wasl/communication/in_process_halo_exchange.hpp"

#include <algorithm>
#include <stdexcept>

namespace wasl {

void exchange_halos_in_process(
    const std::vector<HaloSchedule> &schedules, std::size_t rank_count, std::size_t mode_count,
    std::size_t variable_offset, std::size_t variable_stride,
    const std::vector<const double *> &rank_coefficients,
    const std::vector<std::map<std::int32_t, std::int32_t>> &rank_global_to_local) {
  if (schedules.size() != rank_count || rank_coefficients.size() != rank_count ||
      rank_global_to_local.size() != rank_count)
    throw std::runtime_error("halo exchange: per-rank array size mismatch");

  // The global element index pairs each send with its receive slot.
  for (std::size_t sender = 0; sender < rank_count; ++sender) {
    for (const auto &send : schedules[sender].sends) {
      const auto receiver = static_cast<std::size_t>(send.to_rank);
      const auto &receiver_slots = rank_global_to_local[receiver];
      for (const auto global : send.elements) {
        const auto source_local = rank_global_to_local[sender].at(global);
        const auto slot = receiver_slots.find(global);
        if (slot == receiver_slots.end())
          throw std::runtime_error("halo exchange: send has no matching recv");
        const auto *source = rank_coefficients[sender] + variable_offset +
                             static_cast<std::size_t>(source_local) * variable_stride;
        // The oracle API passes read-only views of arrays it is allowed to fill.
        auto *destination =
            const_cast<double *>(rank_coefficients[receiver] + variable_offset +
                                 static_cast<std::size_t>(slot->second) * variable_stride);
        std::copy_n(source, mode_count, destination);
      }
    }
  }
}

} // namespace wasl
