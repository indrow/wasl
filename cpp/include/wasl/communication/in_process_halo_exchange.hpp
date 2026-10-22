// In-process halo exchange: the deterministic oracle the MPI transport must
// reproduce, runnable in the default CPU CI without MPI.

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "wasl/communication/halo_schedule.hpp"

namespace wasl {

// Copy every rank's owned coefficients of one variable into the matching ghost
// slots of the receiving ranks. Each rank's array has its element axis laid out
// as owned elements then ghosts; an element's modes start at
// variable_offset + local * variable_stride. rank_global_to_local[r] maps every
// global element rank r touches to its local index. Throws when a send has no
// matching receive slot.
void exchange_halos_in_process(
    const std::vector<HaloSchedule> &schedules, std::size_t rank_count, std::size_t mode_count,
    std::size_t variable_offset, std::size_t variable_stride,
    const std::vector<const double *> &rank_coefficients,
    const std::vector<std::map<std::int32_t, std::int32_t>> &rank_global_to_local);

} // namespace wasl
