// Host-staged nonblocking MPI halo exchange over a prepared-mesh partition.
//
// The production transport for halo_schedule.hpp (ARCHITECTURE.md). Each stage
// posts nonblocking receives, packs and posts nonblocking sends of owned
// elements that peers import, lets the caller do local work while the exchange
// is in flight, and finally unpacks the received ghost coefficients.
//
// The header is MPI-only and is compiled solely into targets built with
// WASL_ENABLE_MPI; in_process_halo_exchange.hpp remains the CPU-CI oracle.
//
// Coefficient layout is U[variable][mode][element] with the element axis laid
// out as owned elements then ghost elements, matching the oracle. The residual
// reads a neighbor's modes directly for a partition-boundary face, so this
// exchange only has to fill ghost coefficients; no residual change is needed.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include <mpi.h>

#include "wasl/communication/halo_schedule.hpp"

namespace wasl {

// One round of host-staged nonblocking halo exchange for a single variable.
// Posts Irecv first, then Isend. The caller may do local work between start()
// and finish(), which waits and unpacks received buffers into ghost slots.
//
// ranks is the MPI communicator (typically MPI_COMM_WORLD). tag disambiguates
// concurrent exchanges (one variable at a time here, so a fixed tag is fine).
class HostStagedHaloExchange {
public:
  HostStagedHaloExchange(MPI_Comm comm, const HaloSchedule &schedule, std::size_t mode_count,
                         std::size_t variable_offset, std::size_t variable_stride,
                         const std::map<std::int32_t, std::int32_t> &global_to_local)
      : comm_(comm), schedule_(schedule), mode_count_(mode_count),
        variable_offset_(variable_offset), variable_stride_(variable_stride),
        global_to_local_(global_to_local) {
    int flag = 0;
    MPI_Initialized(&flag);
    mpi_ready_ = flag != 0;
    // Size persistent send/recv buffers from the schedule.
    for (const auto &send : schedule_.sends) {
      send_buffers_.push_back(std::vector<double>(send.elements.size() * mode_count_));
    }
    for (const auto &recv : schedule_.recvs) {
      recv_buffers_.push_back(std::vector<double>(recv.elements.size() * mode_count_));
    }
  }

  // Start the exchange: post receives, pack owned coefficients, post sends.
  // coefficients is this rank's full U[variable][mode][element] array.
  void start(const double *coefficients) {
    if (!mpi_ready_)
      return;
    requests_.clear();
    // Post receives first to guarantee send buffers have a destination.
    for (std::size_t i = 0; i < schedule_.recvs.size(); ++i) {
      const auto &recv = schedule_.recvs[i];
      auto &buf = recv_buffers_[i];
      requests_.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(buf.data(), static_cast<int>(buf.size()), MPI_DOUBLE, recv.from_rank, kTag, comm_,
                &requests_.back());
    }
    // Pack owned element modes into send buffers and post sends.
    for (std::size_t i = 0; i < schedule_.sends.size(); ++i) {
      const auto &send = schedule_.sends[i];
      auto &buf = send_buffers_[i];
      for (std::size_t k = 0; k < send.elements.size(); ++k) {
        // send.elements[k] is a global index; translate to this rank's local
        // (owned) layout before reading the source coefficients.
        const auto local = global_to_local_.at(send.elements[k]);
        const auto *src =
            coefficients + variable_offset_ + static_cast<std::size_t>(local) * variable_stride_;
        std::copy_n(src, mode_count_, buf.data() + k * mode_count_);
      }
      requests_.push_back(MPI_REQUEST_NULL);
      MPI_Isend(buf.data(), static_cast<int>(buf.size()), MPI_DOUBLE, send.to_rank, kTag, comm_,
                &requests_.back());
    }
  }

  // Wait for all sends/receives, then scatter received buffers into ghost slots.
  void finish(double *coefficients) {
    if (!mpi_ready_ || requests_.empty())
      return;
    MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE);
    // Unpack each recv buffer into the ghost region of the coefficient array.
    for (std::size_t i = 0; i < schedule_.recvs.size(); ++i) {
      const auto &recv = schedule_.recvs[i];
      const auto &buf = recv_buffers_[i];
      for (std::size_t k = 0; k < recv.elements.size(); ++k) {
        // A missing ghost slot means the schedule and local layout disagree;
        // fail loudly instead of leaving stale ghost coefficients.
        const auto local = global_to_local_.at(recv.elements[k]);
        auto *dst =
            coefficients + variable_offset_ + static_cast<std::size_t>(local) * variable_stride_;
        std::copy_n(buf.data() + k * mode_count_, mode_count_, dst);
      }
    }
  }

private:
  static constexpr int kTag = 42;
  MPI_Comm comm_;
  const HaloSchedule &schedule_;
  std::size_t mode_count_;
  std::size_t variable_offset_;
  std::size_t variable_stride_;
  const std::map<std::int32_t, std::int32_t> &global_to_local_;
  std::vector<std::vector<double>> send_buffers_;
  std::vector<std::vector<double>> recv_buffers_;
  std::vector<MPI_Request> requests_;
  bool mpi_ready_ = false;
};

} // namespace wasl
