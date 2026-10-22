// Run length and output cadence shared by every case type.

#pragma once

#include <cstdint>

namespace wasl {

// Intervals count steps. The final step always writes history, a snapshot, and
// a checkpoint. At most keep_checkpoints checkpoints remain on disk.
struct RunControl {
  double start_time = 0.0;
  double end_time = 0.0;
  std::uint64_t max_steps = 0;
  std::uint64_t history_interval = 0;
  std::uint64_t visualization_interval = 0;
  std::uint64_t checkpoint_interval = 0;
  std::uint64_t keep_checkpoints = 0;
};

} // namespace wasl
