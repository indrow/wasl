// Outcome of a completed case run.

#pragma once

#include <cstdint>
#include <filesystem>

namespace wasl {

struct SimulationSummary {
  std::uint64_t step;
  double time;
  std::filesystem::path checkpoint; // the last checkpoint written
};

} // namespace wasl
