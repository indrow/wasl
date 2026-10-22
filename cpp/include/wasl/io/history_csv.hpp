// Per-step run history as CSV.

#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

namespace wasl {

struct HistoryRow {
  std::uint64_t step;
  double time;
  double time_step;
  double mass;
};

void write_csv_history(const std::filesystem::path &path, const std::vector<HistoryRow> &rows);

} // namespace wasl
