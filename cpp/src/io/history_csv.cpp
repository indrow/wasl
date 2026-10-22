#include "wasl/io/history_csv.hpp"

#include "io/text_output.hpp"

namespace wasl {

void write_csv_history(const std::filesystem::path &path, const std::vector<HistoryRow> &rows) {
  auto output = open_text_output(path);
  output << "step,time,dt,mass\n";
  for (const auto &row : rows)
    output << row.step << ',' << row.time << ',' << row.time_step << ',' << row.mass << '\n';
}

} // namespace wasl
