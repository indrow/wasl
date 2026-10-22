// ParaView collection of time-stamped snapshot files.

#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace wasl {

// snapshots holds (time, file name relative to the collection).
void write_pvd(const std::filesystem::path &path,
               const std::vector<std::pair<double, std::string>> &snapshots);

} // namespace wasl
