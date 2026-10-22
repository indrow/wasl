// Parallel VTU manifest referencing rank-local pieces.

#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace wasl {

void write_pvtu(const std::filesystem::path &path, const std::vector<std::string> &pieces,
                std::size_t component_count = 1);

} // namespace wasl
