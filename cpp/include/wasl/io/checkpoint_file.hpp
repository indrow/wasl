// Versioned HDF5 checkpoint files.

#pragma once

#include <filesystem>

#include "wasl/io/checkpoint.hpp"

namespace wasl {

void write_checkpoint(const std::filesystem::path &path, const Checkpoint &checkpoint);

// Reads format versions 1 and 2; version 1 has no partition metadata.
Checkpoint read_checkpoint(const std::filesystem::path &path);

} // namespace wasl
