// Deterministic output file names inside a run directory.

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace wasl {

// "<prefix>-<step:06>" followed by extension, for example solution-000010.vtu.
std::string numbered_step_name(const std::string &prefix, std::uint64_t step,
                               const std::string &extension);

// "<base>-rank-<rank:06>" followed by extension.
std::string rank_piece_name(const std::string &base, int rank, const std::string &extension);

// The sibling of a "...-rank-NNNNNN.ext" checkpoint for another rank. A path
// without a rank suffix is accepted for rank 0 only.
std::filesystem::path rank_checkpoint_path(const std::filesystem::path &path, int rank);

} // namespace wasl
