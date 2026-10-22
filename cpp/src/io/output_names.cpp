#include "wasl/io/output_names.hpp"

#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace wasl {
namespace {

std::string zero_padded(std::uint64_t value) {
  std::ostringstream text;
  text << std::setfill('0') << std::setw(6) << value;
  return text.str();
}

constexpr const char *rank_marker = "-rank-";

} // namespace

std::string numbered_step_name(const std::string &prefix, std::uint64_t step,
                               const std::string &extension) {
  return prefix + '-' + zero_padded(step) + extension;
}

std::string rank_piece_name(const std::string &base, int rank, const std::string &extension) {
  return base + rank_marker + zero_padded(static_cast<std::uint64_t>(rank)) + extension;
}

std::filesystem::path rank_checkpoint_path(const std::filesystem::path &path, int rank) {
  auto name = path.filename().string();
  const auto marker = name.rfind(rank_marker);
  const auto extension = name.rfind('.');
  if (marker == std::string::npos || extension == std::string::npos || marker >= extension) {
    if (rank == 0)
      return path;
    throw std::runtime_error("simulation restart: MPI checkpoint must have a -rank-000000 suffix");
  }
  const auto digits = marker + std::char_traits<char>::length(rank_marker);
  name.replace(digits, extension - digits, zero_padded(static_cast<std::uint64_t>(rank)));
  return path.parent_path() / name;
}

} // namespace wasl
