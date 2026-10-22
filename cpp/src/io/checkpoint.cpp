#include "wasl/io/checkpoint.hpp"

#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace wasl {
namespace {

// variable_count * mode_count * element_count == coefficients.size(), without
// overflow.
bool coefficient_layout_matches(const Checkpoint &checkpoint) {
  if (checkpoint.variable_count == 0 || checkpoint.mode_count == 0 ||
      checkpoint.element_count == 0 || checkpoint.prepared_mesh.empty() ||
      checkpoint.variable_count > std::numeric_limits<std::uint64_t>::max() / checkpoint.mode_count)
    return false;
  const auto variable_modes = checkpoint.variable_count * checkpoint.mode_count;
  return variable_modes <= std::numeric_limits<std::uint64_t>::max() / checkpoint.element_count &&
         variable_modes * checkpoint.element_count == checkpoint.coefficients.size();
}

} // namespace

std::uint64_t checkpoint_global_element_count(const Checkpoint &checkpoint) {
  return checkpoint.global_element_count == 0 ? checkpoint.element_count
                                              : checkpoint.global_element_count;
}

std::vector<std::int64_t> checkpoint_element_indices(const Checkpoint &checkpoint) {
  if (!checkpoint.element_indices.empty())
    return checkpoint.element_indices;
  std::vector<std::int64_t> indices(checkpoint.element_count);
  std::iota(indices.begin(), indices.end(), 0);
  return indices;
}

void validate_checkpoint(const Checkpoint &checkpoint) {
  if (!coefficient_layout_matches(checkpoint))
    throw std::runtime_error("checkpoint: coefficient size does not match declared layout");
  if (checkpoint.partition_rank_count < 1 || checkpoint.partition_rank < 0 ||
      checkpoint.partition_rank >= checkpoint.partition_rank_count)
    throw std::runtime_error("checkpoint: invalid partition rank");
  const auto global_count = checkpoint_global_element_count(checkpoint);
  const auto indices = checkpoint_element_indices(checkpoint);
  if (indices.size() != checkpoint.element_count || global_count < checkpoint.element_count)
    throw std::runtime_error("checkpoint: invalid partition element count");
  std::vector<unsigned char> seen(global_count, 0);
  for (const auto index : indices) {
    if (index < 0 || static_cast<std::uint64_t>(index) >= global_count ||
        seen[static_cast<std::size_t>(index)] != 0)
      throw std::runtime_error("checkpoint: invalid partition element index");
    seen[static_cast<std::size_t>(index)] = 1;
  }
}

} // namespace wasl
