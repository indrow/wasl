#include "wasl/io/checkpoint_reassembly.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace wasl {
namespace {

bool consistent_with(const Checkpoint &part, const Checkpoint &first) {
  return part.time == first.time && part.step == first.step &&
         part.variable_count == first.variable_count && part.mode_count == first.mode_count &&
         checkpoint_global_element_count(part) == checkpoint_global_element_count(first) &&
         part.partition_rank_count == first.partition_rank_count &&
         part.prepared_mesh == first.prepared_mesh;
}

} // namespace

Checkpoint reassemble_checkpoints(const std::vector<Checkpoint> &parts) {
  if (parts.empty())
    throw std::runtime_error("checkpoint reassembly: no rank-local checkpoints");
  const auto &first = parts.front();
  validate_checkpoint(first);
  const auto global_count = checkpoint_global_element_count(first);
  if (parts.size() != static_cast<std::size_t>(first.partition_rank_count))
    throw std::runtime_error("checkpoint reassembly: rank count mismatch");

  Checkpoint result;
  result.time = first.time;
  result.step = first.step;
  result.variable_count = first.variable_count;
  result.mode_count = first.mode_count;
  result.element_count = global_count;
  result.prepared_mesh = first.prepared_mesh;
  result.coefficients.resize(first.variable_count * first.mode_count * global_count);

  std::vector<unsigned char> seen_rank(static_cast<std::size_t>(first.partition_rank_count), 0);
  std::vector<unsigned char> seen_element(global_count, 0);
  for (const auto &part : parts) {
    validate_checkpoint(part);
    auto &rank_seen = seen_rank[static_cast<std::size_t>(part.partition_rank)];
    if (!consistent_with(part, first) || rank_seen != 0)
      throw std::runtime_error("checkpoint reassembly: inconsistent rank-local checkpoints");
    rank_seen = 1;

    const auto indices = checkpoint_element_indices(part);
    for (std::size_t local = 0; local < indices.size(); ++local) {
      const auto global = static_cast<std::size_t>(indices[local]);
      if (seen_element[global] != 0)
        throw std::runtime_error("checkpoint reassembly: duplicate owned element");
      seen_element[global] = 1;
      for (std::uint64_t variable = 0; variable < part.variable_count; ++variable) {
        for (std::uint64_t mode = 0; mode < part.mode_count; ++mode) {
          const auto source = (variable * part.mode_count + mode) * part.element_count + local;
          const auto target = (variable * part.mode_count + mode) * global_count + global;
          result.coefficients[target] = part.coefficients[source];
        }
      }
    }
  }
  if (std::find(seen_element.begin(), seen_element.end(), 0) != seen_element.end())
    throw std::runtime_error("checkpoint reassembly: missing owned element");
  return result;
}

} // namespace wasl
