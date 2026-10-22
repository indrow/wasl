// A restart checkpoint: modal coefficients U[variable][mode][element] with the
// time, step, layout, and prepared-mesh provenance needed to resume.
//
// A rank-local checkpoint stores only the elements its rank owns, listed by
// global index in element_indices. A serial checkpoint may leave
// element_indices empty and global_element_count zero; both then default to
// the identity over element_count.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace wasl {

struct Checkpoint {
  double time = 0.0;
  std::uint64_t step = 0;
  std::uint64_t variable_count = 0;
  std::uint64_t mode_count = 0;
  std::uint64_t element_count = 0;
  std::string prepared_mesh;
  std::vector<double> coefficients;
  int partition_rank = 0;
  int partition_rank_count = 1;
  std::uint64_t global_element_count = 0;
  std::vector<std::int64_t> element_indices;
};

std::uint64_t checkpoint_global_element_count(const Checkpoint &checkpoint);
std::vector<std::int64_t> checkpoint_element_indices(const Checkpoint &checkpoint);

// Throws std::runtime_error when the coefficient size, partition rank, or
// owned element indices are inconsistent.
void validate_checkpoint(const Checkpoint &checkpoint);

} // namespace wasl
