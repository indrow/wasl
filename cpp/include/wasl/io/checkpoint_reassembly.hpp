// Merge rank-local checkpoints into one global checkpoint.

#pragma once

#include <vector>

#include "wasl/io/checkpoint.hpp"

namespace wasl {

// Requires exactly one consistent part per rank that together own every
// element once.
Checkpoint reassemble_checkpoints(const std::vector<Checkpoint> &parts);

} // namespace wasl
