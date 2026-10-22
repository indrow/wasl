// Conversion between the scalar solver layout U[element][mode] and the
// persisted layout U[variable][mode][element] with one variable.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/io/checkpoint.hpp"

namespace wasl {

std::vector<double> checkpoint_coefficients(const std::vector<double> &coefficients,
                                            std::size_t element_count, std::size_t mode_count);
std::vector<double> solver_coefficients(const std::vector<double> &coefficients,
                                        std::size_t element_count, std::size_t mode_count);
std::vector<double> solver_coefficients(const Checkpoint &checkpoint);

} // namespace wasl
