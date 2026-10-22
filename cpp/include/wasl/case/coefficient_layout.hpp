#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {

// Convert between solver layout U[element][mode] and checkpoint layout
// U[variable][mode][element] for the scalar (one-variable) path.
std::vector<double> checkpoint_coefficients(const std::vector<double> &coefficients,
                                            std::size_t element_count, std::size_t mode_count);
std::vector<double> solver_coefficients(const std::vector<double> &coefficients,
                                        std::size_t element_count, std::size_t mode_count);
std::vector<double> solver_coefficients(const Checkpoint &checkpoint);

std::vector<double> scalar_cell_averages(const std::vector<double> &coefficients,
                                         std::size_t element_count, std::size_t mode_count);

double scalar_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count);
double scalar_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count);

EulerConserved average_euler_state(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count,
                                   std::size_t element);
std::vector<double> euler_averages(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count);
double euler_density_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                          std::size_t mode_count);
double euler_density_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                          std::size_t mode_count);

} // namespace wasl
