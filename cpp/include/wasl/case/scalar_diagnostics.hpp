// Cell averages and total mass of scalar solver coefficients U[element][mode].

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

std::vector<double> scalar_cell_averages(const std::vector<double> &coefficients,
                                         std::size_t element_count, std::size_t mode_count);

double scalar_mass(const LineMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count);
double scalar_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients,
                   std::size_t mode_count);

} // namespace wasl
