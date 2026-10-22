// Cell averages and total density of Euler coefficients U[variable][mode][element].

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// Packed [element][variable] for visualization.
std::vector<double> euler_averages(const std::vector<double> &coefficients,
                                   std::size_t element_count, std::size_t mode_count);

// Density mass reads mode zero of variable zero, which is independent of
// mode_count in this layout.
double euler_density_mass(const LineMesh &mesh, const std::vector<double> &coefficients);
double euler_density_mass(const TriangleMesh &mesh, const std::vector<double> &coefficients);

} // namespace wasl
