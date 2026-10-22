#pragma once

#include <vector>

#include "wasl/core/simulation.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {

std::vector<double> initial_coefficients(const LineMesh &mesh, const ScalarAdvectionCase &config);
std::vector<double> initial_coefficients(const TriangleMesh &mesh,
                                         const ScalarAdvectionCase &config);
std::vector<double> initial_coefficients(const LineMesh &mesh, const EulerCase &config,
                                         const Euler &equation);
std::vector<double> initial_coefficients(const TriangleMesh &mesh, const EulerCase &config,
                                         const Euler &equation);

} // namespace wasl
