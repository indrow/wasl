// L2 projection of a case's initial condition.

#pragma once

#include <vector>

#include "wasl/case/euler_case.hpp"
#include "wasl/case/scalar_advection_case.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// Scalar results use the solver layout U[element][mode].
std::vector<double> initial_coefficients(const LineMesh &mesh, const ScalarAdvectionCase &config);
std::vector<double> initial_coefficients(const TriangleMesh &mesh,
                                         const ScalarAdvectionCase &config);

// Euler results use U[variable][mode][element].
std::vector<double> initial_coefficients(const LineMesh &mesh, const EulerCase &config,
                                         const Euler &equation);
std::vector<double> initial_coefficients(const TriangleMesh &mesh, const EulerCase &config,
                                         const Euler &equation);

} // namespace wasl
