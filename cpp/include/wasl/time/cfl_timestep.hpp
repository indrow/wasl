// Stable explicit time steps for modal DG with SSP-RK3.
//
// Every limit has the form dt = cfl * h / (lambda * (2p + 1)), where h is the
// smallest element length (line) or altitude (triangle), lambda the largest
// characteristic speed, and (2p + 1) the DG spectral-radius scaling. A zero
// speed returns 0.0, meaning no CFL restriction.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

double advection_cfl_timestep(const LineMesh &mesh, const ScalarAdvection &eq, int polynomial_order,
                              double cfl);

double triangle_scalar_cfl_timestep(const TriangleMesh &mesh, double velocity_x, double velocity_y,
                                    int polynomial_order, double cfl);

// Euler limits use cell-average states from U[variable][mode][element] and
// throw when an average is not physical.
double euler_cfl_timestep(const LineMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl);
double euler_cfl_timestep(const TriangleMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl);

} // namespace wasl
