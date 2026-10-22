// SSP-RK3 step for 1D scalar advection with a boundary policy.

#pragma once

#include <vector>

#include "wasl/boundaries/exterior_function.hpp"
#include "wasl/dg/scalar_advection/line_residual.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

// exterior recomputes face exterior traces from each stage's coefficients.
inline void ssp_rk3_step(const LineMesh &mesh, const ScalarAdvection &eq, int polynomial_order,
                         double dt, const ExteriorFunction &exterior,
                         std::vector<double> &coefficients, SspRk3Workspace &workspace) {
  const auto residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    workspace.exterior = exterior(stage);
    scalar_advection_residual(mesh, eq, polynomial_order, stage.data(), workspace.exterior.data(),
                              rhs.data());
  };
  ssp_rk3_step(dt, residual, coefficients, workspace);
}

} // namespace wasl
