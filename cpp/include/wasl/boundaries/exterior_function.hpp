// Boundary policies for 1D scalar advection produce the exterior trace of every
// face from the current stage coefficients (U[element][mode]). The result is
// indexed by face and read by scalar_advection_residual. Policies recompute
// traces from whatever coefficients they receive, so they are valid at every
// Runge-Kutta stage.

#pragma once

#include <functional>
#include <vector>

namespace wasl {

using ExteriorFunction = std::function<std::vector<double>(const std::vector<double> &)>;

} // namespace wasl
