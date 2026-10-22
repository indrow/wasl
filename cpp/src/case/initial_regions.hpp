// Point evaluation of configured initial conditions.

#pragma once

#include "wasl/case/euler_case.hpp"
#include "wasl/case/scalar_advection_case.hpp"

namespace wasl {

double scalar_initial_value(const ScalarAdvectionCase &config, double x);
double scalar_initial_value(const ScalarAdvectionCase &config, double x, double y);

// Throws when no default region exists.
EulerPrimitive euler_initial_value(const EulerCase &config, double x, double y = 0.0);

} // namespace wasl
