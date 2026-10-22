// Validation of case configurations at the bridge boundary.

#pragma once

#include "wasl/case/euler_case.hpp"
#include "wasl/case/scalar_advection_case.hpp"

namespace wasl {

// Throw std::runtime_error describing the first invalid setting.
void validate_scalar_advection_case(const ScalarAdvectionCase &config);
void validate_euler_case(const EulerCase &config);

} // namespace wasl
