#pragma once

#include "wasl/core/simulation.hpp"

namespace wasl {

void validate_scalar_advection_case(const ScalarAdvectionCase &config);
void validate_euler_case(const EulerCase &config);

} // namespace wasl
