// Troubled-cell characteristic WENO plus positivity limiting for Euler.

#pragma once

#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/limiting/euler_limiter_stencil.hpp"

namespace wasl {

// Limit coefficients (U[variable][mode][element]) in place. An element is
// troubled when the modal-decay ratio of its density exceeds threshold. Each
// troubled element blends its own polynomial with transformed neighbor
// candidates in characteristic variables, then applies the positivity limiter.
// All candidates read the unlimited input state, and every cell average is
// preserved. P0 is returned unchanged.
void limit_euler_modal(const EulerLimiterStencil &stencil, const Euler &equation, double threshold,
                       double density_floor, double pressure_floor,
                       std::vector<double> &coefficients);

} // namespace wasl
