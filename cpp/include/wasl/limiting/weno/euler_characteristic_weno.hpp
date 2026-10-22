// Characteristic-wise simple WENO for Euler.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/limiting/weno/euler_characteristic.hpp"
#include "wasl/limiting/weno/simple_weno.hpp"

namespace wasl {

// Limit in characteristic variables of the given normal, then transform back.
// Arrays use [variable][mode] and candidates [candidate][variable][mode]. Mode
// zero of every variable is copied from the target to preserve cell averages.
inline void euler_characteristic_weno(const Euler &equation, const EulerConserved &average,
                                      double normal_x, double normal_y, const double *target,
                                      const double *candidates, std::size_t candidate_count,
                                      std::size_t mode_count, double *limited) {
  std::vector<double> target_characteristic(Euler::variable_count * mode_count);
  std::vector<double> candidate_characteristic(candidate_count * Euler::variable_count *
                                               mode_count);
  std::vector<double> limited_characteristic(Euler::variable_count * mode_count);

  const auto project = [&](const double *source, double *destination) {
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      const EulerConserved change{source[mode], source[mode_count + mode],
                                  source[2 * mode_count + mode], source[3 * mode_count + mode]};
      const EulerCharacteristic value =
          euler_to_characteristic(equation, average, change, normal_x, normal_y);
      for (int wave = 0; wave < Euler::variable_count; ++wave)
        destination[static_cast<std::size_t>(wave) * mode_count + mode] = value[wave];
    }
  };
  project(target, target_characteristic.data());
  for (std::size_t candidate = 0; candidate < candidate_count; ++candidate)
    project(candidates + candidate * Euler::variable_count * mode_count,
            candidate_characteristic.data() + candidate * Euler::variable_count * mode_count);

  for (int wave = 0; wave < Euler::variable_count; ++wave) {
    const std::size_t offset = static_cast<std::size_t>(wave) * mode_count;
    std::vector<double> wave_candidates(candidate_count * mode_count);
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        wave_candidates[candidate * mode_count + mode] =
            candidate_characteristic[candidate * Euler::variable_count * mode_count + offset +
                                     mode];
    simple_weno(target_characteristic.data() + offset, wave_candidates.data(), candidate_count,
                mode_count, limited_characteristic.data() + offset);
  }

  for (std::size_t mode = 0; mode < mode_count; ++mode) {
    EulerCharacteristic value{};
    for (int wave = 0; wave < Euler::variable_count; ++wave)
      value[wave] = limited_characteristic[static_cast<std::size_t>(wave) * mode_count + mode];
    const EulerConserved state =
        euler_from_characteristic(equation, average, value, normal_x, normal_y);
    limited[mode] = state.density;
    limited[mode_count + mode] = state.momentum_x;
    limited[2 * mode_count + mode] = state.momentum_y;
    limited[3 * mode_count + mode] = state.total_energy;
  }
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    limited[static_cast<std::size_t>(variable) * mode_count] =
        target[static_cast<std::size_t>(variable) * mode_count];
}

} // namespace wasl
