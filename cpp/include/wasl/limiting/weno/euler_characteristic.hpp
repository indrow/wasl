// Euler characteristic projection and simple WENO in a face-normal direction.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/limiting/weno/simple_weno.hpp"

namespace wasl {

using EulerCharacteristic = std::array<double, Euler::variable_count>;

inline EulerCharacteristic euler_to_characteristic(const Euler &equation,
                                                   const EulerConserved &average,
                                                   const EulerConserved &change, double normal_x,
                                                   double normal_y) noexcept {
  const EulerPrimitive mean = primitive(equation, average);
  const double normal_velocity = mean.velocity_x * normal_x + mean.velocity_y * normal_y;
  const double tangent_velocity = -mean.velocity_x * normal_y + mean.velocity_y * normal_x;
  const double normal_momentum = change.momentum_x * normal_x + change.momentum_y * normal_y;
  const double tangent_momentum = -change.momentum_x * normal_y + change.momentum_y * normal_x;
  const double sound_speed_squared = equation.heat_capacity_ratio * mean.pressure / mean.density;
  const double sound_speed = std::sqrt(sound_speed_squared);
  const double kinetic_energy =
      0.5 * (normal_velocity * normal_velocity + tangent_velocity * tangent_velocity);
  const double factor = (equation.heat_capacity_ratio - 1.0) / sound_speed_squared;

  const double acoustic_common =
      factor * (kinetic_energy * change.density - normal_velocity * normal_momentum -
                tangent_velocity * tangent_momentum + change.total_energy);
  const double acoustic_difference =
      (normal_velocity * change.density - normal_momentum) / sound_speed;
  return {0.5 * (acoustic_common + acoustic_difference), change.density - acoustic_common,
          tangent_momentum - tangent_velocity * change.density,
          0.5 * (acoustic_common - acoustic_difference)};
}

inline EulerConserved euler_from_characteristic(const Euler &equation,
                                                const EulerConserved &average,
                                                const EulerCharacteristic &characteristic,
                                                double normal_x, double normal_y) noexcept {
  const EulerPrimitive mean = primitive(equation, average);
  const double normal_velocity = mean.velocity_x * normal_x + mean.velocity_y * normal_y;
  const double tangent_velocity = -mean.velocity_x * normal_y + mean.velocity_y * normal_x;
  const double sound_speed = std::sqrt(equation.heat_capacity_ratio * mean.pressure / mean.density);
  const double enthalpy = (average.total_energy + mean.pressure) / average.density;
  const double kinetic_energy =
      0.5 * (normal_velocity * normal_velocity + tangent_velocity * tangent_velocity);

  const double density = characteristic[0] + characteristic[1] + characteristic[3];
  const double normal_momentum = (normal_velocity - sound_speed) * characteristic[0] +
                                 normal_velocity * characteristic[1] +
                                 (normal_velocity + sound_speed) * characteristic[3];
  const double tangent_momentum = tangent_velocity * density + characteristic[2];
  const double total_energy = (enthalpy - normal_velocity * sound_speed) * characteristic[0] +
                              kinetic_energy * characteristic[1] +
                              tangent_velocity * characteristic[2] +
                              (enthalpy + normal_velocity * sound_speed) * characteristic[3];
  return {density, normal_x * normal_momentum - normal_y * tangent_momentum,
          normal_y * normal_momentum + normal_x * tangent_momentum, total_energy};
}

// Arrays use [variable][mode], and candidates use [candidate][variable][mode].
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
