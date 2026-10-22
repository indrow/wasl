// Euler characteristic decomposition in a face-normal direction.
//
// Changes of the conserved state are projected onto the eigenvectors of the
// normal flux Jacobian frozen at a cell-average state. Waves are ordered
// (u_n - c, entropy, shear, u_n + c).

#pragma once

#include <array>
#include <cmath>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/equations/euler/primitive.hpp"

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
  const double sound_speed = wasl::sound_speed(equation, mean);
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

} // namespace wasl
