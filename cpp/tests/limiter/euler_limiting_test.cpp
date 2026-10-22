#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/limiting/positivity/euler_positivity.hpp"
#include "wasl/limiting/weno/euler_characteristic.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace {

constexpr wasl::Euler kEuler{1.4};
constexpr int kPolynomialOrder = 2;
constexpr std::size_t kModeCount = wasl::triangle_mode_count(kPolynomialOrder);

std::vector<double> check_basis() {
  const auto volume = wasl::triangle_volume_rule(kPolynomialOrder);
  const auto edges = wasl::triangle_edge_rules(kPolynomialOrder);
  std::vector<double> values;
  std::array<double, kModeCount> basis{};
  const auto append = [&](double xi, double eta) {
    wasl::dubiner_basis(xi, eta, kPolynomialOrder, basis.data());
    values.insert(values.end(), basis.begin(), basis.end());
  };
  for (int point = 0; point < volume.size(); ++point)
    append(volume.xi[point], volume.eta[point]);
  for (const auto &edge : edges)
    for (std::size_t point = 0; point < edge.xi.size(); ++point)
      append(edge.xi[point], edge.eta[point]);
  return values;
}

} // namespace

int main() {
  const wasl::EulerConserved average =
      wasl::conserved(kEuler, wasl::EulerPrimitive{1.2, 0.7, -0.3, 1.1});
  const wasl::EulerConserved change{0.13, -0.21, 0.17, 0.42};
  const double normal_x = 0.6;
  const double normal_y = -0.8;
  const auto characteristic =
      wasl::euler_to_characteristic(kEuler, average, change, normal_x, normal_y);
  const auto recovered =
      wasl::euler_from_characteristic(kEuler, average, characteristic, normal_x, normal_y);
  const std::array<double, 4> original{change.density, change.momentum_x, change.momentum_y,
                                       change.total_energy};
  const std::array<double, 4> round_trip{recovered.density, recovered.momentum_x,
                                         recovered.momentum_y, recovered.total_energy};
  for (std::size_t variable = 0; variable < original.size(); ++variable)
    if (std::fabs(round_trip[variable] - original[variable]) > 2e-15)
      return 1;

  std::array<double, 4 * kModeCount> target{};
  std::array<double, 2 * 4 * kModeCount> candidates{};
  const std::array<double, 4> means{average.density, average.momentum_x, average.momentum_y,
                                    average.total_energy};
  for (std::size_t variable = 0; variable < 4; ++variable) {
    target[variable * kModeCount] = means[variable] * std::sqrt(2.0);
    target[variable * kModeCount + 1] = 0.1 * (variable + 1.0);
    for (std::size_t candidate = 0; candidate < 2; ++candidate)
      std::copy_n(target.data() + variable * kModeCount, kModeCount,
                  candidates.data() + candidate * 4 * kModeCount + variable * kModeCount);
  }
  std::array<double, 4 * kModeCount> limited{};
  wasl::euler_characteristic_weno(kEuler, average, normal_x, normal_y, target.data(),
                                  candidates.data(), 2, kModeCount, limited.data());
  for (std::size_t index = 0; index < target.size(); ++index)
    if (std::fabs(limited[index] - target[index]) > 3e-15)
      return 1;

  for (std::size_t candidate = 0; candidate < 2; ++candidate)
    for (std::size_t variable = 0; variable < 4; ++variable)
      std::fill_n(candidates.data() + candidate * 4 * kModeCount + variable * kModeCount + 1,
                  kModeCount - 1, 0.0);
  wasl::euler_characteristic_weno(kEuler, average, normal_x, normal_y, target.data(),
                                  candidates.data(), 2, kModeCount, limited.data());
  double target_energy = 0.0;
  double limited_energy = 0.0;
  for (std::size_t variable = 0; variable < 4; ++variable) {
    if (limited[variable * kModeCount] != target[variable * kModeCount])
      return 1;
    for (std::size_t mode = 1; mode < kModeCount; ++mode) {
      target_energy += target[variable * kModeCount + mode] * target[variable * kModeCount + mode];
      limited_energy +=
          limited[variable * kModeCount + mode] * limited[variable * kModeCount + mode];
    }
  }
  if (!(limited_energy < 1e-8 * target_energy))
    return 1;

  const std::vector<double> basis = check_basis();
  const std::size_t point_count = basis.size() / kModeCount;
  std::array<double, 4 * kModeCount> polynomial{};
  for (std::size_t variable = 0; variable < 4; ++variable)
    polynomial[variable * kModeCount] = means[variable] * std::sqrt(2.0);
  polynomial[1] = 1.8;
  polynomial[kModeCount + 2] = 2.5;
  polynomial[2 * kModeCount + 1] = -2.0;
  polynomial[3 * kModeCount + 2] = -2.8;
  const std::array<double, 4> mode_zero{polynomial[0], polynomial[kModeCount],
                                        polynomial[2 * kModeCount], polynomial[3 * kModeCount]};
  constexpr double density_floor = 0.1;
  constexpr double pressure_floor = 0.08;
  const double scale =
      wasl::limit_euler_positivity(kEuler, polynomial.data(), kModeCount, basis.data(), point_count,
                                   density_floor, pressure_floor);
  if (!(scale >= 0.0 && scale < 1.0))
    return 1;
  for (std::size_t variable = 0; variable < 4; ++variable)
    if (polynomial[variable * kModeCount] != mode_zero[variable])
      return 1;
  for (std::size_t point = 0; point < point_count; ++point) {
    const auto state = wasl::evaluate_euler_modal(polynomial.data(), kModeCount,
                                                  basis.data() + point * kModeCount);
    if (state.density < density_floor || wasl::primitive(kEuler, state).pressure < pressure_floor)
      return 1;
  }

  std::array<double, 4 * kModeCount> constant{};
  for (std::size_t variable = 0; variable < 4; ++variable)
    constant[variable * kModeCount] = means[variable] * std::sqrt(2.0);
  if (wasl::limit_euler_positivity(kEuler, constant.data(), kModeCount, basis.data(), point_count,
                                   density_floor, pressure_floor) != 1.0)
    return 1;
  return 0;
}
