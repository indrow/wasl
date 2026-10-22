#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/limiting/weno/candidate_transform.hpp"
#include "wasl/limiting/weno/simple_weno.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace {

constexpr int kPolynomialOrder = 2;
constexpr std::size_t kModeCount = wasl::triangle_mode_count(kPolynomialOrder);

std::array<double, kModeCount> project_linear(const wasl::TriangleMap &map) {
  std::array<double, kModeCount> coefficients{};
  std::array<double, kModeCount> basis{};
  const auto quadrature = wasl::triangle_volume_rule(kPolynomialOrder);
  for (int q = 0; q < quadrature.size(); ++q) {
    double x = 0.0, y = 0.0;
    map.reference_to_physical(quadrature.xi[q], quadrature.eta[q], x, y);
    wasl::dubiner_basis(quadrature.xi[q], quadrature.eta[q], kPolynomialOrder, basis.data());
    const double value = 2.0 + 0.5 * x - 0.25 * y;
    for (std::size_t mode = 0; mode < kModeCount; ++mode)
      coefficients[mode] += quadrature.weights[q] * value * basis[mode];
  }
  return coefficients;
}

double high_mode_energy(const double *coefficients) {
  double energy = 0.0;
  for (std::size_t mode = 3; mode < kModeCount; ++mode)
    energy += coefficients[mode] * coefficients[mode];
  return energy;
}

} // namespace

int main() {
  const wasl::TriangleMap target_map{0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
  const wasl::TriangleMap neighbor_map{1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
  const auto target_linear = project_linear(target_map);
  const auto neighbor_linear = project_linear(neighbor_map);
  std::array<double, kModeCount> transformed{};
  wasl::project_neighbor_polynomial(neighbor_map, target_map, kPolynomialOrder,
                                    neighbor_linear.data(), transformed.data());
  for (std::size_t mode = 0; mode < kModeCount; ++mode) {
    if (std::fabs(transformed[mode] - target_linear[mode]) > 2e-13)
      return 1;
  }

  std::array<double, kModeCount> candidate{};
  wasl::make_weno_candidate(neighbor_map, target_map, kPolynomialOrder, neighbor_linear.data(),
                            target_linear[0], candidate.data());
  if (candidate[0] != target_linear[0])
    return 1;

  std::array<double, kModeCount> smooth_limited{};
  std::array<double, 2 * kModeCount> smooth_candidates{};
  std::copy(target_linear.begin(), target_linear.end(), smooth_candidates.begin());
  std::copy(target_linear.begin(), target_linear.end(), smooth_candidates.begin() + kModeCount);
  wasl::simple_weno(target_linear.data(), smooth_candidates.data(), 2, kModeCount,
                    smooth_limited.data());
  for (std::size_t mode = 0; mode < kModeCount; ++mode) {
    if (std::fabs(smooth_limited[mode] - target_linear[mode]) > 2e-15)
      return 1;
  }

  const std::array<double, kModeCount> discontinuous{1.4, 0.1, -0.1, 0.8, -0.7, 0.6};
  const std::array<double, 2 * kModeCount> calm_candidates{1.4, 0.02,  -0.01, 0.0, 0.0, 0.0,
                                                           1.4, -0.01, 0.02,  0.0, 0.0, 0.0};
  std::array<double, kModeCount> limited{};
  wasl::simple_weno(discontinuous.data(), calm_candidates.data(), 2, kModeCount, limited.data());
  if (limited[0] != discontinuous[0] ||
      !(high_mode_energy(limited.data()) < 1e-8 * high_mode_energy(discontinuous.data())))
    return 1;
  return 0;
}
