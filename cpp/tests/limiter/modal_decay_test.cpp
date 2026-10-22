#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "wasl/limiting/troubled_cell/modal_decay.hpp"

int main() {
  // P2 modes are ordered by degree: 1 constant, 2 linear, 3 quadratic.
  const std::array<double, 6> smooth{1.0, 0.1, -0.08, 1e-3, -8e-4, 6e-4};
  const std::array<double, 6> discontinuous{1.0, 0.2, -0.1, 0.4, -0.35, 0.3};
  constexpr double threshold = 1e-3;
  if (wasl::modal_decay_troubled(smooth.data(), 2, threshold) ||
      !wasl::modal_decay_troubled(discontinuous.data(), 2, threshold))
    return 1;

  const std::array<double, 1> constant{2.0};
  if (wasl::modal_decay_troubled(constant.data(), 0, threshold))
    return 1;

  const std::array<std::uint8_t, 7> flags{0, 1, 0, 1, 1, 0, 1};
  const std::vector<std::uint32_t> expected{1, 3, 4, 6};
  if (wasl::compact_troubled_elements(flags.data(), flags.size()) != expected ||
      !wasl::compact_troubled_elements(flags.data(), 0).empty())
    return 1;
  return 0;
}
