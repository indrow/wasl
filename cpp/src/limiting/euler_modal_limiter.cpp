#include "wasl/limiting/euler_modal_limiter.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "wasl/dg/modal_layout.hpp"
#include "wasl/limiting/positivity/euler_positivity.hpp"
#include "wasl/limiting/troubled_cell/modal_decay.hpp"
#include "wasl/limiting/weno/euler_characteristic_weno.hpp"

namespace wasl {
namespace {

constexpr std::size_t variable_count = Euler::variable_count;

// The sensor reads the density modes only.
bool is_troubled(const EulerLimiterStencil &stencil, const ModalLayout &layout,
                 const std::vector<double> &source, std::size_t element, double threshold) {
  const auto energy =
      modal_energy(source.data() + layout.index(0, 0, element), layout.element_count,
                   stencil.mode_count, stencil.highest_degree_begin);
  return !(energy.total == 0.0 || energy.highest_degree / energy.total <= threshold);
}

// Copy one element into packed [variable][mode].
void gather_element(const ModalLayout &layout, const std::vector<double> &source,
                    std::size_t element, double *packed) {
  for (std::size_t variable = 0; variable < variable_count; ++variable)
    for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
      packed[variable * layout.mode_count + mode] = source[layout.index(variable, mode, element)];
}

void scatter_element(const ModalLayout &layout, const double *packed, std::size_t element,
                     std::vector<double> &target) {
  for (std::size_t variable = 0; variable < variable_count; ++variable)
    for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
      target[layout.index(variable, mode, element)] = packed[variable * layout.mode_count + mode];
}

// Transform every neighbor into the element's basis, packed
// [candidate][variable][mode], with the element's own mean. Returns the count.
std::size_t build_candidates(const EulerLimiterStencil &stencil, const ModalLayout &layout,
                             const std::vector<double> &source, std::size_t element,
                             const double *target, double *candidates) {
  const std::size_t mode_count = layout.mode_count;
  std::size_t count = 0;
  for (int candidate = 0; candidate < stencil.maximum_neighbors; ++candidate) {
    const auto neighbor = stencil.neighbor(element, candidate);
    if (neighbor == invalid_element)
      continue;
    const double *matrix = stencil.transform(element, candidate);
    for (std::size_t variable = 0; variable < variable_count; ++variable) {
      double *result = candidates + (count * variable_count + variable) * mode_count;
      std::fill(result, result + mode_count, 0.0);
      for (std::size_t target_mode = 0; target_mode < mode_count; ++target_mode)
        for (std::size_t source_mode = 0; source_mode < mode_count; ++source_mode)
          result[target_mode] +=
              matrix[target_mode * mode_count + source_mode] *
              source[layout.index(variable, source_mode, static_cast<std::size_t>(neighbor))];
      result[0] = target[variable * mode_count];
    }
    ++count;
  }
  return count;
}

} // namespace

void limit_euler_modal(const EulerLimiterStencil &stencil, const Euler &equation, double threshold,
                       double density_floor, double pressure_floor,
                       std::vector<double> &coefficients) {
  if (stencil.polynomial_order == 0)
    return;
  const ModalLayout layout{stencil.mode_count, stencil.element_count};
  const auto source = coefficients;
  std::vector<double> target(variable_count * layout.mode_count);
  std::vector<double> candidates(static_cast<std::size_t>(stencil.maximum_neighbors) *
                                 variable_count * layout.mode_count);
  std::vector<double> limited(target.size());
  for (std::size_t element = 0; element < layout.element_count; ++element) {
    if (!is_troubled(stencil, layout, source, element, threshold))
      continue;
    gather_element(layout, source, element, target.data());
    const auto candidate_count =
        build_candidates(stencil, layout, source, element, target.data(), candidates.data());
    const EulerConserved average{target[0] / std::sqrt(2.0),
                                 target[layout.mode_count] / std::sqrt(2.0),
                                 target[2 * layout.mode_count] / std::sqrt(2.0),
                                 target[3 * layout.mode_count] / std::sqrt(2.0)};
    euler_characteristic_weno(equation, average, stencil.normal_x[element],
                              stencil.normal_y[element], target.data(), candidates.data(),
                              candidate_count, layout.mode_count, limited.data());
    limit_euler_positivity(equation, limited.data(), layout.mode_count, stencil.check_basis.data(),
                           stencil.check_point_count(), density_floor, pressure_floor);
    scatter_element(layout, limited.data(), element, coefficients);
  }
}

} // namespace wasl
