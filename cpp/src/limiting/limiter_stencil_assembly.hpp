// Helpers shared by the line and triangle limiter-stencil builders.

#pragma once

#include <cstddef>
#include <cstdint>

#include "wasl/limiting/euler_limiter_stencil.hpp"
#include "wasl/mesh/mesh_index.hpp"

namespace wasl {

inline void allocate_neighbors(EulerLimiterStencil &stencil) {
  stencil.neighbors.assign(
      stencil.element_count * static_cast<std::size_t>(stencil.maximum_neighbors), invalid_element);
  stencil.normal_x.assign(stencil.element_count, 1.0);
  stencil.normal_y.assign(stencil.element_count, 0.0);
}

inline void allocate_transforms(EulerLimiterStencil &stencil) {
  stencil.transforms.assign(stencil.element_count *
                                static_cast<std::size_t>(stencil.maximum_neighbors) *
                                stencil.mode_count * stencil.mode_count,
                            0.0);
}

// Record a neighbor once, in the first free candidate slot.
inline void append_unique_neighbor(EulerLimiterStencil &stencil, std::size_t element,
                                   std::int32_t neighbor) {
  auto *values =
      stencil.neighbors.data() + element * static_cast<std::size_t>(stencil.maximum_neighbors);
  for (int index = 0; index < stencil.maximum_neighbors; ++index) {
    if (values[index] == neighbor)
      return;
    if (values[index] == invalid_element) {
      values[index] = neighbor;
      return;
    }
  }
}

// Offset that brings a neighbor whose centre lies more than half a period away
// next to its target, so the transform evaluates the neighbor polynomial in its
// own element.
inline double periodic_offset(double gap, double period) noexcept {
  return gap > 0.5 * period ? period : gap < -0.5 * period ? -period : 0.0;
}

} // namespace wasl
