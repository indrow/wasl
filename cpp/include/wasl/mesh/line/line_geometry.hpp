// Affine mapping between a physical line element and the reference line.
//
// A line element has physical endpoints x_left < x_right after orientation
// normalization. The reference coordinate xi in [-1, 1] maps to physical space
// by
//   x(xi) = x_centre + (length / 2) * xi.
// The Jacobian dx / d xi = length / 2 is constant, so the element is affine.
// The left face (xi = -1) carries outward normal -1 and the right face
// (xi = +1) carries outward normal +1.

#pragma once

namespace wasl {

struct LineMap {
  double x_left;
  double x_right;

  double length() const noexcept { return x_right - x_left; }
  double jacobian() const noexcept { return length() / 2.0; }
  double centre() const noexcept { return (x_left + x_right) / 2.0; }

  double reference_to_physical(double xi) const noexcept { return centre() + jacobian() * xi; }
  double physical_to_reference(double x) const noexcept { return (x - centre()) / jacobian(); }

  // Translate the element by a periodic image offset.
  void shift(double offset) noexcept {
    x_left += offset;
    x_right += offset;
  }

  static constexpr double left_face_normal = -1.0;
  static constexpr double right_face_normal = 1.0;
};

} // namespace wasl
