// Affine mapping between a physical line element and the reference line.
//
// A line element has physical endpoints x_left and x_right (x_left < x_right
// after orientation normalization). The reference coordinate xi in [-1, 1] maps
// to physical space by
//   x(xi) = x_centre + (length / 2) * xi,
//   x_centre = (x_left + x_right) / 2,  length = x_right - x_left.
// The scalar Jacobian dx / d xi = length / 2 is constant, so the element is
// affine. The left face (xi = -1) carries outward normal -1 and the right face
// (xi = +1) carries outward normal +1.

#pragma once

namespace wasl {

struct LineMap {
  double x_left;
  double x_right;

  // Physical element length (always nonnegative for a normalized element).
  double length() const noexcept { return x_right - x_left; }

  // Reference-to-physical Jacobian dx / d xi, constant over the element.
  double jacobian() const noexcept { return length() / 2.0; }

  // Physical centre of the element.
  double centre() const noexcept { return (x_left + x_right) / 2.0; }

  // Map a reference coordinate to physical space.
  double reference_to_physical(double xi) const noexcept { return centre() + jacobian() * xi; }

  // Map a physical coordinate to reference space (inverse of the affine map).
  double physical_to_reference(double x) const noexcept { return (x - centre()) / jacobian(); }

  // Outward normal at a face: -1 at the left face (xi = -1), +1 at the right
  // face (xi = +1).
  static constexpr double left_face_normal = -1.0;
  static constexpr double right_face_normal = 1.0;
};

} // namespace wasl
