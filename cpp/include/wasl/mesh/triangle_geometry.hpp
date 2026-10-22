// Affine mapping between a physical straight-sided triangle and the reference
// triangle with vertices
//   v0 = (-1, -1),  v1 = ( 1, -1),  v2 = (-1,  1).
//
// The affine map sends reference vertices to physical vertices in order:
//   x(-1,-1) = p0,  x(1,-1) = p1,  x(-1,1) = p2.
// Writing the map as x(xi, eta) = p0 + A [xi+1, eta+1]^T / 2 with
//   A = [p1 - p0 | p2 - p0],
// the constant 2x2 Jacobian dX/d(xi,eta) is A/2 and its determinant is the
// signed physical area: |det A| / 2 is the area, sign(det A) is +1 for
// counterclockwise (positive) vertex ordering and -1 for clockwise.
//
// Element connectivity must be normalized to positive physical orientation
// (NUMERICAL_CONVENTIONS.md). The reference edges carry outward normals
//   edge 0 (v0->v1): (0, -1)
//   edge 1 (v1->v2): ( 1,  1)/sqrt(2)
//   edge 2 (v2->v0): (-1,  0)
// which map to physical space through the inverse-transpose Jacobian.

#pragma once

#include <array>
#include <cmath>
#include <numbers>

namespace wasl {

struct TriangleMap {
  double x0, y0; // physical vertex 0 (image of (-1, -1))
  double x1, y1; // physical vertex 1 (image of ( 1, -1))
  double x2, y2; // physical vertex 2 (image of (-1,  1))

  // Columns of A = [p1 - p0 | p2 - p0].
  double ax() const noexcept { return x1 - x0; } // A_{0,0}
  double ay() const noexcept { return y1 - y0; } // A_{1,0}
  double bx() const noexcept { return x2 - x0; } // A_{0,1}
  double by() const noexcept { return y2 - y0; } // A_{1,1}

  // det(A) = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0). Positive iff the
  // physical vertices are ordered counterclockwise.
  double det_a() const noexcept { return ax() * by() - bx() * ay(); }

  // Reference-to-physical Jacobian determinant det(dX/d(xi,eta)). The Jacobian
  // matrix is A/2, so its determinant is det(A)/4. For the reference triangle
  // det(A) = 4, giving the identity Jacobian determinant 1. Its sign gives the
  // physical orientation; its magnitude scales reference integrals to physical.
  double jacobian_determinant() const noexcept { return det_a() / 4.0; }

  // Physical area = |det(J)| * area(reference triangle) = |det(A)| / 2.
  double area() const noexcept { return std::fabs(det_a()) / 2.0; }

  // +1 for positive (counterclockwise) physical orientation, -1 otherwise.
  int orientation() const noexcept { return det_a() > 0.0 ? 1 : -1; }

  // Map a reference coordinate to physical space.
  void reference_to_physical(double xi, double eta, double &x, double &y) const noexcept {
    x = x0 + 0.5 * (ax() * (xi + 1.0) + bx() * (eta + 1.0));
    y = y0 + 0.5 * (ay() * (xi + 1.0) + by() * (eta + 1.0));
  }

  // Map a physical coordinate back to reference space (inverse of the affine
  // map). Undefined for a degenerate (zero-area) triangle.
  void physical_to_reference(double x, double y, double &xi, double &eta) const noexcept {
    const double dx = x - x0;
    const double dy = y - y0;
    const double inv_det = 1.0 / det_a();
    // Solve A [s; t] = [dx; dy] with s = (xi+1)/2, t = (eta+1)/2.
    const double s = (by() * dx - bx() * dy) * inv_det;
    const double t = (-ay() * dx + ax() * dy) * inv_det;
    xi = 2.0 * s - 1.0;
    eta = 2.0 * t - 1.0;
  }

  // The three physical edge tangent directions, following the reference edge
  // vertex pairs v0->v1, v1->v2, v2->v0. Each tangent is unit length.
  void edge_tangent(int edge, double &tx, double &ty) const noexcept {
    double ex = 0.0, ey = 0.0;
    if (edge == 0) {
      ex = x1 - x0;
      ey = y1 - y0;
    } else if (edge == 1) {
      ex = x2 - x1;
      ey = y2 - y1;
    } else {
      ex = x0 - x2;
      ey = y0 - y2;
    }
    const double len = std::hypot(ex, ey);
    tx = ex / len;
    ty = ey / len;
  }

  // Physical outward unit normal at edge k. For a positively-oriented triangle
  // the outward normal is the 90-degree clockwise rotation of the (reference)
  // edge tangent v(i)->v(i+1): (ty, -tx). This points away from the interior.
  void physical_edge_normal(int edge, double &nx, double &ny) const noexcept {
    double tx = 0.0, ty = 0.0;
    edge_tangent(edge, tx, ty);
    const int sgn = orientation();
    // Clockwise rotation of the tangent is (ty, -tx); flip for negative orient.
    nx = sgn * ty;
    ny = sgn * (-tx);
  }

  // Physical edge length.
  double edge_length(int edge) const noexcept {
    double ex = 0.0, ey = 0.0;
    if (edge == 0) {
      ex = x1 - x0;
      ey = y1 - y0;
    } else if (edge == 1) {
      ex = x2 - x1;
      ey = y2 - y1;
    } else {
      ex = x0 - x2;
      ey = y0 - y2;
    }
    return std::hypot(ex, ey);
  }
};

// Reference outward unit normals for the three reference edges.
inline constexpr std::array<double, 6> reference_edge_normals() {
  // edge 0: (0, -1), edge 1: (1,1)/sqrt2, edge 2: (-1, 0)
  return {0.0, -1.0, 1.0 / std::numbers::sqrt2, 1.0 / std::numbers::sqrt2, -1.0, 0.0};
}

} // namespace wasl
