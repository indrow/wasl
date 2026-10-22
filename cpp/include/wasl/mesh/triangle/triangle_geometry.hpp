// Affine mapping between a physical straight-sided triangle and the reference
// triangle with vertices
//   v0 = (-1, -1),  v1 = ( 1, -1),  v2 = (-1,  1).
//
// The map sends reference vertices to physical vertices in order. Writing it as
// x(xi, eta) = p0 + A [xi+1, eta+1]^T / 2 with A = [p1 - p0 | p2 - p0], the
// constant Jacobian dX/d(xi,eta) is A/2 and det(A)/2 is the signed physical
// area: positive for counterclockwise vertex order.
//
// Element connectivity is normalized to positive orientation
// (NUMERICAL_CONVENTIONS.md). Reference edges carry outward normals
//   edge 0 (v0->v1): (0, -1)
//   edge 1 (v1->v2): ( 1,  1)/sqrt(2)
//   edge 2 (v2->v0): (-1,  0)
// which map to physical space through the inverse-transpose Jacobian.

#pragma once

#include <cmath>

namespace wasl {

// Length of reference edge 0, 1, or 2.
inline double reference_edge_length(int edge) noexcept {
  return edge == 1 ? 2.0 * std::sqrt(2.0) : 2.0;
}

struct TriangleMap {
  double x0, y0; // image of (-1, -1)
  double x1, y1; // image of ( 1, -1)
  double x2, y2; // image of (-1,  1)

  // Columns of A = [p1 - p0 | p2 - p0].
  double ax() const noexcept { return x1 - x0; }
  double ay() const noexcept { return y1 - y0; }
  double bx() const noexcept { return x2 - x0; }
  double by() const noexcept { return y2 - y0; }

  double det_a() const noexcept { return ax() * by() - bx() * ay(); }

  // det(dX/d(xi,eta)) = det(A)/4; one for the reference triangle itself.
  double jacobian_determinant() const noexcept { return det_a() / 4.0; }

  double area() const noexcept { return std::fabs(det_a()) / 2.0; }

  // +1 for counterclockwise vertex order, -1 otherwise.
  int orientation() const noexcept { return det_a() > 0.0 ? 1 : -1; }

  void reference_to_physical(double xi, double eta, double &x, double &y) const noexcept {
    x = x0 + 0.5 * (ax() * (xi + 1.0) + bx() * (eta + 1.0));
    y = y0 + 0.5 * (ay() * (xi + 1.0) + by() * (eta + 1.0));
  }

  // Inverse of the affine map; undefined for a degenerate triangle.
  void physical_to_reference(double x, double y, double &xi, double &eta) const noexcept {
    const double dx = x - x0;
    const double dy = y - y0;
    const double inverse_det = 1.0 / det_a();
    // Solve A [s; t] = [dx; dy] with s = (xi+1)/2, t = (eta+1)/2.
    const double s = (by() * dx - bx() * dy) * inverse_det;
    const double t = (-ay() * dx + ax() * dy) * inverse_det;
    xi = 2.0 * s - 1.0;
    eta = 2.0 * t - 1.0;
  }

  // det(J) J^-T applied to a reference gradient. Volume integrals use it so the
  // physical gradient and the area element combine without a division.
  void scaled_gradient(double gradient_xi, double gradient_eta, double &gradient_x,
                       double &gradient_y) const noexcept {
    gradient_x = 0.5 * (by() * gradient_xi - ay() * gradient_eta);
    gradient_y = 0.5 * (-bx() * gradient_xi + ax() * gradient_eta);
  }

  // Edge vector v(k) -> v(k+1) for local edge k.
  void edge_vector(int edge, double &ex, double &ey) const noexcept {
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
  }

  double edge_length(int edge) const noexcept {
    double ex = 0.0, ey = 0.0;
    edge_vector(edge, ex, ey);
    return std::hypot(ex, ey);
  }

  // Outward unit normal of edge k: the clockwise rotation (ty, -tx) of the unit
  // tangent for a positively oriented triangle, flipped otherwise.
  void physical_edge_normal(int edge, double &nx, double &ny) const noexcept {
    double ex = 0.0, ey = 0.0;
    edge_vector(edge, ex, ey);
    const double length = std::hypot(ex, ey);
    const double tx = ex / length;
    const double ty = ey / length;
    const int sign = orientation();
    nx = sign * ty;
    ny = sign * (-tx);
  }

  // Translate the triangle by a periodic image offset.
  void shift(double offset_x, double offset_y) noexcept {
    x0 += offset_x;
    x1 += offset_x;
    x2 += offset_x;
    y0 += offset_y;
    y1 += offset_y;
    y2 += offset_y;
  }
};

} // namespace wasl
