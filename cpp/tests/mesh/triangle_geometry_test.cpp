// Tests for the affine triangle geometry mapping.
//
// Checks: reference-to-physical vertex images, round-trip
// reference->physical->reference, constant Jacobian, positive and negative
// orientation, area matching the shoelace formula, and physical edge outward
// normals and lengths.

#include <cmath>

#include "wasl/mesh/triangle_geometry.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

wasl::TriangleMap reference_triangle() { return {-1.0, -1.0, 1.0, -1.0, -1.0, 1.0}; }

wasl::TriangleMap physical_triangle() {
  // A generic scalene, positively-oriented triangle.
  return {0.0, 0.0, 4.0, 0.0, 1.0, 3.0};
}

// Reference vertices map to the stored physical vertices.
bool check_vertex_images() {
  const auto m = physical_triangle();
  const double pts[][2] = {{-1.0, -1.0}, {1.0, -1.0}, {-1.0, 1.0}};
  const double expected[][2] = {{m.x0, m.y0}, {m.x1, m.y1}, {m.x2, m.y2}};
  for (int i = 0; i < 3; ++i) {
    double x = 0.0, y = 0.0;
    m.reference_to_physical(pts[i][0], pts[i][1], x, y);
    if (!approx(x, expected[i][0]) || !approx(y, expected[i][1]))
      return false;
  }
  return true;
}

// Round-trip physical->reference->physical recovers the input.
bool check_round_trip() {
  const auto m = physical_triangle();
  const double probe[][2] = {{0.5, 0.2}, {1.5, 0.7}, {2.0, 1.0}, {0.2, 1.8}, {1.0, 0.5}};
  for (const auto &p : probe) {
    double xi = 0.0, eta = 0.0;
    m.physical_to_reference(p[0], p[1], xi, eta);
    double x = 0.0, y = 0.0;
    m.reference_to_physical(xi, eta, x, y);
    if (!approx(x, p[0]) || !approx(y, p[1]))
      return false;
  }
  return true;
}

// Jacobian determinant is constant (2x signed area) and matches the shoelace
// formula. Orientation is +1 for counterclockwise vertex order.
bool check_jacobian_and_orientation() {
  const auto m = physical_triangle();
  // Shoelace double-area = x0(y1-y2) + x1(y2-y0) + x2(y0-y1) equals det(A).
  // The metric Jacobian det(dX/d(xi,eta)) = det(A)/4; area = |det(A)|/2.
  const double signed_double = m.x0 * (m.y1 - m.y2) + m.x1 * (m.y2 - m.y0) + m.x2 * (m.y0 - m.y1);
  if (!approx(m.jacobian_determinant(), signed_double / 4.0))
    return false;
  if (!approx(m.area(), std::fabs(signed_double) / 2.0))
    return false;
  if (m.orientation() != 1)
    return false;

  // Flipped (clockwise) triangle: negative orientation, negative det.
  wasl::TriangleMap flipped{m.x0, m.y0, m.x2, m.y2, m.x1, m.y1};
  if (flipped.orientation() != -1)
    return false;
  if (flipped.jacobian_determinant() >= 0.0)
    return false;
  // Areas agree in magnitude.
  if (!approx(flipped.area(), m.area()))
    return false;
  return true;
}

// The reference triangle is the identity map: det = 1, area = 2.
bool check_reference_identity() {
  const auto m = reference_triangle();
  if (!approx(m.jacobian_determinant(), 1.0))
    return false;
  if (!approx(m.area(), 2.0))
    return false;
  if (m.orientation() != 1)
    return false;
  double x = 0.0, y = 0.0;
  m.reference_to_physical(0.0, -0.5, x, y);
  if (!approx(x, 0.0) || !approx(y, -0.5))
    return false;
  return true;
}

// Physical edge outward normals are unit length and point outward. For the
// positively oriented triangle, the bottom edge normal points in -y.
bool check_edge_normals() {
  const auto m = physical_triangle();
  for (int e = 0; e < 3; ++e) {
    double nx = 0.0, ny = 0.0;
    m.physical_edge_normal(e, nx, ny);
    if (!approx(nx * nx + ny * ny, 1.0))
      return false;
  }
  // Bottom edge (v0->v1 along +x): outward normal is (0, -1).
  double nx = 0.0, ny = 0.0;
  m.physical_edge_normal(0, nx, ny);
  if (!approx(nx, 0.0) || !approx(ny, -1.0))
    return false;
  return true;
}

// Edge lengths match the Euclidean distance between their physical endpoints.
bool check_edge_lengths() {
  const auto m = physical_triangle();
  const double expected[] = {
      std::hypot(m.x1 - m.x0, m.y1 - m.y0),
      std::hypot(m.x2 - m.x1, m.y2 - m.y1),
      std::hypot(m.x0 - m.x2, m.y0 - m.y2),
  };
  for (int e = 0; e < 3; ++e) {
    if (!approx(m.edge_length(e), expected[e]))
      return false;
  }
  return true;
}

} // namespace

int main() {
  if (!check_reference_identity())
    return 1;
  if (!check_vertex_images())
    return 1;
  if (!check_round_trip())
    return 1;
  if (!check_jacobian_and_orientation())
    return 1;
  if (!check_edge_normals())
    return 1;
  if (!check_edge_lengths())
    return 1;
  return 0;
}
