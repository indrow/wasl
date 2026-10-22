#pragma once

#include <cstdint>
#include <vector>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {

// Precomputed periodic donor elements and basis values for triangle boundary faces.
struct PeriodicTriangleTrace {
  std::vector<std::int32_t> elements;
  std::vector<double> basis;
  int point_count = 0;
};

PeriodicTriangleTrace periodic_triangle_trace(const TriangleMesh &mesh, int polynomial_order);

std::vector<double> periodic_triangle_exterior(const PeriodicTriangleTrace &trace,
                                               const std::vector<double> &coefficients,
                                               std::size_t mode_count);

std::vector<double> periodic_line_euler_exterior(const LineMesh &mesh, int polynomial_order,
                                                 const std::vector<double> &coefficients);

std::vector<double> periodic_triangle_euler_exterior(const PeriodicTriangleTrace &trace,
                                                     const std::vector<double> &coefficients,
                                                     std::size_t element_count,
                                                     std::size_t mode_count,
                                                     std::size_t face_count);

double triangle_scalar_cfl_timestep(const TriangleMesh &mesh, double velocity_x, double velocity_y,
                                    int polynomial_order, double cfl);

double euler_cfl_timestep(const LineMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl);
double euler_cfl_timestep(const TriangleMesh &mesh, const Euler &equation,
                          const std::vector<double> &coefficients, std::size_t mode_count,
                          int polynomial_order, double cfl);

} // namespace wasl
