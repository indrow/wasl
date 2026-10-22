// Periodic exterior traces for triangle meshes of an axis-aligned box.
//
// Every boundary edge quadrature point is translated by one box period across
// the domain and located in its donor element once, during setup. Stage
// evaluation then only contracts donor coefficients with stored basis values.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// Donor element and basis values per [face][edge quadrature point]; interior
// faces keep invalid_element donors.
struct PeriodicTriangleTrace {
  std::vector<std::int32_t> elements;
  std::vector<double> basis;
  int point_count = 0;
};

// Throws when the boundary is not an axis-aligned box or a donor is missing.
PeriodicTriangleTrace periodic_triangle_trace(const TriangleMesh &mesh, int polynomial_order);

// Scalar exterior per [face][point] from solver coefficients U[element][mode].
std::vector<double> periodic_triangle_exterior(const PeriodicTriangleTrace &trace,
                                               const std::vector<double> &coefficients,
                                               std::size_t mode_count);

// Euler exterior per [variable][face][point] from U[variable][mode][element].
std::vector<double> periodic_triangle_euler_exterior(const PeriodicTriangleTrace &trace,
                                                     const std::vector<double> &coefficients,
                                                     std::size_t element_count,
                                                     std::size_t mode_count,
                                                     std::size_t face_count);

} // namespace wasl
