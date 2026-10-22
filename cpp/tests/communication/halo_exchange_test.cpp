// Host-staged halo exchange oracle: a two-rank split reproduces the serial
// residual, proving the imported one-layer ghost trace is correct.
//
// This is a single-process test (no mpirun). The serial residual is the oracle.
// The mesh is partitioned into two ranks; each rank holds its owned elements and
// the imported one-layer ghost elements in one local array; exchange_halos_in
// _process fills the ghost coefficients; the SAME residual function runs on each
// rank-local slice. Each rank's owned-element residual must match the serial
// residual for those elements.
//
// Both 1D line and 2D triangle residuals are checked, since both read a
// neighbor's modes directly from the coefficient array at an interior face's
// right element. That read is the seam the ghost layer fills.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/basis/reference_element.hpp"
#include "wasl/communication/halo_exchange.hpp"
#include "wasl/dg/residual_operator.hpp"
#include "wasl/dg/triangle_residual.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/partition.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

// Build the exterior trace vector for a line mesh with given boundary values.
// Interior faces take the right element's trace at xi = -1; boundary faces take
// the supplied physical value.
std::vector<double> line_exterior(const wasl::LineMesh &mesh, int polynomial_order,
                                  const std::vector<double> &coefficients, double left_boundary,
                                  double right_boundary) {
  const int mode_count = static_cast<int>(wasl::line_mode_count(polynomial_order));
  std::vector<double> exterior(mesh.faces.size(), 0.0);
  std::vector<double> basis(mode_count);
  wasl::legendre_basis(-1.0, polynomial_order, basis.data());
  for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
    const auto right = mesh.faces[f].right_element;
    if (right == wasl::invalid_element) {
      exterior[f] = (mesh.faces[f].outward_normal < 0.0) ? left_boundary : right_boundary;
    } else {
      const double *modes = coefficients.data() + static_cast<std::size_t>(right) * mode_count;
      double trace = 0.0;
      for (int i = 0; i < mode_count; ++i) {
        trace += modes[i] * basis[i];
      }
      exterior[f] = trace;
    }
  }
  return exterior;
}

wasl::LineMesh make_uniform_line(double length, int n) {
  wasl::LineMesh mesh;
  const double h = length / n;
  for (int e = 0; e < n; ++e) {
    mesh.elements.push_back(wasl::LineElement{e, e + 1, e * h, (e + 1) * h, 1});
  }
  mesh.faces.push_back(wasl::LineFace{0, wasl::invalid_element, 0.0, -1.0, wasl::invalid_physical});
  for (int e = 0; e + 1 < n; ++e) {
    mesh.faces.push_back(wasl::LineFace{e, e + 1, (e + 1) * h, 1.0, wasl::invalid_physical});
  }
  mesh.faces.push_back(
      wasl::LineFace{n - 1, wasl::invalid_element, length, 1.0, wasl::invalid_physical});
  return mesh;
}

// Rank-local line mesh: owned elements followed by ghost elements, faces owned
// by this rank remapped into local indices. global_to_local maps every element
// the rank touches to its local index.
struct RankLine {
  wasl::LineMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

RankLine slice_line(const wasl::LineMesh &serial, const std::vector<std::int32_t> &owned,
                    const std::vector<std::int32_t> &ghosts) {
  RankLine out;
  std::int32_t local = 0;
  for (const auto global : owned) {
    out.mesh.elements.push_back(serial.elements.at(static_cast<std::size_t>(global)));
    out.global_to_local[global] = local++;
  }
  for (const auto global : ghosts) {
    out.mesh.elements.push_back(serial.elements.at(static_cast<std::size_t>(global)));
    out.global_to_local[global] = local++;
  }
  // Keep a face only if every element it touches is present locally (owned or
  // ghost). A partition-boundary face has the owned element on its right and the
  // imported ghost on its left, so both sides are local and it is kept; a face
  // whose other side is a second-ring neighbor (not imported) is dropped.
  for (const auto &face : serial.faces) {
    const bool left_local = out.global_to_local.count(face.left_element) != 0;
    const bool right_local = face.right_element == wasl::invalid_element ||
                             out.global_to_local.count(face.right_element);
    if (!left_local || !right_local) {
      continue;
    }
    wasl::LineFace local_face = face;
    local_face.left_element = out.global_to_local.at(face.left_element);
    if (face.right_element == wasl::invalid_element) {
      local_face.right_element = wasl::invalid_element;
    } else {
      local_face.right_element = out.global_to_local.at(face.right_element);
    }
    out.mesh.faces.push_back(local_face);
  }
  return out;
}

// Smooth P1 coefficients on a line: constant plus a sine in the cell average.
std::vector<double> line_coefficients(const wasl::LineMesh &mesh, int polynomial_order) {
  const std::size_t mode_count = wasl::line_mode_count(polynomial_order);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const double x = 0.5 * (mesh.elements[e].x_left + mesh.elements[e].x_right);
    coefficients[e * mode_count] = 1.0 + 0.4 * std::sin(x);
    if (mode_count > 1) {
      coefficients[e * mode_count + 1] = 0.2 * std::cos(x);
    }
  }
  return coefficients;
}

// The line acceptance check: split a uniform line into two ranks, exchange the
// ghost layer, run the residual per rank, and compare owned elements to serial.
bool check_line_residual_matches_serial() {
  constexpr int n = 6;
  constexpr int p = 1;
  const auto serial = make_uniform_line(1.0, n);
  const wasl::ScalarAdvection equation{1.0};
  const std::size_t mode_count = wasl::line_mode_count(p);

  // Partition and slice.
  wasl::PreparedMesh prepared = wasl::prepare_mesh(serial);
  const auto partition = wasl::partition_mesh(prepared, 2);

  std::vector<std::int32_t> owned[2];
  for (std::size_t e = 0; e < n; ++e) {
    owned[partition.element_owners[e]].push_back(static_cast<std::int32_t>(e));
  }
  const auto ghosts0 = wasl::rank_ghost_elements(prepared, partition, 0);
  const auto ghosts1 = wasl::rank_ghost_elements(prepared, partition, 1);
  const RankLine rank0 = slice_line(serial, owned[0], ghosts0);
  const RankLine rank1 = slice_line(serial, owned[1], ghosts1);

  // Serial residual as oracle.
  const auto serial_coef = line_coefficients(serial, p);
  const auto serial_ext = line_exterior(serial, p, serial_coef, 0.0, 0.0);
  std::vector<double> serial_rhs(serial_coef.size());
  wasl::scalar_advection_residual(serial, equation, p, serial_coef.data(), serial_ext.data(),
                                  serial_rhs.data());

  // Rank-local coefficient arrays (owned then ghost), seeded from serial.
  const std::size_t local_total0 = rank0.mesh.num_elements();
  const std::size_t local_total1 = rank1.mesh.num_elements();
  std::vector<double> coef0(local_total0 * mode_count), coef1(local_total1 * mode_count);
  for (const auto global : owned[0]) {
    const auto loc = rank0.global_to_local.at(global);
    std::copy_n(serial_coef.data() + global * mode_count, mode_count,
                coef0.data() + loc * mode_count);
  }
  for (const auto global : owned[1]) {
    const auto loc = rank1.global_to_local.at(global);
    std::copy_n(serial_coef.data() + global * mode_count, mode_count,
                coef1.data() + loc * mode_count);
  }

  // Exchange ghosts in process.
  std::vector<wasl::HaloSchedule> schedules = {
      wasl::build_halo_schedule(prepared, partition, 0),
      wasl::build_halo_schedule(prepared, partition, 1),
  };
  std::vector<const double *> src = {coef0.data(), coef1.data()};
  std::vector<std::map<std::int32_t, std::int32_t>> maps = {rank0.global_to_local,
                                                            rank1.global_to_local};
  wasl::exchange_halos_in_process(schedules, 2, mode_count, 0, mode_count, src, maps);

  // Per-rank residuals. True boundaries get exterior 0; partition-boundary
  // faces are interior in the local mesh and read the ghost trace directly.
  const auto ext0 = line_exterior(rank0.mesh, p, coef0, 0.0, 0.0);
  const auto ext1 = line_exterior(rank1.mesh, p, coef1, 0.0, 0.0);
  std::vector<double> rhs0(coef0.size()), rhs1(coef1.size());
  wasl::scalar_advection_residual(rank0.mesh, equation, p, coef0.data(), ext0.data(), rhs0.data());
  wasl::scalar_advection_residual(rank1.mesh, equation, p, coef1.data(), ext1.data(), rhs1.data());

  // Compare each owned element's residual to the serial residual.
  for (const auto global : owned[0]) {
    const auto loc = rank0.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(rhs0[loc * mode_count + m], serial_rhs[global * mode_count + m]))
        return false;
    }
  }
  for (const auto global : owned[1]) {
    const auto loc = rank1.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(rhs1[loc * mode_count + m], serial_rhs[global * mode_count + m]))
        return false;
    }
  }
  return true;
}

// ---- Triangle path ----

wasl::TriangleMesh make_unit_square_triangles(int cells_per_axis) {
  wasl::MshDocument document;
  document.surface_entities.push_back(
      wasl::MshSurfaceEntity{1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, {1}, {}});
  for (int j = 0; j <= cells_per_axis; ++j) {
    for (int i = 0; i <= cells_per_axis; ++i) {
      const auto tag = static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1);
      document.nodes.push_back(wasl::MshNode{tag, 2, 1, static_cast<double>(i) / cells_per_axis,
                                             static_cast<double>(j) / cells_per_axis, 0.0});
    }
  }
  auto node = [cells_per_axis](int i, int j) {
    return static_cast<std::int64_t>(j * (cells_per_axis + 1) + i + 1);
  };
  std::int64_t element_tag = 1;
  for (int j = 0; j < cells_per_axis; ++j) {
    for (int i = 0; i < cells_per_axis; ++i) {
      const auto lower_left = node(i, j);
      const auto lower_right = node(i + 1, j);
      const auto upper_right = node(i + 1, j + 1);
      const auto upper_left = node(i, j + 1);
      document.elements.push_back(
          wasl::MshElement{element_tag++, 2, 1, 2, {lower_left, lower_right, upper_right}});
      document.elements.push_back(
          wasl::MshElement{element_tag++, 2, 1, 2, {lower_left, upper_right, upper_left}});
    }
  }
  return wasl::build_triangle_mesh(document);
}

struct RankTriangle {
  wasl::TriangleMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

RankTriangle slice_triangle(const wasl::TriangleMesh &serial,
                            const std::vector<std::int32_t> &owned,
                            const std::vector<std::int32_t> &ghosts) {
  RankTriangle out;
  std::int32_t local = 0;
  for (const auto global : owned) {
    out.mesh.elements.push_back(serial.elements.at(static_cast<std::size_t>(global)));
    out.global_to_local[global] = local++;
  }
  for (const auto global : ghosts) {
    out.mesh.elements.push_back(serial.elements.at(static_cast<std::size_t>(global)));
    out.global_to_local[global] = local++;
  }
  for (const auto &face : serial.faces) {
    const bool left_local = out.global_to_local.count(face.left_element) != 0;
    const bool right_local = face.right_element == wasl::invalid_element ||
                             out.global_to_local.count(face.right_element);
    if (!left_local || !right_local) {
      continue;
    }
    wasl::TriangleFace local_face = face;
    local_face.left_element = out.global_to_local.at(face.left_element);
    if (face.right_element == wasl::invalid_element) {
      local_face.right_element = wasl::invalid_element;
    } else {
      local_face.right_element = out.global_to_local.at(face.right_element);
    }
    out.mesh.faces.push_back(local_face);
  }
  return out;
}

std::vector<double> triangle_coefficients(const wasl::TriangleMesh &mesh, int polynomial_order) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const auto &elem = mesh.elements[e];
    const double x = (elem.x0 + elem.x1 + elem.x2) / 3.0;
    const double y = (elem.y0 + elem.y1 + elem.y2) / 3.0;
    coefficients[e * mode_count] = (1.0 + 0.3 * std::sin(x) + 0.2 * std::cos(y)) * std::sqrt(2.0);
    if (mode_count >= 3) {
      coefficients[e * mode_count + 1] = 0.15 * x;
      coefficients[e * mode_count + 2] = 0.15 * y;
    }
  }
  return coefficients;
}

std::vector<double> boundary_zero(const wasl::TriangleMesh &mesh, int polynomial_order) {
  return std::vector<double>(mesh.faces.size() * static_cast<std::size_t>(polynomial_order + 1),
                             0.0);
}

bool check_triangle_residual_matches_serial() {
  constexpr int cells = 2; // 8 triangles, several interior faces per rank
  constexpr int p = 1;
  const auto serial = make_unit_square_triangles(cells);
  const wasl::ScalarAdvection2D equation{0.7, -1.1};
  const std::size_t mode_count = wasl::triangle_mode_count(p);

  wasl::PreparedMesh prepared = wasl::prepare_mesh(serial);
  const auto partition = wasl::partition_mesh(prepared, 2);

  std::vector<std::int32_t> owned[2];
  for (std::size_t e = 0; e < serial.num_elements(); ++e) {
    owned[partition.element_owners[e]].push_back(static_cast<std::int32_t>(e));
  }
  const auto ghosts0 = wasl::rank_ghost_elements(prepared, partition, 0);
  const auto ghosts1 = wasl::rank_ghost_elements(prepared, partition, 1);
  const RankTriangle rank0 = slice_triangle(serial, owned[0], ghosts0);
  const RankTriangle rank1 = slice_triangle(serial, owned[1], ghosts1);

  const auto serial_coef = triangle_coefficients(serial, p);
  const auto serial_ext = boundary_zero(serial, p);
  std::vector<double> serial_rhs(serial_coef.size());
  wasl::triangle_scalar_advection_residual(serial, equation, p, serial_coef.data(),
                                           serial_ext.data(), serial_rhs.data());

  std::vector<double> coef0(rank0.mesh.num_elements() * mode_count);
  std::vector<double> coef1(rank1.mesh.num_elements() * mode_count);
  for (const auto global : owned[0]) {
    const auto loc = rank0.global_to_local.at(global);
    std::copy_n(serial_coef.data() + global * mode_count, mode_count,
                coef0.data() + loc * mode_count);
  }
  for (const auto global : owned[1]) {
    const auto loc = rank1.global_to_local.at(global);
    std::copy_n(serial_coef.data() + global * mode_count, mode_count,
                coef1.data() + loc * mode_count);
  }

  std::vector<wasl::HaloSchedule> schedules = {
      wasl::build_halo_schedule(prepared, partition, 0),
      wasl::build_halo_schedule(prepared, partition, 1),
  };
  std::vector<const double *> src = {coef0.data(), coef1.data()};
  std::vector<std::map<std::int32_t, std::int32_t>> maps = {rank0.global_to_local,
                                                            rank1.global_to_local};
  wasl::exchange_halos_in_process(schedules, 2, mode_count, 0, mode_count, src, maps);

  const auto ext0 = boundary_zero(rank0.mesh, p);
  const auto ext1 = boundary_zero(rank1.mesh, p);
  std::vector<double> rhs0(coef0.size()), rhs1(coef1.size());
  wasl::triangle_scalar_advection_residual(rank0.mesh, equation, p, coef0.data(), ext0.data(),
                                           rhs0.data());
  wasl::triangle_scalar_advection_residual(rank1.mesh, equation, p, coef1.data(), ext1.data(),
                                           rhs1.data());

  for (const auto global : owned[0]) {
    const auto loc = rank0.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(rhs0[loc * mode_count + m], serial_rhs[global * mode_count + m]))
        return false;
    }
  }
  for (const auto global : owned[1]) {
    const auto loc = rank1.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(rhs1[loc * mode_count + m], serial_rhs[global * mode_count + m]))
        return false;
    }
  }
  return true;
}

// The exchange must import exactly the owner's current coefficients into each
// ghost slot. Seed owned coefficients with a distinct value per element, leave
// ghosts zero, exchange, and assert every ghost now equals its owner's value.
// This isolates the halo machinery from the residual.
bool check_ghost_coefficients_imported() {
  constexpr int n = 6;
  constexpr int p = 1;
  const auto serial = make_uniform_line(1.0, n);
  const std::size_t mode_count = wasl::line_mode_count(p);

  wasl::PreparedMesh prepared = wasl::prepare_mesh(serial);
  const auto partition = wasl::partition_mesh(prepared, 2);
  std::vector<std::int32_t> owned[2];
  for (std::size_t e = 0; e < n; ++e) {
    owned[partition.element_owners[e]].push_back(static_cast<std::int32_t>(e));
  }
  const RankLine rank0 =
      slice_line(serial, owned[0], wasl::rank_ghost_elements(prepared, partition, 0));
  const RankLine rank1 =
      slice_line(serial, owned[1], wasl::rank_ghost_elements(prepared, partition, 1));

  // Owned modes carry a unique marker (global index); ghosts start at zero.
  std::vector<double> coef0(rank0.mesh.num_elements() * mode_count, 0.0);
  std::vector<double> coef1(rank1.mesh.num_elements() * mode_count, 0.0);
  for (const auto global : owned[0]) {
    const auto loc = rank0.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m)
      coef0[loc * mode_count + m] = static_cast<double>(global) + 0.1 * static_cast<double>(m);
  }
  for (const auto global : owned[1]) {
    const auto loc = rank1.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m)
      coef1[loc * mode_count + m] = static_cast<double>(global) + 0.1 * static_cast<double>(m);
  }

  std::vector<wasl::HaloSchedule> schedules = {
      wasl::build_halo_schedule(prepared, partition, 0),
      wasl::build_halo_schedule(prepared, partition, 1),
  };
  std::vector<const double *> src = {coef0.data(), coef1.data()};
  std::vector<std::map<std::int32_t, std::int32_t>> maps = {rank0.global_to_local,
                                                            rank1.global_to_local};
  wasl::exchange_halos_in_process(schedules, 2, mode_count, 0, mode_count, src, maps);

  // Rank 0 imports ghost 3; rank 1 imports ghost 2. Each must equal the owner.
  for (const auto global : wasl::rank_ghost_elements(prepared, partition, 0)) {
    const auto loc = rank0.global_to_local.at(global);
    const auto owner = partition.element_owners.at(static_cast<std::size_t>(global));
    const auto &owner_coef = owner == 0 ? coef0 : coef1;
    const auto &owner_map = owner == 0 ? rank0.global_to_local : rank1.global_to_local;
    const auto owner_loc = owner_map.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(coef0[loc * mode_count + m], owner_coef[owner_loc * mode_count + m]))
        return false;
    }
  }
  for (const auto global : wasl::rank_ghost_elements(prepared, partition, 1)) {
    const auto loc = rank1.global_to_local.at(global);
    const auto owner = partition.element_owners.at(static_cast<std::size_t>(global));
    const auto &owner_coef = owner == 0 ? coef0 : coef1;
    const auto &owner_map = owner == 0 ? rank0.global_to_local : rank1.global_to_local;
    const auto owner_loc = owner_map.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(coef1[loc * mode_count + m], owner_coef[owner_loc * mode_count + m]))
        return false;
    }
  }
  return true;
}

} // namespace

int main() {
  try {
    require(check_line_residual_matches_serial(),
            "line two-rank residual does not match serial oracle");
    require(check_triangle_residual_matches_serial(),
            "triangle two-rank residual does not match serial oracle");
    require(check_ghost_coefficients_imported(),
            "halo exchange did not import owner coefficients into ghost slots");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "exception: " << error.what() << '\n';
    return 1;
  }
}
