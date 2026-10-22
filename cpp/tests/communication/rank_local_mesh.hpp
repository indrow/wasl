// Shared rank-local mesh slicing and seeding helpers for halo-exchange tests.
//
// A rank-local mesh holds owned elements followed by imported ghost elements,
// with faces kept only when every element they touch is present locally. The
// existing residual reads a neighbor's modes directly from the coefficient
// array, so this layout lets the residual run unchanged on a rank-local slice.
//
// Used by both the in-process halo oracle (tests/communication/halo_exchange
// _test.cpp) and the MPI transport test (halo_mpi_test.cpp) so they build
// identical rank-local state.

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "wasl/basis/reference_element.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/msh_reader.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl_testing {

// Uniform line mesh on [0, length] with n elements, canonical faces.
inline wasl::LineMesh make_uniform_line(double length, int n) {
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

inline wasl::TriangleMesh make_unit_square_triangles(int cells_per_axis) {
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

struct RankLine {
  wasl::LineMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

struct RankTriangle {
  wasl::TriangleMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

inline RankLine slice_line(const wasl::LineMesh &serial, const std::vector<std::int32_t> &owned,
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
  // Keep a face only if every element it touches is present locally. A
  // partition-boundary face has the owned element on one side and the imported
  // ghost on the other, so both are local; a second-ring face is dropped.
  for (const auto &face : serial.faces) {
    const bool left_local = out.global_to_local.count(face.left_element) != 0;
    const bool right_local = face.right_element == wasl::invalid_element ||
                             out.global_to_local.count(face.right_element);
    if (!left_local || !right_local) {
      continue;
    }
    wasl::LineFace local_face = face;
    local_face.left_element = out.global_to_local.at(face.left_element);
    local_face.right_element = face.right_element == wasl::invalid_element
                                   ? wasl::invalid_element
                                   : out.global_to_local.at(face.right_element);
    out.mesh.faces.push_back(local_face);
  }
  return out;
}

inline RankTriangle slice_triangle(const wasl::TriangleMesh &serial,
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
    local_face.right_element = face.right_element == wasl::invalid_element
                                   ? wasl::invalid_element
                                   : out.global_to_local.at(face.right_element);
    out.mesh.faces.push_back(local_face);
  }
  return out;
}

} // namespace wasl_testing
