// Immediate-neighbor characteristic WENO plus positivity limiting for Euler.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/limiting/positivity/euler_positivity.hpp"
#include "wasl/limiting/weno/candidate_transform.hpp"
#include "wasl/limiting/weno/euler_characteristic.hpp"
#include "wasl/mesh/line_geometry.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/quadrature.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace wasl {

struct EulerLimiterData {
  std::size_t element_count = 0;
  std::size_t mode_count = 0;
  int polynomial_order = 0;
  int maximum_neighbors = 0;
  std::size_t highest_degree_begin = 0;
  std::vector<std::int32_t> neighbors;
  std::vector<double> transforms;
  std::vector<double> normal_x, normal_y;
  std::vector<double> check_basis;
};

inline void append_unique_neighbor(EulerLimiterData &data, std::size_t element,
                                   std::int32_t neighbor) {
  auto *values = data.neighbors.data() + element * data.maximum_neighbors;
  for (int index = 0; index < data.maximum_neighbors; ++index) {
    if (values[index] == neighbor)
      return;
    if (values[index] == invalid_element) {
      values[index] = neighbor;
      return;
    }
  }
}

inline EulerLimiterData make_euler_limiter_data(const LineMesh &mesh, int polynomial_order) {
  EulerLimiterData data;
  data.element_count = mesh.num_elements();
  data.mode_count = line_mode_count(polynomial_order);
  data.polynomial_order = polynomial_order;
  data.maximum_neighbors = 2;
  data.highest_degree_begin = polynomial_order > 0 ? data.mode_count - 1 : data.mode_count;
  data.neighbors.assign(data.element_count * data.maximum_neighbors, invalid_element);
  data.normal_x.assign(data.element_count, 1.0);
  data.normal_y.assign(data.element_count, 0.0);
  std::int32_t left_boundary = invalid_element;
  std::int32_t right_boundary = invalid_element;
  for (const auto &face : mesh.faces) {
    if (face.right_element != invalid_element) {
      append_unique_neighbor(data, face.left_element, face.right_element);
      append_unique_neighbor(data, face.right_element, face.left_element);
    } else if (face.outward_normal < 0.0) {
      left_boundary = face.left_element;
    } else {
      right_boundary = face.left_element;
    }
  }
  if (left_boundary != invalid_element && right_boundary != invalid_element) {
    append_unique_neighbor(data, left_boundary, right_boundary);
    append_unique_neighbor(data, right_boundary, left_boundary);
  }
  const auto matrix_size = data.mode_count * data.mode_count;
  data.transforms.assign(data.element_count * data.maximum_neighbors * matrix_size, 0.0);
  const auto rule = gauss_legendre(std::max(1, polynomial_order + 1));
  // The wrapped neighbor sits one period away; shift it next to the target so
  // the transform evaluates its polynomial inside its own element.
  const double period = mesh.elements.back().x_right - mesh.elements.front().x_left;
  std::vector<double> source_basis(data.mode_count), target_basis(data.mode_count);
  for (std::size_t element = 0; element < data.element_count; ++element) {
    const auto target = element_map(mesh, element);
    for (int candidate = 0; candidate < data.maximum_neighbors; ++candidate) {
      const auto neighbor = data.neighbors[element * data.maximum_neighbors + candidate];
      if (neighbor == invalid_element)
        continue;
      auto source = element_map(mesh, static_cast<std::size_t>(neighbor));
      const double gap = target.centre() - source.centre();
      const double shift = gap > 0.5 * period ? period : gap < -0.5 * period ? -period : 0.0;
      source.x_left += shift;
      source.x_right += shift;
      auto *matrix =
          data.transforms.data() + (element * data.maximum_neighbors + candidate) * matrix_size;
      for (int point = 0; point < rule.size(); ++point) {
        const auto x = target.reference_to_physical(rule.points[point]);
        const auto source_xi = source.physical_to_reference(x);
        legendre_basis(source_xi, polynomial_order, source_basis.data());
        legendre_basis(rule.points[point], polynomial_order, target_basis.data());
        for (std::size_t target_mode = 0; target_mode < data.mode_count; ++target_mode)
          for (std::size_t source_mode = 0; source_mode < data.mode_count; ++source_mode)
            matrix[target_mode * data.mode_count + source_mode] +=
                rule.weights[point] * target_basis[target_mode] * source_basis[source_mode];
      }
    }
  }
  const auto append_basis = [&](double xi) {
    const auto offset = data.check_basis.size();
    data.check_basis.resize(offset + data.mode_count);
    legendre_basis(xi, polynomial_order, data.check_basis.data() + offset);
  };
  for (const auto point : rule.points)
    append_basis(point);
  append_basis(-1.0);
  append_basis(1.0);
  return data;
}

inline EulerLimiterData make_euler_limiter_data(const TriangleMesh &mesh, int polynomial_order,
                                                const std::vector<std::int32_t> &periodic_donors,
                                                int edge_point_count) {
  EulerLimiterData data;
  data.element_count = mesh.num_elements();
  data.mode_count = triangle_mode_count(polynomial_order);
  data.polynomial_order = polynomial_order;
  data.maximum_neighbors = 3;
  data.highest_degree_begin =
      polynomial_order > 0 ? triangle_mode_count(polynomial_order - 1) : data.mode_count;
  data.neighbors.assign(data.element_count * data.maximum_neighbors, invalid_element);
  data.normal_x.assign(data.element_count, 1.0);
  data.normal_y.assign(data.element_count, 0.0);
  std::vector<unsigned char> has_normal(data.element_count, 0);
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    const auto left = static_cast<std::size_t>(face.left_element);
    auto right = face.right_element;
    if (right == invalid_element && !periodic_donors.empty())
      right = periodic_donors[face_index * edge_point_count];
    if (right == invalid_element)
      continue;
    append_unique_neighbor(data, left, right);
    if (!has_normal[left]) {
      data.normal_x[left] = face.outward_normal_x;
      data.normal_y[left] = face.outward_normal_y;
      has_normal[left] = 1;
    }
    if (face.right_element != invalid_element) {
      const auto right_element = static_cast<std::size_t>(right);
      append_unique_neighbor(data, right_element, face.left_element);
      if (!has_normal[right_element]) {
        data.normal_x[right_element] = -face.outward_normal_x;
        data.normal_y[right_element] = -face.outward_normal_y;
        has_normal[right_element] = 1;
      }
    }
  }
  const auto matrix_size = data.mode_count * data.mode_count;
  data.transforms.assign(data.element_count * data.maximum_neighbors * matrix_size, 0.0);
  // Periodic donors sit one box period away; shift them next to the target so
  // the projection evaluates the donor polynomial inside its own element.
  double min_x = mesh.elements.front().x0, max_x = min_x;
  double min_y = mesh.elements.front().y0, max_y = min_y;
  for (const auto &element : mesh.elements) {
    min_x = std::min({min_x, element.x0, element.x1, element.x2});
    max_x = std::max({max_x, element.x0, element.x1, element.x2});
    min_y = std::min({min_y, element.y0, element.y1, element.y2});
    max_y = std::max({max_y, element.y0, element.y1, element.y2});
  }
  const auto periodic_shift = [](double gap, double period) {
    return gap > 0.5 * period ? period : gap < -0.5 * period ? -period : 0.0;
  };
  std::vector<double> source(data.mode_count, 0.0), projected(data.mode_count);
  for (std::size_t element = 0; element < data.element_count; ++element) {
    const auto target_map = triangle_element_map(mesh, element);
    for (int candidate = 0; candidate < data.maximum_neighbors; ++candidate) {
      const auto neighbor = data.neighbors[element * data.maximum_neighbors + candidate];
      if (neighbor == invalid_element)
        continue;
      auto source_map = triangle_element_map(mesh, static_cast<std::size_t>(neighbor));
      const double shift_x = periodic_shift((target_map.x0 + target_map.x1 + target_map.x2 -
                                             source_map.x0 - source_map.x1 - source_map.x2) /
                                                3.0,
                                            max_x - min_x);
      const double shift_y = periodic_shift((target_map.y0 + target_map.y1 + target_map.y2 -
                                             source_map.y0 - source_map.y1 - source_map.y2) /
                                                3.0,
                                            max_y - min_y);
      source_map.x0 += shift_x;
      source_map.x1 += shift_x;
      source_map.x2 += shift_x;
      source_map.y0 += shift_y;
      source_map.y1 += shift_y;
      source_map.y2 += shift_y;
      auto *matrix =
          data.transforms.data() + (element * data.maximum_neighbors + candidate) * matrix_size;
      for (std::size_t source_mode = 0; source_mode < data.mode_count; ++source_mode) {
        source[source_mode] = 1.0;
        project_neighbor_polynomial(source_map, target_map, polynomial_order, source.data(),
                                    projected.data());
        source[source_mode] = 0.0;
        for (std::size_t target_mode = 0; target_mode < data.mode_count; ++target_mode)
          matrix[target_mode * data.mode_count + source_mode] = projected[target_mode];
      }
    }
  }
  const auto volume = triangle_volume_rule(polynomial_order);
  const auto edges = triangle_edge_rules(polynomial_order);
  const auto append_basis = [&](double xi, double eta) {
    const auto offset = data.check_basis.size();
    data.check_basis.resize(offset + data.mode_count);
    dubiner_basis(xi, eta, polynomial_order, data.check_basis.data() + offset);
  };
  for (int point = 0; point < volume.size(); ++point)
    append_basis(volume.xi[point], volume.eta[point]);
  for (const auto &edge : edges)
    for (std::size_t point = 0; point < edge.xi.size(); ++point)
      append_basis(edge.xi[point], edge.eta[point]);
  return data;
}

inline void limit_euler_modal(const EulerLimiterData &data, const Euler &equation, double threshold,
                              double density_floor, double pressure_floor,
                              std::vector<double> &coefficients) {
  if (data.polynomial_order == 0)
    return;
  const auto source = coefficients;
  const auto matrix_size = data.mode_count * data.mode_count;
  std::vector<double> target(Euler::variable_count * data.mode_count);
  std::vector<double> candidates(data.maximum_neighbors * Euler::variable_count * data.mode_count);
  std::vector<double> limited(target.size());
  for (std::size_t element = 0; element < data.element_count; ++element) {
    double total_energy = 0.0;
    double high_energy = 0.0;
    for (std::size_t mode = 0; mode < data.mode_count; ++mode) {
      const auto value = source[mode * data.element_count + element];
      total_energy += value * value;
      if (mode >= data.highest_degree_begin)
        high_energy += value * value;
    }
    if (total_energy == 0.0 || high_energy / total_energy <= threshold)
      continue;
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < data.mode_count; ++mode)
        target[static_cast<std::size_t>(variable) * data.mode_count + mode] =
            source[(static_cast<std::size_t>(variable) * data.mode_count + mode) *
                       data.element_count +
                   element];
    std::size_t candidate_count = 0;
    for (int candidate = 0; candidate < data.maximum_neighbors; ++candidate) {
      const auto neighbor = data.neighbors[element * data.maximum_neighbors + candidate];
      if (neighbor == invalid_element)
        continue;
      const auto *matrix =
          data.transforms.data() + (element * data.maximum_neighbors + candidate) * matrix_size;
      for (int variable = 0; variable < Euler::variable_count; ++variable) {
        auto *result = candidates.data() +
                       (candidate_count * Euler::variable_count + variable) * data.mode_count;
        std::fill(result, result + data.mode_count, 0.0);
        for (std::size_t target_mode = 0; target_mode < data.mode_count; ++target_mode)
          for (std::size_t source_mode = 0; source_mode < data.mode_count; ++source_mode)
            result[target_mode] +=
                matrix[target_mode * data.mode_count + source_mode] *
                source[(static_cast<std::size_t>(variable) * data.mode_count + source_mode) *
                           data.element_count +
                       static_cast<std::size_t>(neighbor)];
        result[0] = target[static_cast<std::size_t>(variable) * data.mode_count];
      }
      ++candidate_count;
    }
    const EulerConserved average{
        target[0] / std::sqrt(2.0), target[data.mode_count] / std::sqrt(2.0),
        target[2 * data.mode_count] / std::sqrt(2.0), target[3 * data.mode_count] / std::sqrt(2.0)};
    euler_characteristic_weno(equation, average, data.normal_x[element], data.normal_y[element],
                              target.data(), candidates.data(), candidate_count, data.mode_count,
                              limited.data());
    limit_euler_positivity(equation, limited.data(), data.mode_count, data.check_basis.data(),
                           data.check_basis.size() / data.mode_count, density_floor,
                           pressure_floor);
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < data.mode_count; ++mode)
        coefficients[(static_cast<std::size_t>(variable) * data.mode_count + mode) *
                         data.element_count +
                     element] =
            limited[static_cast<std::size_t>(variable) * data.mode_count + mode];
  }
}

} // namespace wasl
