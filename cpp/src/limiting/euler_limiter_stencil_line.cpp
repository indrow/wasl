#include "wasl/limiting/euler_limiter_stencil.hpp"

#include <algorithm>
#include <vector>

#include "limiting/limiter_stencil_assembly.hpp"
#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/quadrature/gauss_legendre.hpp"

namespace wasl {
namespace {

void connect_neighbors(EulerLimiterStencil &stencil, const LineMesh &mesh, LineEnds ends) {
  std::int32_t left_boundary = invalid_element;
  std::int32_t right_boundary = invalid_element;
  for (const auto &face : mesh.faces) {
    if (!face.is_boundary()) {
      append_unique_neighbor(stencil, static_cast<std::size_t>(face.left_element),
                             face.right_element);
      append_unique_neighbor(stencil, static_cast<std::size_t>(face.right_element),
                             face.left_element);
    } else if (face.outward_normal < 0.0) {
      left_boundary = face.left_element;
    } else {
      right_boundary = face.left_element;
    }
  }
  if (ends == LineEnds::periodic && left_boundary != invalid_element &&
      right_boundary != invalid_element) {
    append_unique_neighbor(stencil, static_cast<std::size_t>(left_boundary), right_boundary);
    append_unique_neighbor(stencil, static_cast<std::size_t>(right_boundary), left_boundary);
  }
}

// L2 projection of each neighbor basis function onto the target basis.
void build_transforms(EulerLimiterStencil &stencil, const LineMesh &mesh,
                      const QuadratureRule &rule) {
  allocate_transforms(stencil);
  const double period = mesh.elements.back().x_right - mesh.elements.front().x_left;
  std::vector<double> source_basis(stencil.mode_count), target_basis(stencil.mode_count);
  for (std::size_t element = 0; element < stencil.element_count; ++element) {
    const auto target = element_map(mesh, element);
    for (int candidate = 0; candidate < stencil.maximum_neighbors; ++candidate) {
      const auto neighbor = stencil.neighbor(element, candidate);
      if (neighbor == invalid_element)
        continue;
      auto source = element_map(mesh, static_cast<std::size_t>(neighbor));
      source.shift(periodic_offset(target.centre() - source.centre(), period));
      auto *matrix = stencil.transform(element, candidate);
      for (int point = 0; point < rule.size(); ++point) {
        const auto source_xi =
            source.physical_to_reference(target.reference_to_physical(rule.points[point]));
        legendre_basis(source_xi, stencil.polynomial_order, source_basis.data());
        legendre_basis(rule.points[point], stencil.polynomial_order, target_basis.data());
        for (std::size_t target_mode = 0; target_mode < stencil.mode_count; ++target_mode)
          for (std::size_t source_mode = 0; source_mode < stencil.mode_count; ++source_mode)
            matrix[target_mode * stencil.mode_count + source_mode] +=
                rule.weights[point] * target_basis[target_mode] * source_basis[source_mode];
      }
    }
  }
}

void build_check_basis(EulerLimiterStencil &stencil, const QuadratureRule &rule) {
  const auto append = [&](double xi) {
    const auto offset = stencil.check_basis.size();
    stencil.check_basis.resize(offset + stencil.mode_count);
    legendre_basis(xi, stencil.polynomial_order, stencil.check_basis.data() + offset);
  };
  for (const auto point : rule.points)
    append(point);
  append(-1.0);
  append(1.0);
}

} // namespace

EulerLimiterStencil make_euler_limiter_stencil(const LineMesh &mesh, int polynomial_order,
                                               LineEnds ends) {
  EulerLimiterStencil stencil;
  stencil.element_count = mesh.num_elements();
  stencil.mode_count = line_mode_count(polynomial_order);
  stencil.polynomial_order = polynomial_order;
  stencil.maximum_neighbors = 2;
  stencil.highest_degree_begin = polynomial_order > 0 ? stencil.mode_count - 1 : stencil.mode_count;
  allocate_neighbors(stencil);
  connect_neighbors(stencil, mesh, ends);

  const auto rule = gauss_legendre(std::max(1, polynomial_order + 1));
  build_transforms(stencil, mesh, rule);
  build_check_basis(stencil, rule);
  return stencil;
}

} // namespace wasl
