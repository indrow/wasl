// Case-level outflow boundaries for 1D Euler.
//
// On an irregular unit interval whose ends are the physical points "left" and
// "right":
//   - a uniform moving state with outflow ends is preserved to 1e-12, with and
//     without limiting;
//   - in the Sod shock tube both ends are at rest, so while the waves are far
//     from them the mass flux through the outflow ends is zero and total
//     density mass is conserved to round-off at every step;
//   - Euler inflow is rejected naming the boundary.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/case/euler_case.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/prepared_mesh_file.hpp"
#include "wasl/mesh/prepared/prepare_mesh.hpp"

#include "named_line_mesh.hpp"

namespace fs = std::filesystem;

namespace {

// Round-off for O(1) states over a few hundred steps.
constexpr double kTolerance = 1e-12;

wasl::NamedBoundaries outflow_ends() {
  const wasl::BoundaryCondition outflow{wasl::BoundaryKind::outflow, 0.0};
  return {{"left", outflow}, {"right", outflow}};
}

wasl::EulerCase base_case(const fs::path &mesh, const fs::path &output, bool limiting) {
  wasl::EulerCase config;
  config.files = {mesh, output, {}};
  config.run = {0.0, 0.2, 100000, 1, 100000, 100000, 1};
  config.limiting = limiting;
  config.polynomial_order = 2;
  config.cfl = 0.2;
  config.boundaries = outflow_ends();
  return config;
}

double free_stream_deviation(const fs::path &mesh, const fs::path &root, bool limiting) {
  const wasl::Euler equation{1.4};
  const wasl::EulerPrimitive state{1.0, 0.5, 0.0, 1.0};
  auto config = base_case(mesh, root / (limiting ? "free-limited" : "free"), limiting);
  config.initial_regions = {{wasl::EulerRegion::Shape::default_region, 0, 0, 0, 0, 0, state}};
  const auto checkpoint = wasl::read_checkpoint(wasl::run_euler_case(config).checkpoint);
  const auto conserved = wasl::conserved(equation, state);
  double deviation = 0.0;
  for (std::size_t variable = 0; variable < checkpoint.variable_count; ++variable) {
    for (std::size_t mode = 0; mode < checkpoint.mode_count; ++mode) {
      for (std::size_t element = 0; element < checkpoint.element_count; ++element) {
        const double expected =
            mode == 0
                ? wasl::euler_component(conserved, static_cast<int>(variable)) * std::sqrt(2.0)
                : 0.0;
        const double actual =
            checkpoint
                .coefficients[(variable * checkpoint.mode_count + mode) * checkpoint.element_count +
                              element];
        deviation = std::max(deviation, std::fabs(actual - expected));
      }
    }
  }
  return deviation;
}

// The DG precursors ahead of the shock and the rarefaction head decay a fixed
// factor per element. At t = 0.05 the waves are still 0.4 of the domain from
// the ends and the drift is 4e-15; it passes 1e-12 near t = 0.11 and reaches
// 1e-8 by t = 0.2 on this 40-element mesh.
constexpr double kSodConservationEnd = 0.05;

double sod_mass_drift(const fs::path &mesh, const fs::path &root) {
  auto config = base_case(mesh, root / "sod", true);
  config.run.end_time = kSodConservationEnd;
  config.initial_regions = {
      {wasl::EulerRegion::Shape::default_region, 0, 0, 0, 0, 0, {1.0, 0.0, 0.0, 1.0}},
      {wasl::EulerRegion::Shape::interval, 0.5, 1.0, 0, 0, 0, {0.125, 0.0, 0.0, 0.1}}};
  wasl::run_euler_case(config);
  std::ifstream history(config.files.output_directory / "history.csv");
  std::string line;
  std::getline(history, line); // header
  double initial = std::nan("");
  double drift = 0.0;
  while (std::getline(history, line)) {
    const double mass = std::stod(line.substr(line.rfind(',') + 1));
    if (std::isnan(initial))
      initial = mass;
    drift = std::max(drift, std::fabs(mass - initial));
  }
  return drift;
}

bool inflow_is_rejected(const fs::path &mesh, const fs::path &root) {
  auto config = base_case(mesh, root / "inflow", false);
  config.initial_regions = {
      {wasl::EulerRegion::Shape::default_region, 0, 0, 0, 0, 0, {1.0, 0.0, 0.0, 1.0}}};
  config.boundaries[0].condition = {wasl::BoundaryKind::inflow, 1.0};
  try {
    wasl::run_euler_case(config);
  } catch (const std::exception &error) {
    return std::string(error.what()).find("boundary left: Euler inflow") != std::string::npos;
  }
  return false;
}

} // namespace

int main() {
  const auto root = fs::temp_directory_path() / "wasl-euler-boundary-case-test";
  fs::remove_all(root);
  fs::create_directories(root);
  const auto mesh = root / "line.h5";
  wasl::write_prepared_mesh(mesh, wasl::prepare_mesh(wasl_test::named_line(40)));

  const double free = free_stream_deviation(mesh, root, false);
  const double free_limited = free_stream_deviation(mesh, root, true);
  const double sod = sod_mass_drift(mesh, root);
  const bool rejected = inflow_is_rejected(mesh, root);
  std::printf("free stream %.3e, limited %.3e, Sod mass drift %.3e, inflow rejected %d\n", free,
              free_limited, sod, rejected);
  fs::remove_all(root);
  return free <= kTolerance && free_limited <= kTolerance && sod <= kTolerance && rejected ? 0 : 1;
}
