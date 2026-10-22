// Case-level inflow and outflow boundaries for 1D scalar advection.
//
// On an irregular unit-interval mesh whose ends are the physical points "left"
// and "right":
//   - a free stream equal to the inflow value is preserved to 1e-12;
//   - with a zero background and an inflow value c, every history mass equals
//     the exact inflow balance |a| c t to round-off. This holds only if the
//     boundary value enters the flux at every SSP-RK3 stage. It is checked for
//     wind in both directions, so each end acts once as the inflow;
//   - misconfigured boundaries fail with a message naming the boundary.

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/case/scalar_advection_case.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/prepared_mesh_file.hpp"
#include "wasl/mesh/prepared/prepare_mesh.hpp"

#include "named_line_mesh.hpp"

namespace fs = std::filesystem;

namespace {

// Round-off over a few hundred steps of O(1) values.
constexpr double kTolerance = 1e-12;

// The upwind DG solution carries a small precursor ahead of the front that
// decays about tenfold per half element. The balance is exact only while the
// precursor at the outflow end is below round-off: ending at t = 0.2 keeps the
// front 0.8 of the domain away (deviation 4e-15 here), whereas at t = 0.4 the
// precursor has reached the outflow and the deviation is 1e-10.
constexpr double kEndTime = 0.2;

wasl::ScalarAdvectionCase base_case(const fs::path &mesh, const fs::path &output) {
  wasl::ScalarAdvectionCase config;
  config.files = {mesh, output, {}};
  config.run = {0.0, kEndTime, 100000, 1, 100000, 100000, 1};
  config.polynomial_order = 2;
  config.wavespeed = 1.0;
  config.cfl = 0.2;
  return config;
}

wasl::NamedBoundary inflow(const std::string &name, double value) {
  return {name, {wasl::BoundaryKind::inflow, value}};
}

wasl::NamedBoundary outflow(const std::string &name) {
  return {name, {wasl::BoundaryKind::outflow, 0.0}};
}

std::vector<std::pair<double, double>> history_time_mass(const fs::path &directory) {
  std::ifstream input(directory / "history.csv");
  std::string line;
  std::getline(input, line); // header
  std::vector<std::pair<double, double>> rows;
  while (std::getline(input, line)) {
    std::istringstream fields(line);
    std::string step, time, time_step, mass;
    std::getline(fields, step, ',');
    std::getline(fields, time, ',');
    std::getline(fields, time_step, ',');
    std::getline(fields, mass, ',');
    rows.emplace_back(std::stod(time), std::stod(mass));
  }
  return rows;
}

bool free_stream_is_preserved(const fs::path &mesh, const fs::path &root) {
  constexpr double state = 0.7;
  auto config = base_case(mesh, root / "free-stream");
  config.initial_default = state;
  config.boundaries = {inflow("left", state), outflow("right")};
  const auto checkpoint = wasl::read_checkpoint(wasl::run_scalar_advection_case(config).checkpoint);
  // Persisted layout U[mode][element] for one variable.
  for (std::size_t mode = 0; mode < checkpoint.mode_count; ++mode) {
    for (std::size_t element = 0; element < checkpoint.element_count; ++element) {
      const double expected = mode == 0 ? state * std::sqrt(2.0) : 0.0;
      const double actual = checkpoint.coefficients[mode * checkpoint.element_count + element];
      if (std::fabs(actual - expected) > kTolerance) {
        std::cerr << "free stream: mode " << mode << " element " << element << " is " << actual
                  << '\n';
        return false;
      }
    }
  }
  return true;
}

bool inflow_mass_balances(const fs::path &mesh, const fs::path &root, double wavespeed) {
  constexpr double inflow_value = 1.5;
  const bool left_inflow = wavespeed > 0.0;
  auto config = base_case(mesh, root / (left_inflow ? "inflow-left" : "inflow-right"));
  config.wavespeed = wavespeed;
  config.boundaries = left_inflow ? std::vector{inflow("left", inflow_value), outflow("right")}
                                  : std::vector{outflow("left"), inflow("right", inflow_value)};
  wasl::run_scalar_advection_case(config);
  const auto rows = history_time_mass(config.files.output_directory);
  if (rows.size() < 10)
    return false;
  for (const auto &[time, mass] : rows) {
    const double expected = std::fabs(wavespeed) * inflow_value * time;
    if (std::fabs(mass - expected) > kTolerance) {
      std::cerr << "inflow balance: mass " << mass << " expected " << expected << " at t=" << time
                << '\n';
      return false;
    }
  }
  return true;
}

bool fails_naming(const wasl::ScalarAdvectionCase &config, const std::string &expected) {
  try {
    wasl::run_scalar_advection_case(config);
  } catch (const std::exception &error) {
    if (std::string(error.what()).find(expected) != std::string::npos)
      return true;
    std::cerr << "unexpected error: " << error.what() << '\n';
    return false;
  }
  std::cerr << "expected an error naming " << expected << '\n';
  return false;
}

bool misconfigurations_are_named(const fs::path &mesh, const fs::path &root) {
  auto config = base_case(mesh, root / "invalid");
  config.boundaries = {inflow("left", 1.0)};
  if (!fails_naming(config, "boundaries.right"))
    return false;
  config.boundaries = {inflow("left", 1.0), outflow("right"), outflow("top")};
  if (!fails_naming(config, "boundaries.top"))
    return false;
  config.boundaries = {{"left", {wasl::BoundaryKind::periodic, 0.0}}, outflow("right")};
  if (!fails_naming(config, "periodic"))
    return false;
  config.boundaries = {inflow("left", 1.0), inflow("left", 2.0)};
  if (!fails_naming(config, "unique"))
    return false;
  config.boundaries = {inflow("left", std::nan("")), outflow("right")};
  return fails_naming(config, "boundary left value");
}

} // namespace

int main() {
  const auto root = fs::temp_directory_path() / "wasl-scalar-boundary-case-test";
  fs::remove_all(root);
  fs::create_directories(root);
  const auto mesh = root / "line.h5";
  wasl::write_prepared_mesh(mesh, wasl::prepare_mesh(wasl_test::named_line(24)));

  const bool ok = free_stream_is_preserved(mesh, root) && inflow_mass_balances(mesh, root, 1.0) &&
                  inflow_mass_balances(mesh, root, -1.0) && misconfigurations_are_named(mesh, root);
  fs::remove_all(root);
  return ok ? 0 : 1;
}
