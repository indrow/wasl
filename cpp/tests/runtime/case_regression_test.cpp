// Case-level numerical regression over every serial driver.
//
// Scalar advection (1D and 2D, P0-P3, with restart) and Euler (1D and 2D,
// P0-P2, with and without limiting) run on irregular meshes built in code. Each
// run is reduced to compact metrics: step count, final time, the sum and L2
// norm of the final checkpoint coefficients, and the final history mass. The
// metrics are compared against a stored reference file.
//
// Usage: wasl-case-regression-test <reference> [--print]
// --print writes the current metrics in reference format instead of checking.
// Regenerate the reference only for an intended numerical change, and say why
// in the commit.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/case/euler_case.hpp"
#include "wasl/case/scalar_advection_case.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/prepared_mesh_file.hpp"
#include "wasl/mesh/line/build_line_mesh.hpp"
#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/prepared/prepare_mesh.hpp"
#include "wasl/mesh/triangle/build_triangle_mesh.hpp"

namespace fs = std::filesystem;

namespace {

// Double precision with a deterministic CPU path reproduces these metrics to
// round-off. The allowance absorbs compiler and optimization-level reordering
// across about a hundred steps, yet any change to a flux, basis, quadrature,
// limiter, or time step moves them by many orders of magnitude more.
constexpr double kRelativeTolerance = 1e-10;

struct Metrics {
  std::uint64_t step = 0;
  double time = 0.0;
  double coefficient_sum = 0.0;
  double coefficient_l2 = 0.0;
  double final_mass = 0.0;
};

double last_history_mass(const fs::path &history) {
  std::ifstream input(history);
  std::string line, last;
  while (std::getline(input, line))
    if (!line.empty())
      last = line;
  return std::stod(last.substr(last.rfind(',') + 1));
}

Metrics measure(const wasl::SimulationSummary &summary) {
  const auto checkpoint = wasl::read_checkpoint(summary.checkpoint);
  Metrics metrics{summary.step, summary.time};
  for (const double value : checkpoint.coefficients) {
    metrics.coefficient_sum += value;
    metrics.coefficient_l2 += value * value;
  }
  metrics.coefficient_l2 = std::sqrt(metrics.coefficient_l2);
  metrics.final_mass = last_history_mass(summary.checkpoint.parent_path() / "history.csv");
  return metrics;
}

// Nodes perturbed off a uniform grid; elements recorded in shuffled order and
// orientation to exercise deterministic mesh normalization.
wasl::LineMesh irregular_line(int count) {
  wasl::MshDocument document;
  document.point_entities.push_back({1, 0.0, 0.0, 0.0, {}});
  document.curve_entities.push_back({1, 0, 0, 0, 1, 0, 0, {1}, {}});
  for (int node = 0; node <= count; ++node) {
    double x = static_cast<double>(node) / count;
    if (node > 0 && node < count)
      x += 0.3 / count * std::sin(7.0 * node);
    document.nodes.push_back({node + 1, 1, 1, x, 0.0, 0.0});
  }
  for (int element = count - 1; element >= 0; --element) {
    const std::int64_t first = element + 1;
    const std::int64_t second = element + 2;
    document.elements.push_back(
        {element + 1, 1, 1, wasl::msh_line_element,
         element % 2 ? std::vector{second, first} : std::vector{first, second}});
  }
  return wasl::build_line_mesh(document);
}

// Unit box with perturbed interior nodes and alternating diagonals.
wasl::TriangleMesh irregular_box(int n) {
  wasl::MshDocument document;
  document.surface_entities.push_back({1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, {1}, {}});
  const auto id = [n](int i, int j) { return static_cast<std::int64_t>(j * (n + 1) + i + 1); };
  for (int j = 0; j <= n; ++j) {
    for (int i = 0; i <= n; ++i) {
      double x = static_cast<double>(i) / n;
      double y = static_cast<double>(j) / n;
      if (i > 0 && i < n && j > 0 && j < n) {
        x += 0.2 / n * std::sin(3.0 * i + 5.0 * j);
        y += 0.2 / n * std::cos(2.0 * i - 7.0 * j);
      }
      document.nodes.push_back({id(i, j), 2, 1, x, y, 0.0});
    }
  }
  std::int64_t tag = 1;
  const auto add = [&](std::int64_t a, std::int64_t b, std::int64_t c) {
    document.elements.push_back({tag++, 2, 1, wasl::msh_triangle_element, {a, b, c}});
  };
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if ((i + j) % 2 == 0) {
        add(id(i, j), id(i + 1, j), id(i + 1, j + 1));
        add(id(i, j), id(i, j + 1), id(i + 1, j + 1));
      } else {
        add(id(i, j), id(i + 1, j), id(i, j + 1));
        add(id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
      }
    }
  }
  return wasl::build_triangle_mesh(document);
}

wasl::RunControl run_control(double end_time, std::uint64_t keep_checkpoints) {
  return {0.0, end_time, 1000, 3, 4, 5, keep_checkpoints};
}

std::map<std::string, Metrics> run_cases(const fs::path &root) {
  const auto line_path = root / "line.h5";
  const auto box_path = root / "box.h5";
  wasl::write_prepared_mesh(line_path, wasl::prepare_mesh(irregular_line(12)));
  wasl::write_prepared_mesh(box_path, wasl::prepare_mesh(irregular_box(5)));
  std::map<std::string, Metrics> results;

  for (const int dimensions : {1, 2}) {
    for (const int order : {0, 1, 2, 3}) {
      const auto name = "scalar-d" + std::to_string(dimensions) + "-p" + std::to_string(order);
      wasl::ScalarAdvectionCase config;
      config.files = {dimensions == 1 ? line_path : box_path, root / name, {}};
      config.run = run_control(0.2, 2);
      config.dimensions = dimensions;
      config.polynomial_order = order;
      config.wavespeed = 0.7;
      config.velocity_y = dimensions == 2 ? -0.4 : 0.0;
      config.cfl = 0.3;
      config.initial_default = 0.25;
      config.initial_sine = dimensions == 1 && order == 2;
      config.initial_sine_amplitude = 1.5;
      config.initial_intervals.push_back({0.2, 0.55, 1.0});
      config.initial_circles.push_back({0.5, 0.4, 0.25, 2.0});
      const auto summary = wasl::run_scalar_advection_case(config);
      results[name] = measure(summary);

      auto restart = config;
      restart.files = {config.files.prepared_mesh, root / (name + "-restart"), summary.checkpoint};
      restart.run = run_control(0.3, 2);
      results[name + "-restart"] = measure(wasl::run_scalar_advection_case(restart));
    }
  }

  for (const int dimensions : {1, 2}) {
    for (const int order : {0, 1, 2}) {
      for (const bool limiting : {false, true}) {
        const auto name = "euler-d" + std::to_string(dimensions) + "-p" + std::to_string(order) +
                          (limiting ? "-limited" : "");
        wasl::EulerCase config;
        config.files = {dimensions == 1 ? line_path : box_path, root / name, {}};
        config.run = {0.0, 0.03, 1000, 2, 3, 4, 1};
        config.limiting = limiting;
        config.sensor_threshold = 1e-4;
        config.dimensions = dimensions;
        config.polynomial_order = order;
        config.cfl = 0.2;
        config.initial_regions.push_back({wasl::EulerRegion::Shape::default_region,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          {1.0, 0.3, dimensions == 2 ? 0.2 : 0.0, 1.0}});
        // A strong jump needs the limiter; unlimited runs get a mild one.
        const wasl::EulerPrimitive inside =
            limiting ? wasl::EulerPrimitive{0.2, 0.0, 0.0, 0.15}
                     : wasl::EulerPrimitive{0.9, 0.25, dimensions == 2 ? 0.1 : 0.0, 0.95};
        if (dimensions == 1)
          config.initial_regions.push_back(
              {wasl::EulerRegion::Shape::interval, 0.3, 0.6, 0, 0, 0, inside});
        else
          config.initial_regions.push_back(
              {wasl::EulerRegion::Shape::circle, 0, 0, 0.5, 0.5, 0.25, inside});
        results[name] = measure(wasl::run_euler_case(config));
      }
    }
  }
  return results;
}

void print(const std::map<std::string, Metrics> &results) {
  std::printf("# case step time coefficient_sum coefficient_l2 final_mass\n");
  for (const auto &[name, value] : results)
    std::printf("%s %llu %.17g %.17g %.17g %.17g\n", name.c_str(),
                static_cast<unsigned long long>(value.step), value.time, value.coefficient_sum,
                value.coefficient_l2, value.final_mass);
}

std::map<std::string, Metrics> read_reference(const fs::path &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("case regression: cannot open reference " + path.string());
  std::map<std::string, Metrics> reference;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#')
      continue;
    std::istringstream fields(line);
    std::string name;
    Metrics metrics;
    if (!(fields >> name >> metrics.step >> metrics.time >> metrics.coefficient_sum >>
          metrics.coefficient_l2 >> metrics.final_mass))
      throw std::runtime_error("case regression: malformed reference line: " + line);
    reference[name] = metrics;
  }
  return reference;
}

bool close(double actual, double expected) {
  return std::fabs(actual - expected) <= kRelativeTolerance * std::max(1.0, std::fabs(expected));
}

bool matches(const std::string &name, const Metrics &actual, const Metrics &expected) {
  const bool ok = actual.step == expected.step && close(actual.time, expected.time) &&
                  close(actual.coefficient_sum, expected.coefficient_sum) &&
                  close(actual.coefficient_l2, expected.coefficient_l2) &&
                  close(actual.final_mass, expected.final_mass);
  if (!ok)
    std::fprintf(stderr,
                 "%s: step %llu/%llu time %.17g/%.17g sum %.17g/%.17g l2 %.17g/%.17g "
                 "mass %.17g/%.17g (actual/expected)\n",
                 name.c_str(), static_cast<unsigned long long>(actual.step),
                 static_cast<unsigned long long>(expected.step), actual.time, expected.time,
                 actual.coefficient_sum, expected.coefficient_sum, actual.coefficient_l2,
                 expected.coefficient_l2, actual.final_mass, expected.final_mass);
  return ok;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <reference> [--print]\n", argv[0]);
    return 2;
  }
  const auto root = fs::temp_directory_path() / "wasl-case-regression-test";
  fs::remove_all(root);
  fs::create_directories(root);
  try {
    const auto results = run_cases(root);
    fs::remove_all(root);
    if (argc > 2 && std::string(argv[2]) == "--print") {
      print(results);
      return 0;
    }
    const auto reference = read_reference(argv[1]);
    bool ok = reference.size() == results.size();
    if (!ok)
      std::fprintf(stderr, "case regression: %zu cases ran, reference lists %zu\n", results.size(),
                   reference.size());
    for (const auto &[name, expected] : reference) {
      const auto actual = results.find(name);
      if (actual == results.end()) {
        std::fprintf(stderr, "%s: missing from this run\n", name.c_str());
        ok = false;
      } else if (!matches(name, actual->second, expected)) {
        ok = false;
      }
    }
    return ok ? 0 : 1;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "case regression: %s\n", error.what());
    return 1;
  }
}
