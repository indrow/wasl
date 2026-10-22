#!/usr/bin/env python3
"""Verify a 1D Euler shock-tube case against the exact Riemann solution.

Checks on uniform unit-interval meshes:
  refinement   The L1 error of the density cell averages against the exact
               solution decreases by at least --minimum-ratio per doubling,
               and the finest error is at most --maximum-error.
  positivity   Every cell average has positive density and pressure.
  free stream  A uniform moving state with the same outflow ends is
               preserved to round-off.
"""

import argparse
import math
import tempfile
import tomllib
from pathlib import Path

from wasl_case_tools import copy_case, final_solution, line_cell_states, prepare_and_run

# Round-off for O(1) states over a few hundred steps.
ROUND_OFF = 1e-12

# Sub-cell samples for exact cell averages. A discontinuity inside a cell adds
# at most jump * h / SAMPLES to the averaged error, far below the scheme error.
SAMPLES = 400


def load(case, name):
    with (case / name).open("rb") as stream:
        return tomllib.load(stream)


class ExactRiemann:
    """Exact ideal-gas Riemann solution (Toro, Riemann Solvers, chapter 4)."""

    def __init__(self, gamma, left, right):
        self.gamma = gamma
        self.left = left    # (density, velocity, pressure)
        self.right = right
        self.sound_left = math.sqrt(gamma * left[2] / left[0])
        self.sound_right = math.sqrt(gamma * right[2] / right[0])
        self.pressure_star = self._star_pressure()
        self.velocity_star = 0.5 * (left[1] + right[1]) + 0.5 * (
            self._wave(self.pressure_star, right, self.sound_right)
            - self._wave(self.pressure_star, left, self.sound_left)
        )

    def _wave(self, pressure, state, sound):
        """Velocity change across a shock or rarefaction to the star pressure."""
        density, _, side_pressure = state
        gamma = self.gamma
        if pressure > side_pressure:
            a = 2.0 / ((gamma + 1.0) * density)
            b = (gamma - 1.0) / (gamma + 1.0) * side_pressure
            return (pressure - side_pressure) * math.sqrt(a / (pressure + b))
        exponent = (gamma - 1.0) / (2.0 * gamma)
        return 2.0 * sound / (gamma - 1.0) * ((pressure / side_pressure) ** exponent - 1.0)

    def _star_pressure(self):
        def residual(pressure):
            return (self._wave(pressure, self.left, self.sound_left)
                    + self._wave(pressure, self.right, self.sound_right)
                    + self.right[1] - self.left[1])

        lower, upper = 1e-14, 10.0 * max(self.left[2], self.right[2])
        for _ in range(200):
            middle = 0.5 * (lower + upper)
            if residual(middle) > 0.0:
                upper = middle
            else:
                lower = middle
        return 0.5 * (lower + upper)

    def sample(self, speed):
        """(density, velocity, pressure) at similarity coordinate x / t."""
        gamma = self.gamma
        if speed <= self.velocity_star:
            state, sound, sign = self.left, self.sound_left, 1.0
        else:
            state, sound, sign = self.right, self.sound_right, -1.0
        density, velocity, pressure = state
        star_ratio = self.pressure_star / pressure
        if self.pressure_star > pressure:  # shock
            shock = velocity - sign * sound * math.sqrt(
                (gamma + 1.0) / (2.0 * gamma) * star_ratio + (gamma - 1.0) / (2.0 * gamma)
            )
            if sign * (speed - shock) <= 0.0:
                return state
            g = (gamma - 1.0) / (gamma + 1.0)
            return (density * (star_ratio + g) / (g * star_ratio + 1.0),
                    self.velocity_star, self.pressure_star)
        # rarefaction
        star_sound = sound * star_ratio ** ((gamma - 1.0) / (2.0 * gamma))
        head = velocity - sign * sound
        tail = self.velocity_star - sign * star_sound
        if sign * (speed - head) <= 0.0:
            return state
        if sign * (speed - tail) >= 0.0:
            return (density * star_ratio ** (1.0 / gamma), self.velocity_star, self.pressure_star)
        factor = 2.0 / (gamma + 1.0) + sign * (gamma - 1.0) / ((gamma + 1.0) * sound) * (
            velocity - speed
        )
        return (density * factor ** (2.0 / (gamma - 1.0)),
                2.0 / (gamma + 1.0) * (sign * sound + (gamma - 1.0) / 2.0 * velocity + speed),
                pressure * factor ** (2.0 * gamma / (gamma - 1.0)))


def primitive(region):
    return (region["density"], region["velocity"][0], region["pressure"])


def case_parameters(case):
    gamma = load(case, "physics.toml")["equation"]["heat_capacity_ratio"]
    regions = load(case, "initial/state.toml")["regions"]
    driver = [region for region in regions if region["where"] == "default"]
    driven = [region for region in regions if region["where"] == "interval"]
    if len(driver) != 1 or len(driven) != 1 or driven[0]["upper"][0] < 1.0:
        raise SystemExit("case must be a default region plus one interval reaching x = 1")
    interface = driven[0]["lower"][0]
    end_time = load(case, "solver.toml")["run"]["end_time"]
    return gamma, primitive(driver[0]), primitive(driven[0]), interface, end_time


def density_l1_error(vtu, exact, interface, time):
    error = 0.0
    for left, right, state in line_cell_states(vtu):
        width = right - left
        average = sum(
            exact.sample((left + (sample + 0.5) * width / SAMPLES - interface) / time)[0]
            for sample in range(SAMPLES)
        ) / SAMPLES
        error += width * abs(state[0] - average)
    return error


def minimum_density_and_pressure(vtu, gamma):
    densities, pressures = [], []
    for _, _, (density, momentum_x, momentum_y, energy) in line_cell_states(vtu):
        densities.append(density)
        pressures.append((gamma - 1.0) * (energy - 0.5 * (momentum_x**2 + momentum_y**2) / density))
    return min(densities), min(pressures)


def write_uniform_state(case, density, velocity, pressure):
    (case / "initial/state.toml").write_text(
        'format_version = 1\nrepresentation = "primitive"\n\n'
        f'[[regions]]\nname = "free_stream"\nwhere = "default"\n'
        f"density = {density!r}\nvelocity = [{velocity!r}]\npressure = {pressure!r}\n"
    )


def free_stream_error(wasl, case, directory, gamma, elements):
    density, velocity, pressure = 1.0, 0.5, 1.0
    work = copy_case(case, directory)
    write_uniform_state(work, density, velocity, pressure)
    result = prepare_and_run(wasl, work, elements)
    energy = pressure / (gamma - 1.0) + 0.5 * density * velocity**2
    expected = (density, density * velocity, 0.0, energy)
    return max(
        abs(value - target)
        for _, _, state in line_cell_states(final_solution(result))
        for value, target in zip(state, expected)
    )


def main():
    parser = argparse.ArgumentParser(description="Verify a 1D Euler shock-tube case")
    parser.add_argument("wasl", type=Path)
    parser.add_argument("case", type=Path)
    parser.add_argument("--levels", type=int, nargs="+", default=[50, 100, 200])
    # Shocks and contacts limit L1 convergence to first order, an error ratio of
    # 2 per doubling. On the tutorial case the observed ratios are 1.91, 1.90,
    # and 2.08 from 50 to 400 elements, with an L1 error of 3.0e-3 at 200.
    parser.add_argument("--minimum-ratio", type=float, default=1.6)
    parser.add_argument("--maximum-error", type=float, default=4e-3)
    args = parser.parse_args()

    case = args.case.resolve()
    wasl = args.wasl.resolve()
    gamma, left, right, interface, end_time = case_parameters(case)
    exact = ExactRiemann(gamma, left, right)

    errors, minima = [], []
    with tempfile.TemporaryDirectory(prefix="wasl-sod-") as temporary:
        work = copy_case(case, Path(temporary) / "case")
        for elements in args.levels:
            vtu = final_solution(prepare_and_run(wasl, work, elements))
            errors.append(density_l1_error(vtu, exact, interface, end_time))
            minima.append(minimum_density_and_pressure(vtu, gamma))
        free_stream = free_stream_error(
            wasl, case, Path(temporary) / "free-stream", gamma, args.levels[0]
        )

    ratios = [errors[index - 1] / errors[index] for index in range(1, len(errors))]
    print("elements,density_L1_error,ratio,minimum_density,minimum_pressure")
    for index, (elements, error) in enumerate(zip(args.levels, errors)):
        ratio = "" if index == 0 else f"{ratios[index - 1]:.3f}"
        print(f"{elements},{error:.6e},{ratio},{minima[index][0]:.6e},{minima[index][1]:.6e}")
    print(f"free_stream,{free_stream:.3e}")

    failures = []
    if len(ratios) < 1 or min(ratios) < args.minimum_ratio:
        failures.append(f"refinement ratio below {args.minimum_ratio}")
    if not errors[-1] <= args.maximum_error:
        failures.append(f"finest L1 error above {args.maximum_error:.1e}")
    if not all(density > 0.0 and pressure > 0.0 for density, pressure in minima):
        failures.append("nonpositive cell average")
    if not free_stream <= ROUND_OFF:
        failures.append("free stream not preserved")
    if failures:
        raise SystemExit("shock-tube verification failed: " + "; ".join(failures))
    print("Shock tube verified")


if __name__ == "__main__":
    main()
