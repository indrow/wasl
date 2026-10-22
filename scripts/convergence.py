#!/usr/bin/env python3
import argparse
import math
import tempfile
import tomllib
from pathlib import Path

from wasl_case_tools import copy_case, final_solution, line_cells, prepare_and_run


def cell_average_error(vtu, amplitude, velocity, time):
    error_squared = 0.0
    for left, right, value in line_cells(vtu):
        width = right - left
        exact = amplitude * (
            math.cos(2 * math.pi * (left - velocity * time))
            - math.cos(2 * math.pi * (right - velocity * time))
        ) / (2 * math.pi * width)
        error_squared += width * (value - exact) ** 2
    return math.sqrt(error_squared)


def main():
    parser = argparse.ArgumentParser(description="Run a scalar-advection convergence check")
    parser.add_argument("wasl", type=Path)
    parser.add_argument("case", type=Path)
    parser.add_argument("--levels", type=int, nargs="+", default=[8, 16, 32])
    parser.add_argument("--minimum-rate", type=float, default=2.7)
    args = parser.parse_args()

    case = args.case.resolve()
    with (case / "physics.toml").open("rb") as stream:
        velocity = tomllib.load(stream)["equation"]["velocity"][0]
    with (case / "solver.toml").open("rb") as stream:
        end_time = tomllib.load(stream)["run"]["end_time"]
    with (case / "initial/state.toml").open("rb") as stream:
        regions = tomllib.load(stream)["regions"]
    sine = [region for region in regions if region["where"] == "sine"]
    if len(sine) != 1:
        parser.error("case must contain exactly one sine initial region")
    amplitude = sine[0]["value"]

    errors = []
    with tempfile.TemporaryDirectory(prefix="wasl-convergence-") as temporary:
        work = copy_case(case, Path(temporary) / "case")
        for element_count in args.levels:
            result = prepare_and_run(args.wasl.resolve(), work, element_count)
            errors.append(cell_average_error(final_solution(result), amplitude, velocity, end_time))

    rates = [
        math.log(errors[index - 1] / errors[index])
        / math.log(args.levels[index] / args.levels[index - 1])
        for index in range(1, len(errors))
    ]
    print("elements,L2_cell_average_error,rate")
    for index, (element_count, error) in enumerate(zip(args.levels, errors)):
        rate = "" if index == 0 else f"{rates[index - 1]:.3f}"
        print(f"{element_count},{error:.6e},{rate}")
    if len(rates) < 2 or min(rates) < args.minimum_rate:
        raise SystemExit(f"convergence failed: minimum rate {min(rates, default=0):.3f} < {args.minimum_rate}")
    print(f"P2 convergence verified: minimum rate {min(rates):.3f}")


if __name__ == "__main__":
    main()
