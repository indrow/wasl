#!/usr/bin/env python3
"""Verify a 1D scalar inflow/outflow advection case through the CLI.

Three checks on a uniform unit-interval mesh:
  balance      With the front still far from the outflow end, the total mass
               equals background + |a| (inflow - background) t to round-off,
               which holds only if the inflow value enters every SSP-RK3 stage.
  transport    After the front has left through the outflow end, every cell
               average equals the exact transported state, the inflow value.
  free stream  Starting from the inflow value, the state never changes.
"""

import argparse
import tempfile
import tomllib
from pathlib import Path

from wasl_case_tools import copy_case, final_solution, history, line_cells, prepare_and_run

# Round-off for O(1) states over a few thousand steps.
ROUND_OFF = 1e-12

# The upwind DG solution carries a precursor ahead of the front that decays
# about tenfold per half element. Mass balance is exact while the front stays
# at least this fraction of the domain away from the outflow end.
BALANCE_CLEARANCE = 0.8


def load(case, name):
    with (case / name).open("rb") as stream:
        return tomllib.load(stream)


def case_parameters(case):
    velocity = load(case, "physics.toml")["equation"]["velocity"][0]
    boundaries = load(case, "mesh/boundaries.toml")["boundaries"]
    inflow_end, outflow_end = ("left", "right") if velocity > 0 else ("right", "left")
    if boundaries[inflow_end]["type"] != "inflow" or boundaries[outflow_end]["type"] != "outflow":
        raise SystemExit(f"case must have inflow on the {inflow_end} end and outflow on the other")
    regions = load(case, "initial/state.toml")["regions"]
    if len(regions) != 1 or regions[0]["where"] != "default":
        raise SystemExit("case must start from a single constant default region")
    end_time = load(case, "solver.toml")["run"]["end_time"]
    if abs(velocity) * end_time <= 1.0:
        raise SystemExit("run.end_time must let the front leave the unit interval")
    return velocity, boundaries[inflow_end]["value"], regions[0]["value"]


def write_constant_initial_state(case, value):
    (case / "initial/state.toml").write_text(
        'format_version = 1\nrepresentation = "scalar"\n\n'
        f'[[regions]]\nname = "free_stream"\nwhere = "default"\nvalue = {value!r}\n'
    )


def maximum_deviation(run_directory, state):
    return max(abs(value - state) for _, _, value in line_cells(final_solution(run_directory)))


def main():
    parser = argparse.ArgumentParser(description="Verify a 1D inflow/outflow advection case")
    parser.add_argument("wasl", type=Path)
    parser.add_argument("case", type=Path)
    parser.add_argument("--elements", type=int, default=40)
    # The front's numerical tail behind it must also have left: at 40 elements
    # and t = 1.5 the deviation is 8e-14, while 20 elements still leave 3e-9.
    parser.add_argument("--transport-tolerance", type=float, default=ROUND_OFF)
    args = parser.parse_args()

    case = args.case.resolve()
    wasl = args.wasl.resolve()
    velocity, inflow, background = case_parameters(case)
    speed = abs(velocity)

    with tempfile.TemporaryDirectory(prefix="wasl-inflow-") as temporary:
        work = copy_case(case, Path(temporary) / "case")
        result = prepare_and_run(wasl, work, args.elements)
        balance_rows = [
            (time, mass) for time, mass in history(result)
            if speed * time <= 1.0 - BALANCE_CLEARANCE
        ]
        balance_error = max(
            abs(mass - (background + speed * (inflow - background) * time))
            for time, mass in balance_rows
        )
        transport_error = maximum_deviation(result, inflow)

        steady = copy_case(case, Path(temporary) / "free-stream")
        write_constant_initial_state(steady, inflow)
        steady_result = prepare_and_run(wasl, steady, args.elements)
        free_stream_error = max(
            maximum_deviation(steady_result, inflow),
            max(abs(mass - inflow) for _, mass in history(steady_result)),
        )

    print("check,error,tolerance")
    checks = [
        ("balance", balance_error, ROUND_OFF),
        ("transport", transport_error, args.transport_tolerance),
        ("free_stream", free_stream_error, ROUND_OFF),
    ]
    for name, error, tolerance in checks:
        print(f"{name},{error:.3e},{tolerance:.0e}")
    if len(balance_rows) < 10:
        raise SystemExit("balance check needs at least 10 history rows before the front nears outflow")
    failed = [name for name, error, tolerance in checks if not error <= tolerance]
    if failed:
        raise SystemExit("inflow verification failed: " + ", ".join(failed))
    print("Inflow advection verified")


if __name__ == "__main__":
    main()
