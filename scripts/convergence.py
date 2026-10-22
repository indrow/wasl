#!/usr/bin/env python3
import argparse
import math
import shutil
import subprocess
import tempfile
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path


def write_line_mesh(path, element_count):
    interior = list(range(2, element_count + 1))
    node_blocks = 3 if interior else 2
    lines = [
        "$MeshFormat", "4.1 0 8", "$EndMeshFormat",
        "$PhysicalNames", "3", '0 2 "left"', '0 3 "right"', '1 1 "fluid"',
        "$EndPhysicalNames", "$Entities", "2 1 0 0", "1 0 0 0 1 2",
        "2 1 0 0 1 3", "1 0 0 0 1 0 0 1 1 2 1 -2", "$EndEntities", "$Nodes",
        f"{node_blocks} {element_count + 1} 1 {element_count + 1}",
        "0 1 0 1", "1", "0 0 0", "0 2 0 1", str(element_count + 1), "1 0 0",
    ]
    if interior:
        lines += ["1 1 0 " + str(len(interior))]
        lines += [str(tag) for tag in interior]
        lines += [f"{(tag - 1) / element_count:.17g} 0 0" for tag in interior]
    lines += [
        "$EndNodes", "$Elements", f"3 {element_count + 2} 1 {element_count + 2}",
        "0 1 15 1", "1 1", "0 2 15 1", f"2 {element_count + 1}",
        f"1 1 1 {element_count}",
    ]
    lines += [f"{index + 3} {index + 1} {index + 2}" for index in range(element_count)]
    lines += ["$EndElements"]
    path.write_text("\n".join(lines) + "\n")


def cell_average_error(vtu, amplitude, velocity, time):
    piece = ET.parse(vtu).find(".//Piece")
    points = [float(value) for value in piece.find("./Points/DataArray").text.split()][::3]
    arrays = {array.get("Name"): array for array in piece.findall("./Cells/DataArray")}
    connectivity = [int(value) for value in arrays["connectivity"].text.split()]
    values = [float(value) for value in piece.find("./CellData/DataArray").text.split()]
    error_squared = 0.0
    for index, value in enumerate(values):
        left, right = (points[node] for node in connectivity[2 * index:2 * index + 2])
        width = right - left
        exact = amplitude * (
            math.cos(2 * math.pi * (left - velocity * time))
            - math.cos(2 * math.pi * (right - velocity * time))
        ) / (2 * math.pi * width)
        error_squared += width * (value - exact) ** 2
    return math.sqrt(error_squared)


def run(command):
    subprocess.run(command, check=True, text=True, capture_output=True)


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
        work = Path(temporary) / "case"
        shutil.copytree(case, work, ignore=shutil.ignore_patterns("results", "*.h5", "*.msh"))
        for element_count in args.levels:
            source = work / "mesh/source.msh"
            prepared = work / "mesh/prepared.h5"
            write_line_mesh(source, element_count)
            run([str(args.wasl.resolve()), "mesh", "prepare", str(source), str(prepared)])
            run([str(args.wasl.resolve()), "run", str(work)])
            result = max((work / "results").glob("run-*"))
            final_vtu = max(result.glob("solution-*.vtu"))
            errors.append(cell_average_error(final_vtu, amplitude, velocity, end_time))

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
