"""Shared helpers for tutorial verifiers: meshes, CLI runs, and VTU output."""

import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def write_line_mesh(path, element_count):
    """Write a uniform unit-interval MSH 4.1 mesh with points "left" and "right"."""
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


def run(command):
    subprocess.run([str(part) for part in command], check=True, text=True, capture_output=True)


def copy_case(case, destination):
    """Copy a case without generated meshes or results."""
    shutil.copytree(case, destination, ignore=shutil.ignore_patterns("results", "*.h5", "*.msh"))
    return destination


def prepare_and_run(wasl, case, element_count):
    """Mesh the unit interval, run the case, and return its run directory."""
    source = case / "mesh/source.msh"
    write_line_mesh(source, element_count)
    run([wasl, "mesh", "prepare", source, case / "mesh/prepared.h5"])
    run([wasl, "run", case])
    return max((case / "results").glob("run-*"))


def line_cell_states(vtu):
    """Return (left, right, [components]) for every cell of a 1D VTU file."""
    piece = ET.parse(vtu).find(".//Piece")
    points = [float(value) for value in piece.find("./Points/DataArray").text.split()][::3]
    arrays = {array.get("Name"): array for array in piece.findall("./Cells/DataArray")}
    connectivity = [int(value) for value in arrays["connectivity"].text.split()]
    data = piece.find("./CellData/DataArray")
    width = int(data.get("NumberOfComponents"))
    values = [float(value) for value in data.text.split()]
    return [
        (
            points[connectivity[2 * cell]],
            points[connectivity[2 * cell + 1]],
            values[width * cell:width * (cell + 1)],
        )
        for cell in range(len(values) // width)
    ]


def line_cells(vtu):
    """Return (left, right, cell average) for every cell of a scalar 1D VTU file."""
    return [(left, right, state[0]) for left, right, state in line_cell_states(vtu)]


def final_solution(run_directory):
    return max(Path(run_directory).glob("solution-*.vtu"))


def history(run_directory):
    """Return (time, mass) for every history row."""
    rows = (Path(run_directory) / "history.csv").read_text().splitlines()[1:]
    return [(float(row.split(",")[1]), float(row.split(",")[3])) for row in rows if row]
