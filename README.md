# Wasl

Wasl is a CPU research solver for high-order discontinuous Galerkin methods
with WENO limiting on unstructured meshes.

It runs modal RKDG with SSP-RK3 on unstructured 1D Gmsh line meshes and
straight-sided 2D triangles, for scalar advection and compressible Euler with
characteristic WENO and positivity limiting. Two-rank MPI execution is
available for 1D scalar advection.

## Quick start

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
```

Or build the same smoke target in Docker:

```bash
docker compose build
docker compose run --rm dev
```

The Rust CLI validates cases, prepares meshes, creates run manifests, inspects
and restarts checkpoints, and exports checkpoint cell averages to VTU:

```bash
cargo run -p wasl-cli -- check cases/tutorials/advection/line-pulse
cargo run -p wasl-cli -- mesh prepare source.msh prepared.h5
```

Two-rank CPU cases use the opt-in MPI bridge and `parallel.type = "mpi"`:

```bash
cargo build -p wasl-cli --features mpi
mpirun -np 2 target/debug/wasl-cli run path/to/case
```

## Project boundaries

- C++20 owns numerical state and runs complete simulations.
- Rust is the CLI and crosses into C++ only for coarse operations.
- Gmsh MSH 4.1 is input; prepared HDF5 is the runtime mesh.
- TOML configures cases; HDF5 checkpoints state; VTK XML visualizes it.
- Distributed CPU execution currently uses exactly two MPI ranks.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a pull request.

## Documentation

- [`docs/BUILDING.md`](docs/BUILDING.md): dependencies, build, and presets.
- [`docs/CASE_FORMAT.md`](docs/CASE_FORMAT.md): case directory and TOML files.
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md): boundaries, data flow, and
  design decisions.
- [`docs/NUMERICAL_CONVENTIONS.md`](docs/NUMERICAL_CONVENTIONS.md): bases,
  orientation, and state layout.
- [`docs/CODING_STYLE.md`](docs/CODING_STYLE.md): language, test, and
  verification conventions.
- [`docs/PROJECT_STATE.md`](docs/PROJECT_STATE.md) and
  [`docs/ROADMAP.md`](docs/ROADMAP.md): current capability and next work.
