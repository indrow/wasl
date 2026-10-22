# Current Project State

Last updated: 2026-10-22
Current milestone: M9 (reproducible 1D and 2D tutorials), active

## Working

### Repository and tooling

- Repository instructions, architecture, conventions, and roadmap exist.
- CMake builds the C++20 core library, a smoke executable, and registered
  numerical tests through the `cpu-debug` and `cpu-mpi` presets.
- Docker provides a reproducible development environment with HDF5 and
  Open MPI. CI runs the C++ and Rust tests and checks formatting.
- The C++ core is organized one concept per file in subsystem subdirectories
  (see `ARCHITECTURE.md`), with private headers beside their sources.

### Numerics

- 1D Legendre and 2D Dubiner modal bases, line and triangle quadrature, affine
  geometry, Gmsh MSH 4.1 import, and canonical line and triangle meshes.
- Deterministic scalar-advection and Euler DG residuals on lines and straight
  triangles, with upwind, Rusanov, and HLLC fluxes.
- SSP-RK3 time integration with CFL time steps. Periodic, inflow, outflow, and
  slip-wall boundary policies exist at the operator level.
- Euler limiting: Persson-style modal-decay sensor, stable troubled-cell
  compaction, neighbor polynomial transform, characteristic simple WENO, and
  density and pressure positivity limiting. Mode-zero coefficients are
  preserved.
- Irregular-mesh convergence and conservation are covered by regressions.
- `wasl.runtime.case-regression` pins every serial case driver (scalar 1D/2D
  P0-P3 with restart, Euler 1D/2D P0-P2 with and without limiting) to stored
  step counts, final times, coefficient norms, and masses at relative `1e-10`.

### Case workflow

- The Rust case loader merges the case TOML files into one checked model,
  names the owning field on every validation failure, and writes a
  deterministic `effective.toml`.
- The CLI provides `check`, `run`, `restart`, `mesh prepare`, `inspect`, and
  `export`. A coarse `cxx` bridge keeps numerical work in C++; Rust passes
  paths and summaries only.
- Prepared HDF5 meshes carry connectivity, metrics, face ownership, physical
  tags, and partition metadata. Versioned HDF5 checkpoints, CSV history, and
  VTU/PVD output are written without linking VTK.
- `wasl run` and `wasl restart` dispatch periodic 1D and 2D scalar-advection
  and Euler cases, including Euler limiting when enabled. Serial 1D cases also
  accept outflow boundaries, and scalar advection inflow boundaries, matched to
  mesh physical names.

### Distributed execution

- Deterministic contiguous-block partitioning with one face-adjacent ghost
  layer (see `ARCHITECTURE.md`).
- A transport-agnostic halo schedule with an in-process oracle, and a
  host-staged nonblocking MPI transport behind `WASL_ENABLE_MPI`.
- Two-rank 1D scalar-advection cases run with `parallel.type = "mpi"` and the
  Cargo `mpi` feature. Ranks write owned VTU and checkpoint pieces; rank 0
  writes PVTU, PVD, history, and the effective configuration. Restart
  reassembles rank-local checkpoints.

### Tutorials

- `line-pulse`: 1D periodic scalar advection with run and restart coverage.
- `rotating-cone`: 2D periodic P2 scalar advection with a conservation check.
- `smooth-periodic`: 1D P2 convergence study; every observed rate must exceed
  2.7.
- `line-inflow`: 1D P2 inflow/outflow advection; mass balance, transport of
  the inflow state, and free-stream preservation are verified to `1e-12`.
- `sod-shock-tube`: 1D P2 limited Euler with outflow ends; density L1 error
  against the exact Riemann solution halves per mesh doubling (3.0e-3 at 200
  elements), and a uniform moving state is preserved to `1e-12`.

## Not implemented

- Euler inflow and slip-wall boundaries; non-periodic boundaries for 2D cases
  and MPI runs; scalar-advection limiting; and Euler with MPI at case level.
- 2D case-level MPI.

## Known constraints

- Distributed CPU cases use exactly two MPI ranks.
- Gmsh is not included in the Docker image.
- Generated `.msh`, `.h5`, build, target, and result files are ignored.

## Next recommended task

Add a smooth 1D Euler convergence tutorial: an advected density wave on the
periodic unit interval.

Acceptance criteria:

- Euler initial regions accept a 1D `sine` density perturbation,
  `density = rho0 + value * sin(2 pi x)` at uniform velocity and pressure, and
  the configuration rejects it for 2D cases.
- The exact solution is the translated wave; P2 density cell-average L2 rates
  on 8, 16, and 32 elements are at least 2.7 without limiting.
- `wasl.runtime.case-regression` passes unchanged.
- The tutorial README provides reproducible mesh preparation, run, and
  verification commands covered by an automated CLI check.
