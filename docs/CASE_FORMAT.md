# Case Format

A case directory has one entry point and small included TOML files:

```text
case/
├── case.toml
├── physics.toml
├── numerics.toml
├── solver.toml
├── mesh/
│   ├── geometry.geo
│   ├── source.msh
│   ├── prepared.h5
│   └── boundaries.toml
├── initial/state.toml
└── results/
```

Every TOML file starts with `format_version = 1`. Paths are relative to the
case directory. `case.toml` names the prepared mesh and included files. A run
must validate all files before allocating numerical state and must write the
merged effective configuration to a new result directory.

Normal users select compiled equations, fluxes, basis/order, time integrator,
boundaries, and output intervals by name. v0.1 does not accept arbitrary PDE
expressions or scripts. Initial regions use a small controlled set: `interval`,
`circle`, `sine`, and `default`. The 1D `sine` profile is
`value * sin(2 pi x)`.

See `cases/tutorials/advection/line-pulse` for the first bootstrap example,
`smooth-periodic` for executable convergence verification, `line-inflow` for
inflow and outflow boundaries, `rotating-cone` for the later 2D case, and
`cases/tutorials/euler/sod-shock-tube` for limited Euler with outflow ends.
Generate ignored `prepared.h5` runtime caches with `wasl mesh prepare` after
producing a Gmsh MSH 4.1 file.

## Boundaries

`mesh/boundaries.toml` assigns one condition to each boundary physical group,
keyed by its Gmsh physical name:

```toml
format_version = 1

[boundaries.left]
type = "inflow"
value = 1.0

[boundaries.right]
type = "outflow"
```

- `periodic` wraps to the opposite boundary; name it with `pair`. Periodic
  conditions apply to both ends of a line together.
- `inflow` prescribes a fixed exterior state. Scalar advection requires a
  finite `value`; where the flow leaves the domain the value is unused.
- `outflow` copies the interior trace, so waves leave without reflection.
- `slip_wall` is accepted by validation for Euler but not yet dispatched.

Serial 1D cases dispatch `outflow`; scalar advection also dispatches `inflow`.
Other cases support `periodic` only.

At run time every boundary face's physical group must have a condition and
every configured name must match a boundary face; either mistake fails with the
boundary name.

## Commands

```bash
wasl check path/to/case
wasl mesh prepare path/to/source.msh path/to/prepared.h5
wasl run path/to/case
wasl restart path/to/case path/to/checkpoint.h5
wasl inspect path/to/checkpoint.h5
wasl export path/to/checkpoint.h5 path/to/prepared.h5 output.vtu
```

`check` validates configuration without requiring the prepared mesh. `run` and
`restart` require a readable prepared mesh and create the next
`results/run-NNNNNN` directory atomically. Each directory contains the fully
resolved `effective.toml`; restart directories also record checkpoint
provenance for serial runs. A completed serial run also contains `history.csv`,
numbered HDF5 checkpoints and VTU snapshots, and `solution.pvd`.

The serial case driver supports periodic 1D line and 2D straight-triangle
meshes, outflow boundaries on 1D lines, and inflow boundaries for 1D scalar
advection:

- Scalar advection uses `equation.type = "scalar_advection"`, one velocity
  component per dimension, `flux.type = "upwind"`, and scalar `default`,
  `interval`, `circle`, or 1D `sine` regions as appropriate.
- Euler uses `equation.type = "euler"`,
  `heat_capacity_ratio = 1.4`, `flux.type = "hllc"`, and primitive regions
  with positive `density`, a velocity vector, and positive `pressure`.
- Euler limiting is selected with `limiting.enabled = true`. Optional
  `sensor_threshold`, `density_floor`, and `pressure_floor` values default to
  `1e-3`, `1e-10`, and `1e-10`.

CPU 1D scalar advection also supports exactly two ranks with
`parallel.type = "mpi"`. Build and run the MPI variant with:

```bash
cargo build -p wasl-cli --features mpi
mpirun -np 2 target/debug/wasl-cli run path/to/case
mpirun -np 2 target/debug/wasl-cli restart path/to/case \
  path/to/checkpoint-NNNNNN-rank-000000.h5
```

The MPI restart command derives the other rank-local checkpoint path from the
rank-0 filename. Non-periodic boundaries outside serial 1D cases, Euler inflow,
scalar limiting, and Euler with MPI fail explicitly at dispatch until their
numerical drivers are connected.
