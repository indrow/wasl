# Architecture

## Scope

Wasl solves conservation laws using modal RKDG on unstructured meshes. WENO
reconstructs only troubled elements. The first usable path is scalar advection
on an unstructured 1D line mesh, followed by straight-sided 2D triangles.
Euler and distributed execution build on proven mesh and polynomial
infrastructure.

## Ownership boundary

```text
Rust CLI
  check/run/restart/inspect/export
                 |
                 v
C++ case engine and orchestration
                 |
                 v
     deterministic CPU runner
```

Rust never calls individual numerical routines or owns numerical arrays. C++
executes a complete operation across the bridge, parses MSH, controls
HDF5/VTK output, and runs all numerical work.

## Runtime data flow

```text
geometry.geo -> Gmsh -> source.msh -> wasl mesh prepare -> prepared.h5
                                                        -> mesh/state
                                                        -> RKDG stages
                                                        -> HDF5/VTK/CSV
```

The MSH parser runs once during preparation. `prepared.h5` contains normalized
IDs, element/face connectivity, orientations, physical groups, metrics,
partition ownership, and one ghost-element layer. The first element type is a
1D line; the next is a straight-sided 2D triangle.

## Numerical composition

Residuals use compile-time equation, flux, basis/order, and scalar policies.
Configuration selects only combinations compiled into the binary. Runtime
polymorphism stays in case orchestration, not in inner numerical loops.

Every SSP-RK3 stage eventually performs:

1. Exchange ghost modal coefficients.
2. Evaluate local volume terms and local faces.
3. Complete exchange and evaluate partition faces.
4. Apply mass inverse/lifting and update the stage.
5. Detect and compact troubled elements.
6. Limit only the compacted elements and apply positivity correction.
7. Reduce timestep and diagnostics.

No global sparse matrix is assembled.

## Source layout

The C++ core is organized one concept per file, grouped by the object it
describes:

```text
basis, quadrature        reference-element polynomials and rules
mesh/{msh,line,triangle} MSH import and canonical topologies
mesh/{prepared,partition} file-ready mesh, ownership, and ghosts
equations, flux          physics and numerical fluxes
dg, boundaries, time     residuals, exterior traces, SSP-RK3 and CFL
limiting                 sensors, WENO, positivity, Euler limiter stencil
communication            halo schedule, in-process oracle, MPI transport
io                       one file per format (HDF5, CSV, VTU, PVTU, PVD)
case                     configuration, validation, and case drivers
```

Dependencies point down this list. `mesh` never depends on `io`; `io` reads and
writes the `mesh` representation.

## Distributed execution

One MPI rank owns one mesh partition and its elements. Host-staged nonblocking
MPI is the transport. Work overlaps naturally by evaluating local volume and face terms while halo data
is in flight.

## Design decisions

- Modal RKDG is the primary discretization. WENO limits only troubled elements,
  using immediate face neighbors and preserving cell averages, so one
  ghost-element layer is enough.
- Gmsh MSH 4.1 is source input only. `wasl mesh prepare` normalizes IDs, builds
  and orients faces, maps physical groups, computes affine metrics, and writes
  `prepared.h5`; runs consume that file.
- One MPI rank owns one mesh partition. Elements are assigned by contiguous
  block, `rank(e) = floor(e * rank_count / element_count)`, and a face is owned
  by its left element's owner. `prepared.h5` (format version 2) stores
  `partition_rank_count`, `element_owners`, and `face_owners`; the ghost layer
  is derived from ownership and topology.
- Halo exchange runs behind a transport-agnostic schedule. A rank-local array
  holds owned elements followed by ghosts, so the residual is unchanged at
  partition faces. An in-process oracle proves the semantics without MPI; the
  host-staged nonblocking MPI transport drives the same schedule.

- Case configurations share `CaseFiles` (mesh, output, restart paths) and
  `RunControl` (time range, step limit, output intervals). A serial case
  describes its state layout, diagnostics, time-step limit, and advance as a
  `SerialCaseModel`, and one runner owns restart, stepping, history, snapshots,
  and checkpoint rotation for every equation and element type.
- Case boundary conditions are keyed by mesh physical name. The case resolves
  each boundary face's physical group to its condition before the run, and
  rejects unconfigured faces, unmatched names, and half-periodic pairs. An empty
  boundary list means periodic, so periodic cases need no names.
- Line limiter stencils join the two boundary elements only for periodic ends.
  Open (outflow) ends leave each boundary element one candidate, so WENO never
  blends a state with the far side of the domain.
- The structure-of-arrays index of `NUMERICAL_CONVENTIONS.md` is spelled in one
  place, `ModalLayout`. The scalar solver path keeps `U[element][mode]` and
  converts at checkpoint boundaries.

## File formats

- Input mesh: Gmsh MSH 4.1, prepared once.
- Runtime mesh and checkpoints: HDF5.
- Visualization: VTU, PVTU, and PVD.
- Histories: CSV.
- Configuration and effective run manifest: TOML.

Each invocation creates a new run directory and never silently overwrites an
older run.
