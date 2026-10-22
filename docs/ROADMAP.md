# Roadmap

## Completed

- M0: repository conventions, build environments, CI, and initial cases.
- M1: 1D basis, quadrature, affine geometry, Gmsh import, and line mesh.
- M2: deterministic 1D scalar RKDG, boundaries, SSP-RK3, CFL, and
  irregular-mesh convergence.
- M4: straight-sided triangle geometry, topology, and scalar RKDG.
- M5: 2D Euler equations, fluxes, boundaries, and triangle DG.
- M6: troubled-cell sensing, simple and characteristic WENO, positivity.
- M7: TOML validation, prepared HDF5 meshes, checkpoints, CSV and VTU/PVD
  output, Rust CLI, and the coarse `cxx` bridge.
- M8: partition metadata, one ghost layer, host-staged MPI halo exchange, PVTU,
  and rank-local checkpoints.

## M9: Reproducible 1D and 2D tutorials

- [x] Case-level execution for 1D and 2D scalar and Euler solvers, including
  limiting where the case requires it.
- [ ] Add executable 1D tutorials for periodic and inflow advection, smooth
  Euler, the Sod shock tube, restart, and post-processing.
- [ ] Add executable 2D tutorials for rotating and diagonal advection, the
  isentropic vortex, a supersonic wedge, and a partitioned vortex.
- [ ] Give each tutorial a concise README, reproducible mesh and configuration,
  and one verification command using shared scripts where possible.
- [ ] Record expected convergence or conservation and MPI restart behavior in
  automated checks.

Exit: a new user can run and verify at least one 1D and one 2D tutorial without
reading the implementation.

## M10: Straight-sided 3D tetrahedral scalar RKDG

- [ ] Add compile-time dimension primitives for vectors,
  normals, gradients, fluxes, and dimension-dependent equation state.
- [ ] Add shared simplex traits for vertex, face, mode, permutation, and
  orientation rules while retaining specialized line, triangle, and
  tetrahedron kernels.
- [ ] Add an orthonormal tetrahedral basis, volume and triangular-face
  quadrature, affine geometry, and face restrictions for P0-P3.
- [ ] Import Gmsh tetrahedra and build canonical oriented triangular faces,
  physical groups, partition faces, and one ghost layer.
- [ ] Implement the deterministic CPU 3D scalar-advection residual, boundary
  conditions, SSP-RK3, and CFL oracle.
- [ ] Reuse the M8 halo exchange for partitioned 3D execution.
- [ ] Verify free streams, distorted meshes, partition equivalence, and P1/P2
  convergence on irregular tetrahedra.

Exit: smooth 3D scalar advection converges at `p + 1`, with serial and
partitioned agreement at documented tolerances.

## M11: Three-dimensional Euler

- [ ] Add the five-component conserved state, three-direction physical flux,
  and arbitrary-normal Rusanov and HLLC fluxes.
- [ ] Add slip-wall, periodic, far-field, and checked inflow/outflow boundaries.
- [ ] Implement tetrahedral Euler residuals with global conserved
  diagnostics.
- [ ] Verify a distorted-mesh free stream, embedded Sod tube, smooth acoustic
  convergence, and coordinate-permutation symmetry.

Exit: smooth P1/P2 convergence and shock-tube behavior agree across serial
and partitioned execution without dependence on face ordering.

## M12: Three-dimensional DG-WENO and positivity

- [ ] Extend modal sensing, neighbor transforms, characteristic limiting, and
  positivity preservation to tetrahedra.
- [ ] Keep stencil construction deterministic across partitions and preserve
  every troubled cell's mode-zero coefficients.
- [ ] Verify an embedded shock tube, spherical blast, partition-crossing shock,
  and near-vacuum positivity case.

Exit: 3D discontinuous Euler cases remain conservative and physical across
partition boundaries at the documented floors and tolerances.

## M13: Three-dimensional scalability and performance

- [ ] Measure and document state, geometry, and workspace memory by polynomial
  order before changing their layouts.
- [ ] Split local and partition-boundary work, then overlap useful computation
  with halo exchange where profiling shows a benefit.
- [ ] Add scalable global CFL and diagnostic reductions for more than two
  ranks.
- [ ] Record strong and weak scaling, checkpoint throughput, and stable
  performance baselines.

Exit: the complete 3D path has reproducible performance baselines and measured
scaling, with no loss of numerical parity.

## M14: Curved and high-order geometry

- [ ] Import high-order Gmsh line, triangle, and tetrahedral elements with
  geometry nodes distinct from solution coefficients.
- [ ] Add isoparametric maps, per-quadrature-point metrics, curved face normals,
  and checked Jacobians.
- [ ] Preserve discrete metric identities and free streams on curved meshes.
- [ ] Verify curved periodic advection and order recovery when geometry and
  solution orders match.

Exit: curved domains recover the expected convergence order without losing
free-stream preservation.

## Deferred

P4+ solution spaces, p-adaptation, mixed elements, AMR, Navier-Stokes, MHD,
reacting flow, dynamic plugins, and Python bindings.
