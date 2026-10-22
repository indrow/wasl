# Inflow advection

This P2 tutorial advects a constant inflow state into an empty unit interval.
The left end prescribes `1.0`, the right end is a zero-gradient outflow, and the
wave speed is `1.0`. The front enters at `t = 0`, leaves through the outflow end
at `t = 1`, and the run ends at `t = 1.5`, when the exact state is `1.0`
everywhere.

From the repository root, generate the mesh with Gmsh, then prepare it and run
the CPU simulation:

```bash
gmsh -1 cases/tutorials/advection/line-inflow/mesh/geometry.geo \
  -format msh41 -o cases/tutorials/advection/line-inflow/mesh/source.msh
cargo run --release -p wasl-cli -- mesh prepare \
  cases/tutorials/advection/line-inflow/mesh/source.msh \
  cases/tutorials/advection/line-inflow/mesh/prepared.h5
cargo run --release -p wasl-cli -- run cases/tutorials/advection/line-inflow
```

Verify the case without Gmsh:

```bash
cargo build -p wasl-cli
python3 scripts/inflow_advection.py target/debug/wasl-cli \
  cases/tutorials/advection/line-inflow
```

The verifier writes a uniform 40-element mesh in a temporary directory and
checks three properties to `1e-12`:

- `balance`: while the front is at least 0.8 of the domain from the outflow
  end, the total mass equals the inflow balance `t`. This requires the inflow
  value at every SSP-RK3 stage.
- `transport`: at `t = 1.5` every cell average equals the inflow state.
- `free_stream`: a run starting from the inflow state never changes.

On coarser meshes the front's numerical tail decays more slowly: 20 elements
leave a `3e-9` transport deviation at `t = 1.5`.
