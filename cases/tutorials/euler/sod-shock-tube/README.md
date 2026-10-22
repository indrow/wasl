# Sod shock tube

This P2 tutorial solves the Sod shock tube with HLLC fluxes, characteristic
WENO limiting of troubled cells, and positivity preservation. Gas at rest with
density and pressure `1.0` fills `x < 0.5`; gas at rest with density `0.125`
and pressure `0.1` fills `x > 0.5`. Both ends are outflow boundaries. At
`t = 0.2` a rarefaction, a contact, and a shock are inside the domain.

From the repository root, generate the mesh with Gmsh, then prepare it and run
the CPU simulation:

```bash
gmsh -1 cases/tutorials/euler/sod-shock-tube/mesh/geometry.geo \
  -format msh41 -o cases/tutorials/euler/sod-shock-tube/mesh/source.msh
cargo run --release -p wasl-cli -- mesh prepare \
  cases/tutorials/euler/sod-shock-tube/mesh/source.msh \
  cases/tutorials/euler/sod-shock-tube/mesh/prepared.h5
cargo run --release -p wasl-cli -- run cases/tutorials/euler/sod-shock-tube
```

Each VTU snapshot stores the cell-average conserved state as density,
x-momentum, y-momentum, and total energy.

Verify the case without Gmsh:

```bash
cargo build -p wasl-cli
python3 scripts/sod_shock_tube.py target/debug/wasl-cli \
  cases/tutorials/euler/sod-shock-tube
```

The verifier runs uniform 50, 100, and 200 element meshes and compares the
density cell averages with the exact Riemann solution:

```text
elements,density_L1_error,ratio
50,1.083e-02,
100,5.673e-03,1.910
200,2.992e-03,1.896
```

Shocks and contacts limit L1 convergence to first order, so the error should
halve per doubling. The verifier requires a ratio of at least 1.6 per doubling,
an L1 error of at most `4e-3` on 200 elements, positive density and pressure in
every cell, and preservation of a uniform moving state to `1e-12`.
