# Periodic cone advection

This P2 tutorial transports a circular scalar profile once across the periodic
unit square. From the repository root, generate the mesh with Gmsh, then prepare
it and run the CPU simulation:

```bash
gmsh -2 cases/tutorials/advection/rotating-cone/mesh/geometry.geo \
  -format msh41 -o cases/tutorials/advection/rotating-cone/mesh/source.msh
cargo run --release -p wasl-cli -- mesh prepare \
  cases/tutorials/advection/rotating-cone/mesh/source.msh \
  cases/tutorials/advection/rotating-cone/mesh/prepared.h5
cargo run --release -p wasl-cli -- run cases/tutorials/advection/rotating-cone
```

Verify conservation in the newest run:

```bash
awk -F, 'NR==2 { initial=$4 } NR>2 && ($4-initial > 1e-10 || initial-$4 > 1e-10) { bad=1 } END { if (bad) exit 1; print "mass conserved" }' \
  "$(ls -dt cases/tutorials/advection/rotating-cone/results/run-* | head -1)/history.csv"
```
