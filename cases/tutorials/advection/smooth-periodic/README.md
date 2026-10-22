# Smooth periodic advection

This CPU P2 tutorial transports `sin(2 pi x)` on the periodic unit interval.
From the repository root, build the CLI and run the shared verifier:

```bash
cargo build -p wasl-cli
python3 scripts/convergence.py target/debug/wasl-cli \
  cases/tutorials/advection/smooth-periodic
```

The command generates uniform 8, 16, and 32 element MSH 4.1 meshes in a
temporary directory, runs the case at each resolution, prints cell-average L2
errors and observed rates, and requires every P2 rate to be at least 2.7.
