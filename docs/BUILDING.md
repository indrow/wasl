# Building

## Dependencies

| Dependency | Version | Purpose |
| --- | --- | --- |
| CMake | 3.25+ | Presets and C++ build |
| Ninja | any | CMake generator |
| C++ compiler | C++20 | Numerical core |
| Rust and Cargo | 2021 edition | CLI and bridge |
| HDF5 (C++) | distribution package | Prepared meshes and checkpoints |
| Open MPI | optional | `cpu-mpi` preset and Cargo `mpi` feature |
| Gmsh | optional | Generating tutorial meshes |
| clang-format, rustfmt | development | Formatting checks |

`docker/Dockerfile` installs everything except Gmsh.

## Native build

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
cargo test --workspace
```

## Docker build

```bash
docker compose build dev
docker compose run --rm dev
```

This is the build used by the initial CI job.

## Formatting

```bash
find cpp crates -type f \( -name '*.hpp' -o -name '*.cpp' \) -print0 \
  | xargs -0 clang-format -i
cargo fmt --all
git diff --check
```

CI checks formatting without rewriting files.

## Presets

| Preset | Backend | Purpose |
| --- | --- | --- |
| `cpu-debug` | Standard C++ | CPU development and CI |
| `cpu-mpi` | Standard C++ and Open MPI | Two-rank halo transport tests |
