# Contributing to Wasl

Thanks for helping improve Wasl. This guide covers setup, the change workflow,
and what a reviewer will check.

## Set up

You need CMake 3.25+, Ninja, a C++20 compiler, Rust (2021 edition), and the
HDF5 C++ development package. Open MPI and Gmsh are optional. See
[`docs/BUILDING.md`](docs/BUILDING.md) for the full list.

The quickest reproducible setup is the development container:

```bash
docker compose build
docker compose run --rm dev
```

## Build and test

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
cargo test --workspace
```

If your change touches partitioning, halo exchange, or MPI dispatch, also run:

```bash
cmake --preset cpu-mpi
cmake --build --preset cpu-mpi
ctest --preset cpu-mpi
cargo build -p wasl-cli --features mpi
```

## Format

CI rejects unformatted code. Before committing:

```bash
find cpp crates -type f \( -name '*.hpp' -o -name '*.cpp' \) -print0 \
  | xargs -0 clang-format -i
cargo fmt --all
git diff --check
```

## Making a change

1. Read [`docs/PROJECT_STATE.md`](docs/PROJECT_STATE.md) and
   [`docs/ROADMAP.md`](docs/ROADMAP.md) to see what exists and what comes next.
2. Keep the change focused on one concern.
3. Put code in the subsystem that owns its invariant, under
   `cpp/include/wasl/<subsystem>` and `cpp/src/<subsystem>`, with a matching
   test in `cpp/tests/<subsystem>` registered in `cpp/CMakeLists.txt`.
4. Follow [`docs/CODING_STYLE.md`](docs/CODING_STYLE.md) for language, error
   handling, and testing rules.
5. Update `docs/PROJECT_STATE.md` when capability changes, and
   `docs/ROADMAP.md` when a milestone item is completed.

## Numerical changes

A numerical change needs the smallest test that exposes its likely failure:
conservation, consistency, constant-state preservation, or measured
convergence. The verification checklist in `docs/CODING_STYLE.md` lists the
expected check for each subsystem.

- Do not change bases, orientation, or state layout without updating
  [`docs/NUMERICAL_CONVENTIONS.md`](docs/NUMERICAL_CONVENTIONS.md).
- Record durable design decisions in
  [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
- Keep the CPU path deterministic: results must not depend on hash-map order,
  thread scheduling, or the standard library implementation.

## Rust and the bridge

Rust owns the CLI and case validation. It calls C++ only for coarse,
case-level operations and never holds solver arrays. New case options must be
validated in `crates/wasl-bridge/src/config.rs` and documented in
[`docs/CASE_FORMAT.md`](docs/CASE_FORMAT.md).

## Commits and pull requests

- Use an imperative subject in the form `<scope>: <summary>`, for example
  `mesh: validate MSH node references`.
- List the validation commands you ran in the commit body for numerical,
  cross-language, or MPI changes.
- Commit source, tests, and required documentation together.
- Do not use em dashes or emoji in code, comments, docs, or commit messages.

Before opening a pull request, confirm that:

- [ ] `ctest --preset cpu-debug` and `cargo test --workspace` pass;
- [ ] the `cpu-mpi` checks pass if MPI code changed;
- [ ] formatting checks pass;
- [ ] new behavior has a focused test;
- [ ] docs reflect the change.

## License

By contributing, you agree that your contributions are licensed under the
terms in [`LICENSE`](LICENSE).
