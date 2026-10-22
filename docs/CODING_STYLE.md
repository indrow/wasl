# Coding Style

This guide defines how Wasl code should look and behave. Automated formatters
settle whitespace. This document covers choices that formatters cannot decide.

## General principles

- Prefer the smallest implementation that satisfies the current milestone.
- Keep ownership visible and lifetimes local.
- Put code in the subsystem that owns its invariant.
- Use standard-library facilities before adding dependencies or helpers.
- Do not introduce an interface with one implementation, a factory for one
  product, or configuration for a value that is not variable.
- Avoid generic `utils`, `common`, `manager`, and `misc` modules.
- Validate input at file, process, FFI, and user-command boundaries.
- Do not silently recover from invalid numerical state or unsupported input.
- Do not use em dashes or emoji in repository text or user-facing output.

## Automated formatting

C++ uses the root `.clang-format` file:

```bash
find cpp -type f \( -name '*.hpp' -o -name '*.cpp' \) -print0 \
  | xargs -0 clang-format -i
```

Rust uses `rustfmt.toml`:

```bash
cargo fmt --all
```

CI checks both formatters and `git diff --check`. Do not hand-align code in a
way the formatter will immediately undo.

## Naming

### Files and directories

- C++ and Rust files use lowercase `snake_case`.
- C++ tests use `<subject>_test.cpp`.
- Case directories and shell entry points may use lowercase hyphenated names.
- Public include paths mirror subsystem ownership, for example
  `wasl/mesh/line_mesh.hpp`.

### C++

- Namespaces, functions, and variables use `snake_case`.
- Types and template parameters use `PascalCase`.
- Public namespace constants use `snake_case`. Internal and test constants may
  use the `kPascalCase` form.
- Private data members use a trailing underscore when a class needs to
  distinguish them from parameters.
- Avoid abbreviations unless they are standard numerical terms such as `dg`,
  `rk`, `cfl`, or `weno`.

### Rust

Follow standard Rust naming enforced by the compiler and Clippy:

- modules, functions, and variables use `snake_case`;
- types and traits use `PascalCase`;
- constants use `SCREAMING_SNAKE_CASE`.

## C++ headers and sources

- Headers use `#pragma once`.
- Include every standard and project header used directly. Do not depend on
  transitive includes.
- Order includes as: matching project header, standard headers, then other
  project headers. Let clang-format sort within each group.
- Public headers contain declarations, small value types, templates, and
  header-only numerical evaluators.
- Host parsing, file I/O, mesh preparation, MPI orchestration, and other
  non-template implementations belong in `cpp/src`.
- Keep implementation details in an unnamed namespace inside `.cpp` files. Use
  a named `detail` namespace only when templates require header visibility.
- Do not expose MPI, HDF5, or `cxx` implementation types through core
  numerical APIs.

## Types and ownership

- Use RAII for owned resources.
- Use values by default. Use `std::unique_ptr` only for dynamic ownership and
  references or pointers for non-owning access.
- A raw pointer in a numerical view is non-owning and must have an obvious owner.
- Use fixed-width integers for persisted formats and cross-language APIs.
- Preserve external Gmsh tags in 64-bit types. Use 32-bit local indices only
  after validating the range.
- Prefer `std::size_t` for container sizes and indexing. Convert explicitly at
  API boundaries.
- Avoid shared ownership unless the lifetime genuinely has multiple owners.

## Functions and classes

- A function should perform one named operation and expose required
  preconditions.
- Prefer free functions and small structs over classes that only group static
  functions.
- Mark functions `const`, `noexcept`, and `constexpr` only when their contracts
  genuinely allow it.
- Do not mark fallible allocation, parsing, or validation paths `noexcept`.
- Use early returns for invalid or completed paths when they reduce nesting.
- Avoid boolean parameters whose meaning is unclear at the call site. Use a
  named enum when two modes are both real requirements.
- Keep public APIs narrow. Add flexibility only when a second concrete use
  requires it.

## Error handling

- Throw exceptions from host-side C++ when an operation cannot return a valid
  result.
- Error messages at trust boundaries must identify the operation and field,
  such as `MSH reader: invalid or missing node tag`.
- Never dereference a failed search or use an unchecked narrowing conversion.
- Reject unsupported file versions, element types, and configurations
  explicitly.
- Rust converts bridge failures into user-facing command errors. It does not
  terminate for normal invalid input.

## Numerical code

- State mathematical assumptions in the owning header or test.
- Follow `NUMERICAL_CONVENTIONS.md` exactly for basis normalization, ordering,
  face orientation, state layout, and precision.
- Use descriptive numerical names such as `jacobian`, `mode_count`, and
  `outward_normal`. Avoid unexplained single-letter variables outside short
  formulas and loops.
- Document tolerances beside the property they test.
- Do not compare floating-point values exactly unless the result is exact by
  construction.
- Preserve constant states, conservation, and cell averages by construction,
  then test those properties.
- Keep the CPU path simple and deterministic. It is the correctness oracle, not
  a second independently designed solver.
- Mark a deliberate algorithmic ceiling with a `ponytail:` comment that names
  the upgrade condition.

Example:

```cpp
// ponytail: linear lookup is sufficient below 100 entries; index by tag when
// mesh metadata grows beyond that limit.
```

## Numerical verification

Compilation is not numerical validation. Match every change to the cheapest
property that can expose its likely failure:

- basis/quadrature: polynomial exactness;
- mappings: round trip and constant metric checks;
- faces: deterministic orientation and flux cancellation;
- DG residual: constant-state preservation;
- time integration: known ODE order;
- smooth PDE path: measured convergence;
- limiter: cell-average preservation and reduced oscillation;
- Euler: positive supported states and conservation;
- partitioning: agreement with the single-rank result.

Convergence records must include mesh sizes, polynomial order, timestep policy,
norms, observed rates, precision, and tolerance. Reduce the timestep enough to
separate spatial error from SSP-RK3 temporal error.

## Rust and the bridge

- Rust owns CLI parsing, terminal presentation, and coarse command flow.
- C++ owns complete simulations and all numerical state.
- Bridge functions operate at case level, for example `check_case` and
  `run_case`. Never expose numerical routines or solver-array mutation through `cxx`.
- Keep `unsafe` code confined to the generated or unavoidable FFI boundary and
  explain its safety contract.
- Avoid duplicating C++ configuration and numerical types in Rust when an
  opaque handle or serialized result is sufficient.

## CMake

- Use target-based commands such as `target_link_libraries` and
  `target_compile_options`.
- Do not set global compiler flags when a target property is sufficient.
- Dependencies are required only by the target that uses them.
- Default builds must not require MPI or Gmsh.
- Name CTest tests under the `wasl.` namespace.
- Add a source file to the narrowest owning target. Do not collect unrelated
  sources into a generic convenience library.

## Shell and automation

- Shell scripts use Bash with `set -euo pipefail`.
- Quote paths and variable expansions.
- Prefer existing tools and direct commands over wrapper scripts.
- Scripts must be non-interactive unless interaction is their explicit purpose.
- Generated artifacts go under ignored build, target, case result, or temporary
  directories.

## Version control

- Treat existing tracked and untracked changes as user work. Do not discard,
  overwrite, or reformat unrelated files.
- Inspect `git status --short` before editing and before reporting completion.
- An agent may create a local commit without a separate request when a coherent
  task checkpoint is complete and relevant checks pass.
- Before committing, inspect the final diff and status, then stage only files
  owned by that change. Do not commit partial or failing work merely to clean
  the working tree.
- Use an imperative commit subject in the form `<scope>: <summary>`, for example
  `mesh: validate MSH node references`.
- Record meaningful validation commands in the commit body when a change is
  numerical, cross-language, or MPI-specific.
- Commit source, tests, configuration, and required documentation together so
  repository state does not describe code that is absent.
- Do not commit build directories, generated meshes, result directories,
  checkpoints, credentials, local editor state, or machine-specific paths.
- Do not commit large numerical outputs when compact metrics can reproduce the
  regression expectation.
- Do not amend, rebase, force-push, delete branches, create tags, or publish
  changes unless the user explicitly requests that operation.
- Do not push merely because a task is complete. Publishing requires an
  explicit user request.
- Never use destructive commands such as `git reset --hard` or discard another
  contributor's changes to make a clean working tree.
- Before handing off, run `git diff --check` and report relevant untracked or
  uncommitted work without claiming it is committed.

## Tests

- Put subsystem tests under the matching `cpp/tests` directory.
- Put cross-subsystem workflows under root `tests`.
- Use the smallest test that fails when the property breaks.
- Tests must report failure through a nonzero exit status.
- Numerical tests should check properties, not only stored output snapshots.
- A parser test includes malformed input and verifies a descriptive error.
- Avoid test frameworks until plain executables become measurably cumbersome.
- Keep fixtures small, local, deterministic, and automatically cleaned up.

## Comments and documentation

- Keep comments brief. Explain why, invariants, ownership, units, and
  mathematical conventions.
- Do not narrate syntax that the code already states.
- Keep comments synchronized with behavior or delete them.
- Use `TODO:` for scoped missing work, `FIXME:` for known incorrect behavior,
  `HACK:` for a deliberate workaround, and `NOTE:` for a non-obvious
  constraint. State the action or reason, the replacement condition for a
  workaround, and an issue reference when one exists. Do not leave vague or
  stale annotations.
- Public user workflows belong in `README.md` or focused documents under
  `docs`.
- Current status belongs only in `PROJECT_STATE.md`.
- `PROJECT_STATE.md` retains current reality, and `ROADMAP.md` retains current
  and future milestones. Summarize completed milestones in one line each.
- Durable architectural decisions are recorded in `ARCHITECTURE.md`.

## Before completing a change

Run:

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
cargo fmt --all --check
cargo test --workspace
git diff --check
```

Also run the relevant convergence, MPI, or file-format check when the
change affects that subsystem. Clearly report checks that were not run.
