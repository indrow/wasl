# Agent Instructions

Wasl is a framework for high-order modal RKDG with WENO limiting on
unstructured Gmsh meshes, running on the CPU.

## Read first

1. `docs/PROJECT_STATE.md`
2. `docs/NUMERICAL_CONVENTIONS.md`
3. `docs/ARCHITECTURE.md`
4. `docs/CODING_STYLE.md`
5. `docs/ROADMAP.md`

Inspect code and tests after reading; they override stale prose.

## Build and check

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
cargo test --workspace
git diff --check
```


## Non-negotiable rules

- C++20 numerical core running on the CPU.
- Rust calls coarse case operations, never numerical routines or solver arrays.
- Use modal RKDG; apply WENO only to troubled elements.
- Start with unstructured 1D line elements, then straight-sided 2D triangles.
- Preserve a deterministic CPU reference path.
- Standard DG and WENO need one ghost-element layer per MPI rank.
- No global sparse matrix or MPI/HDF5 types in public numerical APIs.
- A numerical change needs the smallest test that checks its conservation,
  consistency, or convergence risk.
- Do not change numerical conventions without updating the convention file
  and recording the decision in `docs/ARCHITECTURE.md`.
- Do not scaffold speculative modules. Add files when the milestone needs them.
- Put a new file in the narrowest existing `cpp/include/wasl/<subsystem>`
  directory that owns its invariant, with a matching `cpp/tests` test.
- Follow `docs/CODING_STYLE.md` for language, error, test, and automation rules.
- Keep code comments brief and explain why or an invariant, not obvious syntax.
- Spell debt annotations `TODO:`, `FIXME:`, `HACK:`, or `NOTE:` and give a
  concrete action or reason, plus an issue reference when one exists.
- Keep `docs/PROJECT_STATE.md` and `docs/ROADMAP.md` describing current
  capability; summarize completed milestones in one line each.
- An agent may create a local commit at a coherent, validated task checkpoint.
  Inspect the diff, stage only task-owned files, and report the hash and checks.
- Preserve unrelated working-tree changes. Do not amend, rebase, tag, or push
  unless the user explicitly requests it.

## Writing style

- Do not use em dashes or emoji in repository text, code comments, commits, or
  generated user-facing output.
- Use commas, colons, parentheses, or short sentences instead.

## Before finishing a task

1. Run relevant formatting, builds, and tests.
2. Update `docs/PROJECT_STATE.md` when reality changed.
3. Update `docs/ROADMAP.md` only when milestone status changed.
4. Record durable architectural decisions in `docs/ARCHITECTURE.md`.
5. State checks not run and one concrete next task with acceptance criteria.
