// Case-level public entry points live in focused translation units:
//   scalar_advection_case.cpp  run_scalar_advection_case, prepare_mpi_run_directory
//   euler_case.cpp             run_euler_case
//
// Shared case helpers live under wasl/case/ (prepared mesh, initial state,
// coefficient layout, periodic boundaries, validation, and the serial run loop).
//
// This file remains as a thin ownership anchor for the case subsystem so the
// historical source path stays stable for agents and
// tooling. It intentionally has no additional symbols.
