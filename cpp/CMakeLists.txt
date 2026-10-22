option(WASL_ENABLE_MPI "Enable the MPI halo-exchange transport and mpirun tests" OFF)

find_package(HDF5 REQUIRED COMPONENTS CXX)

if(WASL_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

set(WASL_MESH_SOURCES
  src/mesh/line/build_line_mesh.cpp
  src/mesh/line/line_mesh.cpp
  src/mesh/msh/msh_document.cpp
  src/mesh/msh/msh_node_lookup.cpp
  src/mesh/msh/msh_reader.cpp
  src/mesh/msh/msh_sections.cpp
  src/mesh/partition/partition.cpp
  src/mesh/partition/rank_local_line_mesh.cpp
  src/mesh/prepared/mesh_from_prepared.cpp
  src/mesh/prepared/owned_mesh_slice.cpp
  src/mesh/prepared/prepare_mesh.cpp
  src/mesh/prepared/prepared_mesh_validation.cpp
  src/mesh/prepared/prepared_node_table.cpp
  src/mesh/prepared/prepared_physical_names.cpp
  src/mesh/triangle/boundary_edge_tags.cpp
  src/mesh/triangle/build_triangle_mesh.cpp
  src/mesh/triangle/triangle_mesh.cpp
)
set(WASL_IO_SOURCES
  src/io/checkpoint.cpp
  src/io/checkpoint_file.cpp
  src/io/checkpoint_reassembly.cpp
  src/io/history_csv.cpp
  src/io/output_names.cpp
  src/io/prepared_mesh_file.cpp
  src/io/pvd_writer.cpp
  src/io/pvtu_writer.cpp
  src/io/text_output.cpp
  src/io/vtu_writer.cpp
)
set(WASL_BOUNDARIES_SOURCES
  src/boundaries/line_euler_exterior.cpp
  src/boundaries/periodic_triangle_trace.cpp
)
set(WASL_TIME_SOURCES
  src/time/cfl_timestep.cpp
)
set(WASL_LIMITING_SOURCES
  src/limiting/euler_limiter_stencil_line.cpp
  src/limiting/euler_limiter_stencil_triangle.cpp
  src/limiting/euler_modal_limiter.cpp
)
set(WASL_COMMUNICATION_SOURCES
  src/communication/halo_schedule.cpp
  src/communication/in_process_halo_exchange.cpp
  src/communication/mpi_session.cpp
)
set(WASL_CASE_SOURCES
  src/case/case_mesh.cpp
  src/case/case_validation.cpp
  src/case/euler_case.cpp
  src/case/euler_diagnostics.cpp
  src/case/euler_line.cpp
  src/case/euler_state_operations.cpp
  src/case/euler_triangle.cpp
  src/case/initial_regions.cpp
  src/case/initial_state.cpp
  src/case/line_boundary_resolution.cpp
  src/case/mpi_run_directory.cpp
  src/case/scalar_advection_case.cpp
  src/case/scalar_advection_distributed.cpp
  src/case/scalar_advection_line.cpp
  src/case/scalar_advection_triangle.cpp
  src/case/scalar_coefficient_layout.cpp
  src/case/scalar_diagnostics.cpp
  src/case/scalar_state_operations.cpp
  src/case/serial_case_runner.cpp
)

add_library(wasl_core STATIC
  ${WASL_MESH_SOURCES}
  ${WASL_IO_SOURCES}
  ${WASL_BOUNDARIES_SOURCES}
  ${WASL_TIME_SOURCES}
  ${WASL_LIMITING_SOURCES}
  ${WASL_COMMUNICATION_SOURCES}
  ${WASL_CASE_SOURCES}
)
# Public headers live under include/wasl; private headers sit beside their
# sources and are included relative to src.
target_include_directories(wasl_core PUBLIC include PRIVATE src)
target_link_libraries(wasl_core PUBLIC HDF5::HDF5)
if(WASL_ENABLE_MPI)
  target_compile_definitions(wasl_core PRIVATE WASL_ENABLE_MPI OMPI_SKIP_MPICXX)
  target_link_libraries(wasl_core PUBLIC MPI::MPI_CXX)
endif()
target_compile_features(wasl_core PUBLIC cxx_std_20)
target_compile_options(wasl_core PUBLIC
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wall;-Wextra;-Wpedantic>
)

add_executable(wasl-smoke src/smoke.cpp)
target_link_libraries(wasl-smoke PRIVATE wasl_core)

if(BUILD_TESTING)
  # Each test is one executable named after its source file; extra arguments
  # are passed to it on the command line.
  function(wasl_add_test name source)
    get_filename_component(stem "${source}" NAME_WE)
    set(target "wasl-${stem}")
    string(REPLACE "_" "-" target "${target}")
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE wasl_core)
    add_test(NAME ${name} COMMAND ${target} ${ARGN})
  endfunction()

  # M1 numerical kernels: basis, quadrature, geometry, and mesh.
  wasl_add_test(wasl.core tests/basis/mode_count_test.cpp)
  wasl_add_test(wasl.legendre tests/basis/legendre_test.cpp)
  wasl_add_test(wasl.quadrature tests/quadrature/quadrature_test.cpp)
  wasl_add_test(wasl.line-geometry tests/mesh/line_geometry_test.cpp)
  wasl_add_test(wasl.line-mesh tests/mesh/line_mesh_test.cpp)

  # M4 straight-sided triangle basis, quadrature, geometry, and topology.
  wasl_add_test(wasl.basis.dubiner-triangle tests/basis/dubiner_triangle_test.cpp)
  wasl_add_test(wasl.quadrature.triangle tests/quadrature/triangle_test.cpp)
  wasl_add_test(wasl.triangle-geometry tests/mesh/triangle_geometry_test.cpp)
  wasl_add_test(wasl.triangle-mesh tests/mesh/triangle_mesh_test.cpp)

  # M8 partition ownership, ghost layer, and the MPI-free halo oracle.
  wasl_add_test(wasl.mesh.partition tests/mesh/partition_test.cpp)
  wasl_add_test(wasl.mpi.halo tests/communication/halo_exchange_test.cpp)

  # M2/M5 equations and DG residuals.
  wasl_add_test(wasl.scalar-advection tests/equations/scalar_advection_test.cpp)
  wasl_add_test(wasl.euler tests/equations/euler_test.cpp)
  wasl_add_test(wasl.dg.residual tests/dg/residual_test.cpp)
  wasl_add_test(wasl.dg.triangle-residual tests/dg/triangle_residual_test.cpp)
  wasl_add_test(wasl.dg.euler-triangle-residual tests/dg/euler_triangle_residual_test.cpp)
  wasl_add_test(wasl.dg.euler-line-residual tests/dg/euler_line_residual_test.cpp)

  # M6 limiting.
  wasl_add_test(wasl.limiter.modal-decay tests/limiter/modal_decay_test.cpp)
  wasl_add_test(wasl.limiter.simple-weno tests/limiter/simple_weno_test.cpp)
  wasl_add_test(wasl.limiter.euler tests/limiter/euler_limiting_test.cpp)
  wasl_add_test(wasl.limiter.line-stencil tests/limiter/euler_limiter_stencil_test.cpp)

  # M2/M7 time integration, boundaries, I/O, and case runs.
  wasl_add_test(wasl.time.ssp-rk3 tests/time/ssp_rk3_test.cpp)
  wasl_add_test(wasl.boundaries tests/boundaries/boundaries_test.cpp)
  wasl_add_test(wasl.boundaries.line-euler tests/boundaries/line_euler_exterior_test.cpp)
  wasl_add_test(wasl.io.case tests/io/case_io_test.cpp)
  wasl_add_test(wasl.runtime.euler-case tests/runtime/euler_case_test.cpp)
  wasl_add_test(wasl.runtime.scalar-boundary-case tests/runtime/scalar_boundary_case_test.cpp)
  wasl_add_test(wasl.runtime.euler-boundary-case tests/runtime/euler_boundary_case_test.cpp)
  wasl_add_test(wasl.runtime.case-regression tests/runtime/case_regression_test.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/runtime/case_regression_reference.txt)

  if(WASL_ENABLE_MPI)
    # Two ranks reproduce the serial residual and case loop on a split mesh;
    # wasl.mpi.halo already proves the exchange semantics without MPI.
    add_executable(wasl-halo-mpi-test tests/communication/halo_mpi_test.cpp)
    target_link_libraries(wasl-halo-mpi-test PRIVATE wasl_core MPI::MPI_CXX)
    add_test(NAME wasl.mpi.transport
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:wasl-halo-mpi-test> ${MPIEXEC_POSTFLAGS})
  endif()
endif()
