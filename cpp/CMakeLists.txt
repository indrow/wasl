option(WASL_ENABLE_MPI "Enable the MPI halo-exchange transport and mpirun tests" OFF)

find_package(HDF5 REQUIRED COMPONENTS CXX)

if(WASL_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

add_library(wasl_core STATIC
  src/case/coefficient_layout.cpp
  src/case/euler_case.cpp
  src/case/initial_state.cpp
  src/case/periodic_boundary.cpp
  src/case/prepared_mesh.cpp
  src/case/scalar_advection_case.cpp
  src/case/simulation.cpp
  src/case/validation.cpp
  src/mesh/line_mesh.cpp
  src/mesh/msh_reader.cpp
  src/mesh/partition.cpp
  src/mesh/triangle_mesh.cpp
  src/communication/halo_exchange.cpp
  src/io/case_io.cpp
)
target_include_directories(wasl_core PUBLIC include)
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
  # Core mode-count smoke test.
  add_executable(wasl-core-test tests/basis/reference_element_test.cpp)
  target_link_libraries(wasl-core-test PRIVATE wasl_core)
  add_test(NAME wasl.core COMMAND wasl-core-test)

  # M1 numerical-kernel tests: basis, quadrature, geometry, and mesh.
  add_executable(wasl-legendre-test tests/basis/legendre_test.cpp)
  target_link_libraries(wasl-legendre-test PRIVATE wasl_core)
  add_test(NAME wasl.legendre COMMAND wasl-legendre-test)

  add_executable(wasl-quadrature-test tests/quadrature/quadrature_test.cpp)
  target_link_libraries(wasl-quadrature-test PRIVATE wasl_core)
  add_test(NAME wasl.quadrature COMMAND wasl-quadrature-test)

  add_executable(wasl-line-geometry-test tests/mesh/line_geometry_test.cpp)
  target_link_libraries(wasl-line-geometry-test PRIVATE wasl_core)
  add_test(NAME wasl.line-geometry COMMAND wasl-line-geometry-test)

  add_executable(wasl-line-mesh-test tests/mesh/line_mesh_test.cpp)
  target_link_libraries(wasl-line-mesh-test PRIVATE wasl_core)
  add_test(NAME wasl.line-mesh COMMAND wasl-line-mesh-test)

  # M4 straight-sided 2D triangle basis, quadrature, and geometry.
  add_executable(wasl-dubiner-triangle-test tests/basis/dubiner_triangle_test.cpp)
  target_link_libraries(wasl-dubiner-triangle-test PRIVATE wasl_core)
  add_test(NAME wasl.basis.dubiner-triangle COMMAND wasl-dubiner-triangle-test)

  add_executable(wasl-triangle-test tests/quadrature/triangle_test.cpp)
  target_link_libraries(wasl-triangle-test PRIVATE wasl_core)
  add_test(NAME wasl.quadrature.triangle COMMAND wasl-triangle-test)

  add_executable(wasl-triangle-geometry-test tests/mesh/triangle_geometry_test.cpp)
  target_link_libraries(wasl-triangle-geometry-test PRIVATE wasl_core)
  add_test(NAME wasl.triangle-geometry COMMAND wasl-triangle-geometry-test)

  add_executable(wasl-triangle-mesh-test tests/mesh/triangle_mesh_test.cpp)
  target_link_libraries(wasl-triangle-mesh-test PRIVATE wasl_core)
  add_test(NAME wasl.triangle-mesh COMMAND wasl-triangle-mesh-test)

  # M8 deterministic partition ownership and one ghost-element layer.
  add_executable(wasl-partition-test tests/mesh/partition_test.cpp)
  target_link_libraries(wasl-partition-test PRIVATE wasl_core)
  add_test(NAME wasl.mesh.partition COMMAND wasl-partition-test)

  # M8 host-staged halo exchange oracle (MPI-free; single process).
  add_executable(wasl-halo-exchange-test tests/communication/halo_exchange_test.cpp)
  target_link_libraries(wasl-halo-exchange-test PRIVATE wasl_core)
  add_test(NAME wasl.mpi.halo COMMAND wasl-halo-exchange-test)

  # M2 scalar-advection equation and DG residual tests.
  add_executable(wasl-scalar-advection-test tests/equations/scalar_advection_test.cpp)
  target_link_libraries(wasl-scalar-advection-test PRIVATE wasl_core)
  add_test(NAME wasl.scalar-advection COMMAND wasl-scalar-advection-test)

  add_executable(wasl-euler-test tests/equations/euler_test.cpp)
  target_link_libraries(wasl-euler-test PRIVATE wasl_core)
  add_test(NAME wasl.euler COMMAND wasl-euler-test)

  add_executable(wasl-residual-test tests/dg/residual_test.cpp)
  target_link_libraries(wasl-residual-test PRIVATE wasl_core)
  add_test(NAME wasl.dg.residual COMMAND wasl-residual-test)

  add_executable(wasl-triangle-residual-test tests/dg/triangle_residual_test.cpp)
  target_link_libraries(wasl-triangle-residual-test PRIVATE wasl_core)
  add_test(NAME wasl.dg.triangle-residual COMMAND wasl-triangle-residual-test)

  add_executable(wasl-euler-triangle-residual-test tests/dg/euler_triangle_residual_test.cpp)
  target_link_libraries(wasl-euler-triangle-residual-test PRIVATE wasl_core)
  add_test(NAME wasl.dg.euler-triangle-residual COMMAND wasl-euler-triangle-residual-test)

  add_executable(wasl-euler-line-residual-test tests/dg/euler_line_residual_test.cpp)
  target_link_libraries(wasl-euler-line-residual-test PRIVATE wasl_core)
  add_test(NAME wasl.dg.euler-line-residual COMMAND wasl-euler-line-residual-test)

  add_executable(wasl-euler-case-test tests/runtime/euler_case_test.cpp)
  target_link_libraries(wasl-euler-case-test PRIVATE wasl_core)
  add_test(NAME wasl.runtime.euler-case COMMAND wasl-euler-case-test)

  add_executable(wasl-modal-decay-test tests/limiter/modal_decay_test.cpp)
  target_link_libraries(wasl-modal-decay-test PRIVATE wasl_core)
  add_test(NAME wasl.limiter.modal-decay COMMAND wasl-modal-decay-test)

  add_executable(wasl-simple-weno-test tests/limiter/simple_weno_test.cpp)
  target_link_libraries(wasl-simple-weno-test PRIVATE wasl_core)
  add_test(NAME wasl.limiter.simple-weno COMMAND wasl-simple-weno-test)

  add_executable(wasl-euler-limiting-test tests/limiter/euler_limiting_test.cpp)
  target_link_libraries(wasl-euler-limiting-test PRIVATE wasl_core)
  add_test(NAME wasl.limiter.euler COMMAND wasl-euler-limiting-test)

  add_executable(wasl-case-io-test tests/io/case_io_test.cpp)
  target_link_libraries(wasl-case-io-test PRIVATE wasl_core)
  add_test(NAME wasl.io.case COMMAND wasl-case-io-test)

  # M2 SSP-RK3 time integration and convergence tests.
  add_executable(wasl-ssp-rk3-test tests/time/ssp_rk3_test.cpp)
  target_link_libraries(wasl-ssp-rk3-test PRIVATE wasl_core)
  add_test(NAME wasl.time.ssp-rk3 COMMAND wasl-ssp-rk3-test)

  # M2 boundary-policy tests.
  add_executable(wasl-boundaries-test tests/boundaries/boundaries_test.cpp)
  target_link_libraries(wasl-boundaries-test PRIVATE wasl_core)
  add_test(NAME wasl.boundaries COMMAND wasl-boundaries-test)

  if(WASL_ENABLE_MPI)
    # M8 MPI halo transport: two ranks reproduce the serial residual on a split
    # mesh. Run under mpirun; the in-process oracle (wasl.mpi.halo) already
    # proves the exchange semantics without MPI.
    add_executable(wasl-halo-mpi-test tests/communication/halo_mpi_test.cpp)
    target_link_libraries(wasl-halo-mpi-test PRIVATE wasl_core)
    target_link_libraries(wasl-halo-mpi-test PRIVATE MPI::MPI_CXX)
    add_test(NAME wasl.mpi.transport
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:wasl-halo-mpi-test> ${MPIEXEC_POSTFLAGS})
  endif()
endif()
