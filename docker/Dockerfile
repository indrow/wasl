FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang-format \
    cmake \
    git \
    libhdf5-dev \
    libopenmpi-dev \
    ninja-build \
    openmpi-bin \
    pkg-config \
    python3 \
    rustc \
    rustfmt \
    cargo \
 && rm -rf /var/lib/apt/lists/*

# OpenMPI refuses to run as root unless explicitly allowed. The dev container
# runs as root, so allow it for interactive `mpirun` and for the MPI CTest
# preset.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /work
COPY . .
RUN cmake --preset cpu-debug \
 && cmake --build --preset cpu-debug \
 && ctest --preset cpu-debug \
 && cmake --preset cpu-mpi \
 && cmake --build --preset cpu-mpi \
 && cargo test --workspace

CMD ["ctest", "--preset", "cpu-debug"]
