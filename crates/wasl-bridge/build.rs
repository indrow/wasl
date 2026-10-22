use std::env;
use std::process::Command;

fn main() {
    let mut build = cxx_build::bridge("src/lib.rs");
    build
        .file("src/bridge.cpp")
        .file("../../cpp/src/case/coefficient_layout.cpp")
        .file("../../cpp/src/case/euler_case.cpp")
        .file("../../cpp/src/case/initial_state.cpp")
        .file("../../cpp/src/case/periodic_boundary.cpp")
        .file("../../cpp/src/case/prepared_mesh.cpp")
        .file("../../cpp/src/case/scalar_advection_case.cpp")
        .file("../../cpp/src/case/simulation.cpp")
        .file("../../cpp/src/case/validation.cpp")
        .file("../../cpp/src/mesh/msh_reader.cpp")
        .file("../../cpp/src/mesh/line_mesh.cpp")
        .file("../../cpp/src/mesh/partition.cpp")
        .file("../../cpp/src/communication/halo_exchange.cpp")
        .file("../../cpp/src/mesh/triangle_mesh.cpp")
        .file("../../cpp/src/io/case_io.cpp")
        .include("include")
        .include("../../cpp/include")
        .std("c++20");

    let mpi = env::var_os("CARGO_FEATURE_MPI").is_some();
    let mut mpi_link_flags = Vec::new();
    if mpi {
        build
            .compiler("mpicxx")
            .define("WASL_ENABLE_MPI", None)
            .define("OMPI_SKIP_MPICXX", None);
        let compile = Command::new("mpicxx")
            .arg("--showme:compile")
            .output()
            .expect("wasl-bridge: mpicxx is required for the mpi feature");
        if !compile.status.success() {
            panic!("wasl-bridge: mpicxx could not report compile flags");
        }
        for flag in String::from_utf8(compile.stdout)
            .unwrap()
            .split_whitespace()
        {
            if let Some(path) = flag.strip_prefix("-I") {
                build.include(path);
            }
        }
        let link = Command::new("mpicxx")
            .arg("--showme:link")
            .output()
            .expect("wasl-bridge: mpicxx is required for the mpi feature");
        if !link.status.success() {
            panic!("wasl-bridge: mpicxx could not report link flags");
        }
        mpi_link_flags = String::from_utf8(link.stdout)
            .unwrap()
            .split_whitespace()
            .map(str::to_owned)
            .collect();
    }

    let flags = Command::new("pkg-config")
        .args(["--cflags", "--libs", "hdf5"])
        .output()
        .expect("wasl-bridge: pkg-config is required for HDF5");
    if !flags.status.success() {
        panic!("wasl-bridge: HDF5 C++ development files are required");
    }
    let mut search_paths = Vec::new();
    let mut libraries = vec!["hdf5_cpp".to_string()];
    for flag in String::from_utf8(flags.stdout).unwrap().split_whitespace() {
        if let Some(path) = flag.strip_prefix("-I") {
            build.include(path);
        } else if let Some(path) = flag.strip_prefix("-L") {
            search_paths.push(path.to_string());
        } else if let Some(library) = flag.strip_prefix("-l") {
            libraries.push(library.to_string());
        }
    }
    build.compile("wasl-bridge-cxx");
    for path in search_paths {
        println!("cargo:rustc-link-search=native={path}");
    }
    for library in libraries {
        println!("cargo:rustc-link-lib={library}");
    }
    for flag in mpi_link_flags {
        if let Some(path) = flag.strip_prefix("-L") {
            println!("cargo:rustc-link-search=native={path}");
        } else if let Some(library) = flag.strip_prefix("-l") {
            println!("cargo:rustc-link-lib={library}");
        }
    }

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/config.rs");
    println!("cargo:rerun-if-changed=src/ops.rs");
    println!("cargo:rerun-if-changed=src/bridge.cpp");
    println!("cargo:rerun-if-changed=include/bridge.hpp");
    println!("cargo:rerun-if-changed=../../cpp/src/io/case_io.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/coefficient_layout.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/euler_case.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/initial_state.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/periodic_boundary.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/prepared_mesh.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/scalar_advection_case.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/simulation.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/case/validation.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/mesh/partition.cpp");
    println!("cargo:rerun-if-changed=../../cpp/src/communication/halo_exchange.cpp");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_MPI");
}
