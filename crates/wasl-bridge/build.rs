use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// The bridge compiles the whole C++ core except the standalone smoke program.
const CORE_SOURCE_DIR: &str = "../../cpp/src";
const EXCLUDED_SOURCES: &[&str] = &["smoke.cpp"];

fn core_sources(directory: &Path) -> Vec<PathBuf> {
    let mut sources = Vec::new();
    let mut entries: Vec<PathBuf> = fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("wasl-bridge: cannot read {}: {error}", directory.display()))
        .map(|entry| entry.expect("wasl-bridge: unreadable source entry").path())
        .collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            sources.extend(core_sources(&path));
        } else if path.extension().is_some_and(|extension| extension == "cpp")
            && !EXCLUDED_SOURCES
                .iter()
                .any(|name| path.file_name().is_some_and(|file| file == *name))
        {
            sources.push(path);
        }
    }
    sources
}

fn main() {
    let mut build = cxx_build::bridge("src/lib.rs");
    build
        .file("src/bridge.cpp")
        .files(core_sources(Path::new(CORE_SOURCE_DIR)))
        .include("include")
        .include("../../cpp/include")
        .include(CORE_SOURCE_DIR)
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
    // A directory path makes Cargo rescan every file below it.
    println!("cargo:rerun-if-changed={CORE_SOURCE_DIR}");
    println!("cargo:rerun-if-changed=../../cpp/include");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_MPI");
}
