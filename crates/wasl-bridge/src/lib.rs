//! Coarse Rust/C++ bridge for case loading and case-level operations.
//!
//! Configuration loading lives in `config`. FFI operations live in `ops`.

mod config;
mod ops;

#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("wasl-bridge/include/bridge.hpp");

        fn prepare_mesh(source_path: &str, output_path: &str) -> Result<String>;
        fn inspect_prepared_mesh(mesh_path: &str) -> Result<String>;
        fn inspect_checkpoint(checkpoint_path: &str) -> Result<String>;
        fn export_checkpoint(
            checkpoint_path: &str,
            mesh_path: &str,
            output_path: &str,
        ) -> Result<String>;
        fn prepare_mpi_run_directory(root_path: &str, effective_config: &str) -> Result<String>;
        fn run_scalar_advection_case(
            mesh_path: &str,
            output_path: &str,
            restart_path: &str,
            mpi: bool,
            dimensions: i32,
            polynomial_order: i32,
            velocity_x: f64,
            velocity_y: f64,
            cfl: f64,
            start_time: f64,
            end_time: f64,
            max_steps: u64,
            history_interval: u64,
            visualization_interval: u64,
            checkpoint_interval: u64,
            keep_checkpoints: u64,
            initial_default: f64,
            initial_sine: bool,
            initial_sine_amplitude: f64,
            intervals: &[f64],
            circles: &[f64],
        ) -> Result<String>;
        fn run_euler_case(
            mesh_path: &str,
            output_path: &str,
            restart_path: &str,
            limiting: bool,
            sensor_threshold: f64,
            density_floor: f64,
            pressure_floor: f64,
            dimensions: i32,
            polynomial_order: i32,
            heat_capacity_ratio: f64,
            cfl: f64,
            start_time: f64,
            end_time: f64,
            max_steps: u64,
            history_interval: u64,
            visualization_interval: u64,
            checkpoint_interval: u64,
            keep_checkpoints: u64,
            regions: &[f64],
        ) -> Result<String>;
    }
}

pub use config::{
    effective_toml, load_case, Boundary, CheckpointOutput, Dg, EffectiveConfig, Equation, Flux,
    Interval, Limiting, Mesh, Output, Parallel, Quadrature, Region, Run, Time, Visualization,
};
pub use ops::{
    check_prepared_mesh, check_run_support, create_run_directory, export_checkpoint,
    inspect_checkpoint, inspect_prepared_mesh, prepare_mesh, run_case,
};

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::{Path, PathBuf};

    fn tutorial() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../cases/tutorials/advection/line-pulse")
    }

    #[test]
    fn tutorial_is_checked_and_serialized_deterministically() {
        let config = load_case(&tutorial()).unwrap();
        let first = effective_toml(&config).unwrap();
        assert_eq!(first, effective_toml(&config).unwrap());
        assert!(first.contains("surface_order = 5"));
    }

    #[test]
    fn tutorial_is_supported_by_dispatch() {
        let config = load_case(&tutorial()).unwrap();
        assert!(check_run_support(&config).is_ok());
    }

    #[test]
    fn malformed_and_unsupported_fields_are_named() {
        let mut file: toml::Value = fs::read_to_string(tutorial().join("numerics.toml"))
            .unwrap()
            .parse()
            .unwrap();
        file["dg"]["polynomial_order"] = toml::Value::Integer(9);
        let directory = std::env::temp_dir().join("wasl-config-error-test");
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(directory.join("mesh")).unwrap();
        fs::create_dir_all(directory.join("initial")).unwrap();
        for name in ["case.toml", "physics.toml", "solver.toml"] {
            fs::copy(tutorial().join(name), directory.join(name)).unwrap();
        }
        fs::copy(
            tutorial().join("mesh/boundaries.toml"),
            directory.join("mesh/boundaries.toml"),
        )
        .unwrap();
        fs::copy(
            tutorial().join("initial/state.toml"),
            directory.join("initial/state.toml"),
        )
        .unwrap();
        fs::write(
            directory.join("numerics.toml"),
            toml::to_string(&file).unwrap(),
        )
        .unwrap();
        assert!(load_case(&directory)
            .unwrap_err()
            .contains("dg.polynomial_order"));

        fs::write(directory.join("numerics.toml"), "format_version = [").unwrap();
        assert!(load_case(&directory)
            .unwrap_err()
            .contains("includes.numerics"));

        fs::copy(
            tutorial().join("numerics.toml"),
            directory.join("numerics.toml"),
        )
        .unwrap();
        let mut physics: toml::Value = fs::read_to_string(tutorial().join("physics.toml"))
            .unwrap()
            .parse()
            .unwrap();
        physics["equation"]
            .as_table_mut()
            .unwrap()
            .remove("velocity");
        fs::write(
            directory.join("physics.toml"),
            toml::to_string(&physics).unwrap(),
        )
        .unwrap();
        assert!(load_case(&directory).unwrap_err().contains("velocity"));
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn coarse_bridge_prepares_and_reads_a_line_mesh() {
        const MSH: &str = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
0 2 "left"
0 3 "right"
1 1 "fluid"
$EndPhysicalNames
$Entities
2 1 0 0
1 0 0 0 1 2
2 1 0 0 1 3
1 0 0 0 1 0 0 1 1 2 1 -2
$EndEntities
$Nodes
2 2 1 2
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
$EndNodes
$Elements
1 1 1 1
1 1 1 1
1 1 2
$EndElements
"#;
        let directory = std::env::temp_dir().join("wasl-bridge-prepare-test");
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        let source = directory.join("source.msh");
        let prepared = directory.join("prepared.h5");
        fs::write(&source, MSH).unwrap();
        assert_eq!(
            prepare_mesh(&source, &prepared).unwrap(),
            "prepared 1 elements"
        );
        assert_eq!(
            inspect_prepared_mesh(&prepared).unwrap(),
            "1D mesh with 1 elements"
        );
        fs::remove_dir_all(directory).unwrap();
    }
}
