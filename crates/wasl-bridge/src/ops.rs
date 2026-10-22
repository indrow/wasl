use crate::config::{checked_relative, effective_toml, EffectiveConfig};
use crate::ffi;
use std::fs;
use std::path::{Path, PathBuf};

pub fn prepare_mesh(source: &Path, output: &Path) -> Result<String, String> {
    ffi::prepare_mesh(&source.to_string_lossy(), &output.to_string_lossy())
        .map_err(|e| e.to_string())
}

pub fn inspect_checkpoint(path: &Path) -> Result<String, String> {
    ffi::inspect_checkpoint(&path.to_string_lossy()).map_err(|e| e.to_string())
}

pub fn check_prepared_mesh(case: &Path, config: &EffectiveConfig) -> Result<String, String> {
    let path = checked_relative(case, &config.mesh.file, "mesh.file")?;
    ffi::inspect_prepared_mesh(&path.to_string_lossy()).map_err(|e| format!("mesh.file: {e}"))
}

pub fn export_checkpoint(checkpoint: &Path, mesh: &Path, output: &Path) -> Result<String, String> {
    ffi::export_checkpoint(
        &checkpoint.to_string_lossy(),
        &mesh.to_string_lossy(),
        &output.to_string_lossy(),
    )
    .map_err(|e| e.to_string())
}

pub fn run_case(
    case: &Path,
    config: &EffectiveConfig,
    output: &Path,
    restart: Option<&Path>,
) -> Result<String, String> {
    check_run_support(config)?;
    let mesh = checked_relative(case, &config.mesh.file, "mesh.file")?
        .canonicalize()
        .map_err(|e| format!("mesh.file: {e}"))?;
    if config.equation.kind == "euler" {
        let mut regions = Vec::new();
        for region in &config.regions {
            let shape = match region.where_.as_str() {
                "default" => 0.0,
                "interval" => 1.0,
                "circle" => 2.0,
                _ => unreachable!(),
            };
            let velocity = region.velocity.as_ref().unwrap();
            regions.extend_from_slice(&[
                shape,
                region.lower.as_ref().map_or(0.0, |value| value[0]),
                region.upper.as_ref().map_or(0.0, |value| value[0]),
                region.center.as_ref().map_or(0.0, |value| value[0]),
                region.center.as_ref().map_or(0.0, |value| value[1]),
                region.radius.unwrap_or(0.0),
                region.density.unwrap(),
                velocity[0],
                velocity.get(1).copied().unwrap_or(0.0),
                region.pressure.unwrap(),
            ]);
        }
        return ffi::run_euler_case(
            &mesh.to_string_lossy(),
            &output.to_string_lossy(),
            &restart.unwrap_or(Path::new("")).to_string_lossy(),
            config.limiting.enabled,
            config.limiting.sensor_threshold,
            config.limiting.density_floor,
            config.limiting.pressure_floor,
            config.equation.dimensions as i32,
            config.dg.polynomial_order as i32,
            config.equation.heat_capacity_ratio.unwrap(),
            config.time.cfl,
            config.run.start_time,
            config.run.end_time,
            config.run.max_steps,
            config.output.history.interval_steps,
            config.output.visualization.interval_steps,
            config.output.checkpoint.interval_steps,
            config.output.checkpoint.keep_last,
            &regions,
        )
        .map_err(|e| format!("run: {e}"));
    }

    let mut initial_default = None;
    let mut initial_sine = None;
    let mut intervals = Vec::new();
    let mut circles = Vec::new();
    for region in &config.regions {
        match region.where_.as_str() {
            "default" => match initial_default {
                Some(_) => {
                    return Err("initial.regions: multiple default regions are unsupported".into())
                }
                None => initial_default = region.value,
            },
            "interval" => {
                intervals.extend_from_slice(&[
                    region.lower.as_ref().unwrap()[0],
                    region.upper.as_ref().unwrap()[0],
                    region.value.unwrap(),
                ]);
            }
            "circle" => {
                circles.extend_from_slice(&[
                    region.center.as_ref().unwrap()[0],
                    region.center.as_ref().unwrap()[1],
                    region.radius.unwrap(),
                    region.value.unwrap(),
                ]);
            }
            "sine" => match initial_sine {
                Some(_) => {
                    return Err("initial.regions: multiple sine regions are unsupported".into())
                }
                None => initial_sine = region.value,
            },
            _ => unreachable!(),
        }
    }
    ffi::run_scalar_advection_case(
        &mesh.to_string_lossy(),
        &output.to_string_lossy(),
        &restart.unwrap_or(Path::new("")).to_string_lossy(),
        config.parallel.kind == "mpi",
        config.equation.dimensions as i32,
        config.dg.polynomial_order as i32,
        config.equation.velocity[0],
        config.equation.velocity.get(1).copied().unwrap_or(0.0),
        config.time.cfl,
        config.run.start_time,
        config.run.end_time,
        config.run.max_steps,
        config.output.history.interval_steps,
        config.output.visualization.interval_steps,
        config.output.checkpoint.interval_steps,
        config.output.checkpoint.keep_last,
        initial_default.ok_or("initial.regions: one default region is required")?,
        initial_sine.is_some(),
        initial_sine.unwrap_or(0.0),
        &intervals,
        &circles,
    )
    .map_err(|e| format!("run: {e}"))
}

pub fn check_run_support(config: &EffectiveConfig) -> Result<(), String> {
    if !matches!(config.equation.dimensions, 1 | 2)
        || !matches!(config.equation.kind.as_str(), "scalar_advection" | "euler")
        || !matches!(config.parallel.kind.as_str(), "none" | "mpi")
        || config
            .boundaries
            .values()
            .any(|value| value.kind != "periodic")
    {
        return Err("run: supported dispatch requires a 1D or 2D scalar-advection or Euler case with periodic boundaries".into());
    }
    if config.parallel.kind == "mpi" && !cfg!(feature = "mpi") {
        return Err("run: MPI case requires building wasl-cli with --features mpi".into());
    }
    if config.equation.kind == "euler" && config.parallel.kind == "mpi" {
        return Err("run: case-level Euler with MPI is not implemented".into());
    }
    if config.equation.kind == "scalar_advection" && config.limiting.enabled {
        return Err("run: scalar-advection case limiting is not implemented".into());
    }
    for region in &config.regions {
        match (
            config.equation.kind.as_str(),
            config.equation.dimensions,
            region.where_.as_str(),
        ) {
            ("euler", _, "default") | ("euler", 1, "interval") | ("euler", 2, "circle") => {}
            ("euler", 1, _) => {
                return Err("initial.regions: 1D Euler dispatch supports intervals only".into())
            }
            ("euler", 2, _) => {
                return Err("initial.regions: 2D Euler dispatch supports circles only".into())
            }
            ("scalar_advection", _, "default")
            | ("scalar_advection", 1, "interval")
            | ("scalar_advection", 1, "sine")
            | ("scalar_advection", 2, "circle") => {}
            ("scalar_advection", 1, _) => {
                return Err("initial.regions: 1D dispatch supports intervals and sine only".into())
            }
            ("scalar_advection", 2, _) => {
                return Err("initial.regions: 2D dispatch supports circles only".into())
            }
            _ => unreachable!(),
        }
    }
    if config
        .regions
        .iter()
        .filter(|region| region.where_ == "default")
        .count()
        != 1
    {
        return Err("initial.regions: one default region is required".into());
    }
    Ok(())
}

pub fn inspect_prepared_mesh(path: &Path) -> Result<String, String> {
    ffi::inspect_prepared_mesh(&path.to_string_lossy()).map_err(|e| e.to_string())
}

pub fn prepare_mpi_run_directory(root: &Path, effective_config: &str) -> Result<PathBuf, String> {
    ffi::prepare_mpi_run_directory(&root.to_string_lossy(), effective_config)
        .map(PathBuf::from)
        .map_err(|e| e.to_string())
}

pub fn create_run_directory(case: &Path, config: &EffectiveConfig) -> Result<PathBuf, String> {
    let root = checked_relative(case, &config.output.directory, "output.directory")?;
    if config.parallel.kind == "mpi" {
        let effective = effective_toml(config)?;
        return prepare_mpi_run_directory(&root, &effective)
            .map_err(|e| format!("output.directory: {e}"));
    }
    fs::create_dir_all(&root).map_err(|e| format!("output.directory: {e}"))?;
    for number in 1..=999_999 {
        let directory = root.join(format!("run-{number:06}"));
        match fs::create_dir(&directory) {
            Ok(()) => {
                fs::write(directory.join("effective.toml"), effective_toml(config)?)
                    .map_err(|e| format!("effective configuration: {e}"))?;
                return Ok(directory);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(format!("output.directory: {error}")),
        }
    }
    Err("output.directory: no available run directory".into())
}
