use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Component, Path, PathBuf};

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct CaseFile {
    format_version: u32,
    name: String,
    mesh: Mesh,
    includes: Includes,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Mesh {
    pub file: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Includes {
    physics: String,
    numerics: String,
    solver: String,
    boundaries: String,
    initial: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct PhysicsFile {
    format_version: u32,
    equation: Equation,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Equation {
    #[serde(rename = "type")]
    pub kind: String,
    pub dimensions: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub velocity: Vec<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heat_capacity_ratio: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct NumericsFile {
    format_version: u32,
    dg: Dg,
    flux: Flux,
    limiting: Limiting,
    time: Time,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Dg {
    pub formulation: String,
    pub basis: String,
    pub polynomial_order: u32,
    pub geometry_order: u32,
    pub quadrature: Quadrature,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Quadrature {
    pub volume_order: u32,
    #[serde(default)]
    pub surface_order: Option<u32>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Flux {
    #[serde(rename = "type")]
    pub kind: String,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Limiting {
    pub enabled: bool,
    #[serde(default = "default_sensor_threshold")]
    pub sensor_threshold: f64,
    #[serde(default = "default_positivity_floor")]
    pub density_floor: f64,
    #[serde(default = "default_positivity_floor")]
    pub pressure_floor: f64,
}

fn default_sensor_threshold() -> f64 {
    1e-3
}

fn default_positivity_floor() -> f64 {
    1e-10
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Time {
    pub integrator: String,
    pub cfl: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SolverFile {
    format_version: u32,
    run: Run,
    parallel: Parallel,
    output: Output,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Run {
    pub start_time: f64,
    pub end_time: f64,
    pub max_steps: u64,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Parallel {
    #[serde(rename = "type")]
    pub kind: String,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub directory: String,
    pub history: Interval,
    pub visualization: Visualization,
    pub checkpoint: CheckpointOutput,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Interval {
    pub interval_steps: u64,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Visualization {
    pub format: String,
    pub interval_steps: u64,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct CheckpointOutput {
    pub format: String,
    pub interval_steps: u64,
    pub keep_last: u64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryFile {
    format_version: u32,
    boundaries: BTreeMap<String, Boundary>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Boundary {
    #[serde(rename = "type")]
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pair: Option<String>,
    /// Prescribed exterior state of a scalar inflow boundary.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct InitialFile {
    format_version: u32,
    representation: String,
    regions: Vec<Region>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub name: String,
    #[serde(rename = "where")]
    pub where_: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lower: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub upper: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub center: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub velocity: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct EffectiveConfig {
    pub format_version: u32,
    pub name: String,
    pub mesh: Mesh,
    pub equation: Equation,
    pub dg: Dg,
    pub flux: Flux,
    pub limiting: Limiting,
    pub time: Time,
    pub run: Run,
    pub parallel: Parallel,
    pub output: Output,
    pub boundaries: BTreeMap<String, Boundary>,
    pub initial_representation: String,
    pub regions: Vec<Region>,
}

pub(crate) fn parse<T: for<'de> Deserialize<'de>>(
    case: &Path,
    relative: &str,
    field: &str,
) -> Result<T, String> {
    let path = checked_relative(case, relative, field)?;
    let text = fs::read_to_string(&path)
        .map_err(|e| format!("{field}: cannot read {}: {e}", path.display()))?;
    toml::from_str(&text).map_err(|e| format!("{field}: {e}"))
}

pub(crate) fn checked_relative(case: &Path, value: &str, field: &str) -> Result<PathBuf, String> {
    let path = Path::new(value);
    if path.is_absolute()
        || path
            .components()
            .any(|part| matches!(part, Component::ParentDir))
    {
        return Err(format!("{field}: path must stay inside the case directory"));
    }
    Ok(case.join(path))
}

pub(crate) fn version(value: u32, field: &str) -> Result<(), String> {
    if value != 1 {
        return Err(format!("{field}.format_version: unsupported value {value}"));
    }
    Ok(())
}

pub(crate) fn one_of(value: &str, supported: &[&str], field: &str) -> Result<(), String> {
    if !supported.contains(&value) {
        return Err(format!("{field}: unsupported value {value:?}"));
    }
    Ok(())
}

pub fn load_case(case: &Path) -> Result<EffectiveConfig, String> {
    let root: CaseFile = parse(case, "case.toml", "case")?;
    version(root.format_version, "case")?;
    if root.name.trim().is_empty() {
        return Err("case.name: must not be empty".into());
    }
    if root.mesh.file.is_empty() {
        return Err("mesh.file: must not be empty".into());
    }
    checked_relative(case, &root.mesh.file, "mesh.file")?;
    let physics: PhysicsFile = parse(case, &root.includes.physics, "includes.physics")?;
    let mut numerics: NumericsFile = parse(case, &root.includes.numerics, "includes.numerics")?;
    let solver: SolverFile = parse(case, &root.includes.solver, "includes.solver")?;
    let boundaries: BoundaryFile = parse(case, &root.includes.boundaries, "includes.boundaries")?;
    let initial: InitialFile = parse(case, &root.includes.initial, "includes.initial")?;
    for (value, field) in [
        (physics.format_version, "physics"),
        (numerics.format_version, "numerics"),
        (solver.format_version, "solver"),
        (boundaries.format_version, "boundaries"),
        (initial.format_version, "initial"),
    ] {
        version(value, field)?;
    }

    one_of(
        &physics.equation.kind,
        &["scalar_advection", "euler"],
        "equation.type",
    )?;
    if !(1..=2).contains(&physics.equation.dimensions) {
        return Err("equation.dimensions: supported values are 1 and 2".into());
    }
    if physics.equation.kind == "scalar_advection" {
        if physics.equation.velocity.len() != physics.equation.dimensions
            || physics.equation.velocity.iter().any(|v| !v.is_finite())
        {
            return Err("equation.velocity: must contain one finite value per dimension".into());
        }
        if physics.equation.heat_capacity_ratio.is_some() {
            return Err("equation.heat_capacity_ratio: valid for Euler only".into());
        }
    } else if !physics.equation.velocity.is_empty()
        || physics
            .equation
            .heat_capacity_ratio
            .map_or(true, |value| !value.is_finite() || value <= 1.0)
    {
        return Err(
            "equation.heat_capacity_ratio: Euler requires a finite value greater than 1".into(),
        );
    }
    one_of(&numerics.dg.formulation, &["weak"], "dg.formulation")?;
    let basis = if physics.equation.dimensions == 1 {
        "modal_legendre"
    } else {
        "modal_dubiner"
    };
    one_of(&numerics.dg.basis, &[basis], "dg.basis")?;
    if numerics.dg.polynomial_order > 3 {
        return Err("dg.polynomial_order: supported range is 0 through 3".into());
    }
    if numerics.dg.geometry_order != 1 {
        return Err("dg.geometry_order: only straight-sided order 1 is supported".into());
    }
    if numerics.dg.quadrature.volume_order == 0 {
        return Err("dg.quadrature.volume_order: must be positive".into());
    }
    let surface = numerics
        .dg
        .quadrature
        .surface_order
        .get_or_insert(numerics.dg.quadrature.volume_order);
    if *surface == 0 {
        return Err("dg.quadrature.surface_order: must be positive".into());
    }
    let flux = if physics.equation.kind == "scalar_advection" {
        &["upwind"][..]
    } else {
        &["hllc"][..]
    };
    one_of(&numerics.flux.kind, flux, "flux.type")?;
    one_of(&numerics.time.integrator, &["ssp_rk3"], "time.integrator")?;
    if !numerics.limiting.sensor_threshold.is_finite()
        || numerics.limiting.sensor_threshold < 0.0
        || !numerics.limiting.density_floor.is_finite()
        || numerics.limiting.density_floor < 0.0
        || !numerics.limiting.pressure_floor.is_finite()
        || numerics.limiting.pressure_floor < 0.0
    {
        return Err(
            "limiting: threshold and positivity floors must be finite and nonnegative".into(),
        );
    }
    if !numerics.time.cfl.is_finite() || numerics.time.cfl <= 0.0 {
        return Err("time.cfl: must be finite and positive".into());
    }
    if !solver.run.start_time.is_finite()
        || !solver.run.end_time.is_finite()
        || solver.run.end_time <= solver.run.start_time
    {
        return Err("run.end_time: must be finite and greater than run.start_time".into());
    }
    if solver.run.max_steps == 0 {
        return Err("run.max_steps: must be positive".into());
    }
    one_of(&solver.parallel.kind, &["none", "mpi"], "parallel.type")?;
    if solver.output.directory.is_empty() {
        return Err("output.directory: must not be empty".into());
    }
    checked_relative(case, &solver.output.directory, "output.directory")?;
    for (value, field) in [
        (
            solver.output.history.interval_steps,
            "output.history.interval_steps",
        ),
        (
            solver.output.visualization.interval_steps,
            "output.visualization.interval_steps",
        ),
        (
            solver.output.checkpoint.interval_steps,
            "output.checkpoint.interval_steps",
        ),
        (
            solver.output.checkpoint.keep_last,
            "output.checkpoint.keep_last",
        ),
    ] {
        if value == 0 {
            return Err(format!("{field}: must be positive"));
        }
    }
    one_of(
        &solver.output.visualization.format,
        &["vtk_xml"],
        "output.visualization.format",
    )?;
    one_of(
        &solver.output.checkpoint.format,
        &["hdf5"],
        "output.checkpoint.format",
    )?;
    if boundaries.boundaries.is_empty() {
        return Err("boundaries: must contain at least one boundary".into());
    }
    for (name, boundary) in &boundaries.boundaries {
        one_of(
            &boundary.kind,
            &["periodic", "inflow", "outflow", "slip_wall"],
            &format!("boundaries.{name}.type"),
        )?;
        if boundary.kind == "periodic" {
            if let Some(pair) = &boundary.pair {
                if !boundaries.boundaries.contains_key(pair) {
                    return Err(format!("boundaries.{name}.pair: unknown boundary {pair:?}"));
                }
            }
        }
        let scalar_inflow =
            boundary.kind == "inflow" && physics.equation.kind == "scalar_advection";
        match boundary.value {
            Some(_) if !scalar_inflow => {
                return Err(format!(
                    "boundaries.{name}.value: valid for scalar-advection inflow only"
                ));
            }
            Some(value) if !value.is_finite() => {
                return Err(format!("boundaries.{name}.value: must be finite"));
            }
            None if scalar_inflow => {
                return Err(format!(
                    "boundaries.{name}.value: scalar-advection inflow requires a value"
                ));
            }
            _ => {}
        }
    }
    let representation = if physics.equation.kind == "scalar_advection" {
        "scalar"
    } else {
        "primitive"
    };
    one_of(
        &initial.representation,
        &[representation],
        "initial.representation",
    )?;
    if initial.regions.is_empty() {
        return Err("initial.regions: must not be empty".into());
    }
    for (index, region) in initial.regions.iter().enumerate() {
        one_of(
            &region.where_,
            &["default", "interval", "circle", "sine"],
            &format!("initial.regions[{index}].where"),
        )?;
        if physics.equation.kind == "scalar_advection" {
            if region.value.map_or(true, |value| !value.is_finite()) {
                return Err(format!("initial.regions[{index}].value: must be finite"));
            }
            if region.density.is_some() || region.velocity.is_some() || region.pressure.is_some() {
                return Err(format!(
                    "initial.regions[{index}]: scalar state accepts value only"
                ));
            }
        } else {
            let velocity = region.velocity.as_deref();
            if region.value.is_some()
                || region
                    .density
                    .map_or(true, |value| !value.is_finite() || value <= 0.0)
                || velocity.map_or(true, |values| {
                    values.len() != physics.equation.dimensions
                        || values.iter().any(|value| !value.is_finite())
                })
                || region
                    .pressure
                    .map_or(true, |value| !value.is_finite() || value <= 0.0)
            {
                return Err(format!("initial.regions[{index}]: Euler primitive state requires positive density and pressure plus one finite velocity per dimension"));
            }
        }
        match region.where_.as_str() {
            "interval" => {
                let Some([lower]) = region.lower.as_deref() else {
                    return Err(format!(
                        "initial.regions[{index}]: interval requires one lower and upper value"
                    ));
                };
                let Some([upper]) = region.upper.as_deref() else {
                    return Err(format!(
                        "initial.regions[{index}]: interval requires one lower and upper value"
                    ));
                };
                if physics.equation.dimensions != 1
                    || !lower.is_finite()
                    || !upper.is_finite()
                    || lower >= upper
                {
                    return Err(format!(
                        "initial.regions[{index}]: interval requires 1D finite bounds with lower < upper"
                    ));
                }
            }
            "circle" => {
                let Some([x, y]) = region.center.as_deref() else {
                    return Err(format!(
                        "initial.regions[{index}]: circle requires a 2D center and positive radius"
                    ));
                };
                if physics.equation.dimensions != 2
                    || !x.is_finite()
                    || !y.is_finite()
                    || region.radius.map_or(true, |r| !r.is_finite() || r <= 0.0)
                {
                    return Err(format!(
                        "initial.regions[{index}]: circle requires a finite 2D center and positive radius"
                    ));
                }
            }
            "sine" if physics.equation.dimensions != 1 => {
                return Err(format!(
                    "initial.regions[{index}]: sine requires a 1D equation"
                ));
            }
            _ => {}
        }
    }

    Ok(EffectiveConfig {
        format_version: 1,
        name: root.name,
        mesh: root.mesh,
        equation: physics.equation,
        dg: numerics.dg,
        flux: numerics.flux,
        limiting: numerics.limiting,
        time: numerics.time,
        run: solver.run,
        parallel: solver.parallel,
        output: solver.output,
        boundaries: boundaries.boundaries,
        initial_representation: initial.representation,
        regions: initial.regions,
    })
}

pub fn effective_toml(config: &EffectiveConfig) -> Result<String, String> {
    toml::to_string_pretty(config).map_err(|e| format!("effective configuration: {e}"))
}
