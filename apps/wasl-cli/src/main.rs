use std::env;
use std::path::Path;

fn usage() -> String {
    "usage:\n  wasl check CASE\n  wasl run CASE\n  wasl restart CASE CHECKPOINT\n  wasl mesh prepare SOURCE.msh PREPARED.h5\n  wasl inspect CHECKPOINT.h5\n  wasl export CHECKPOINT.h5 PREPARED.h5 OUTPUT.vtu".into()
}

fn command(arguments: &[String]) -> Result<String, String> {
    match arguments {
        [command, case] if command == "check" => {
            let config = wasl_bridge::load_case(Path::new(case))?;
            wasl_bridge::effective_toml(&config)
        }
        [command, case] if command == "run" => {
            let config = wasl_bridge::load_case(Path::new(case))?;
            wasl_bridge::check_prepared_mesh(Path::new(case), &config)?;
            wasl_bridge::check_run_support(&config)?;
            let directory = wasl_bridge::create_run_directory(Path::new(case), &config)?;
            let summary = wasl_bridge::run_case(Path::new(case), &config, &directory, None)?;
            if summary.is_empty() {
                return Ok(String::new());
            }
            Ok(format!("created {} {summary}", directory.display()))
        }
        [command, case, checkpoint] if command == "restart" => {
            let summary = wasl_bridge::inspect_checkpoint(Path::new(checkpoint))?;
            let config = wasl_bridge::load_case(Path::new(case))?;
            wasl_bridge::check_prepared_mesh(Path::new(case), &config)?;
            wasl_bridge::check_run_support(&config)?;
            let directory = wasl_bridge::create_run_directory(Path::new(case), &config)?;
            if config.parallel.kind == "none" {
                std::fs::write(
                    directory.join("restart.txt"),
                    format!("{checkpoint}\n{summary}\n"),
                )
                .map_err(|e| format!("restart: {e}"))?;
            }
            let run = wasl_bridge::run_case(
                Path::new(case),
                &config,
                &directory,
                Some(Path::new(checkpoint)),
            )?;
            if run.is_empty() {
                return Ok(String::new());
            }
            Ok(format!(
                "created {} from {summary} {run}",
                directory.display()
            ))
        }
        [mesh, prepare, source, output] if mesh == "mesh" && prepare == "prepare" => {
            wasl_bridge::prepare_mesh(Path::new(source), Path::new(output))
        }
        [command, checkpoint] if command == "inspect" => {
            wasl_bridge::inspect_checkpoint(Path::new(checkpoint))
        }
        [command, checkpoint, mesh, output] if command == "export" => {
            wasl_bridge::export_checkpoint(
                Path::new(checkpoint),
                Path::new(mesh),
                Path::new(output),
            )
        }
        _ => Err(usage()),
    }
}

fn main() {
    let arguments: Vec<_> = env::args().skip(1).collect();
    match command(&arguments) {
        Ok(message) if !message.is_empty() => println!("{message}"),
        Ok(_) => {}
        Err(error) => {
            eprintln!("wasl: {error}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

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
3 2 0 0
1 0 0 0 1 2
2 0.5 0 0 0
3 1 0 0 1 3
1 0 0 0 0.5 0 0 1 1 2 1 -2
2 0.5 0 0 1 0 0 1 1 2 2 -3
$EndEntities
$Nodes
3 3 1 3
0 1 0 1
1
0 0 0
0 2 0 1
2
0.5 0 0
0 3 0 1
3
1 0 0
$EndNodes
$Elements
4 4 1 4
0 1 15 1
1 1
0 3 15 1
2 3
1 1 1 1
3 1 2
1 2 1 1
4 2 3
$EndElements
"#;

    const TRIANGLE_MSH: &str = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
1 2 "domain"
2 1 "fluid"
$EndPhysicalNames
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 2 2 1 -2
2 1 0 0 1 1 0 1 2 2 2 -3
3 0 1 0 1 1 0 1 2 2 4 -3
4 0 0 0 0 1 0 1 2 2 1 -4
1 0 0 0 1 1 0 1 1 4 1 2 3 4
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
1 1 2 3
2 1 3 4
$EndElements
"#;

    fn copy_tutorial(destination: &Path) {
        let source = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../cases/tutorials/advection/line-pulse");
        fs::create_dir_all(destination.join("mesh")).unwrap();
        fs::create_dir_all(destination.join("initial")).unwrap();
        for name in ["case.toml", "physics.toml", "numerics.toml"] {
            fs::copy(source.join(name), destination.join(name)).unwrap();
        }
        fs::copy(
            source.join("mesh/boundaries.toml"),
            destination.join("mesh/boundaries.toml"),
        )
        .unwrap();
        fs::copy(
            source.join("initial/state.toml"),
            destination.join("initial/state.toml"),
        )
        .unwrap();
        let solver = fs::read_to_string(source.join("solver.toml"))
            .unwrap()
            .replace("end_time = 1.0", "end_time = 0.08")
            .replace("interval_steps = 1000", "interval_steps = 1")
            .replace("interval_steps = 100", "interval_steps = 1");
        fs::write(destination.join("solver.toml"), solver).unwrap();
        fs::write(destination.join("mesh/source.msh"), MSH).unwrap();
        wasl_bridge::prepare_mesh(
            &destination.join("mesh/source.msh"),
            &destination.join("mesh/prepared.h5"),
        )
        .unwrap();
    }

    fn copy_rotating_cone(destination: &Path) {
        let source = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../cases/tutorials/advection/rotating-cone");
        fs::create_dir_all(destination.join("mesh")).unwrap();
        fs::create_dir_all(destination.join("initial")).unwrap();
        for name in ["case.toml", "physics.toml", "numerics.toml"] {
            fs::copy(source.join(name), destination.join(name)).unwrap();
        }
        fs::copy(
            source.join("mesh/boundaries.toml"),
            destination.join("mesh/boundaries.toml"),
        )
        .unwrap();
        fs::copy(
            source.join("initial/state.toml"),
            destination.join("initial/state.toml"),
        )
        .unwrap();
        let solver = fs::read_to_string(source.join("solver.toml"))
            .unwrap()
            .replace("end_time = 1.0", "end_time = 0.02")
            .replace("interval_steps = 1000", "interval_steps = 1")
            .replace("interval_steps = 100", "interval_steps = 1");
        fs::write(destination.join("solver.toml"), solver).unwrap();
        fs::write(destination.join("mesh/source.msh"), TRIANGLE_MSH).unwrap();
        wasl_bridge::prepare_mesh(
            &destination.join("mesh/source.msh"),
            &destination.join("mesh/prepared.h5"),
        )
        .unwrap();
    }

    fn copy_euler(destination: &Path) {
        copy_tutorial(destination);
        fs::write(
            destination.join("physics.toml"),
            "format_version = 1\n\n[equation]\ntype = \"euler\"\ndimensions = 1\nheat_capacity_ratio = 1.4\n",
        )
        .unwrap();
        let numerics = fs::read_to_string(destination.join("numerics.toml"))
            .unwrap()
            .replace("type = \"upwind\"", "type = \"hllc\"");
        fs::write(destination.join("numerics.toml"), numerics).unwrap();
        fs::write(
            destination.join("initial/state.toml"),
            "format_version = 1\nrepresentation = \"primitive\"\n\n[[regions]]\nname = \"free-stream\"\nwhere = \"default\"\ndensity = 1.1\nvelocity = [0.4]\npressure = 0.9\n",
        )
        .unwrap();
    }

    fn only_checkpoint(directory: &Path) -> PathBuf {
        fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| {
                path.file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("checkpoint-")
            })
            .max()
            .unwrap()
    }

    fn masses(directory: &Path) -> Vec<f64> {
        fs::read_to_string(directory.join("history.csv"))
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| line.split(',').nth(3).unwrap().parse().unwrap())
            .collect()
    }

    #[test]
    fn usage_is_returned_for_unknown_command() {
        assert!(command(&["unknown".into()])
            .unwrap_err()
            .starts_with("usage:"));
    }

    #[test]
    fn run_rejects_missing_prepared_mesh_before_creating_results() {
        let case = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../cases/tutorials/advection/line-pulse");
        let error = command(&["run".into(), case.to_string_lossy().into_owned()]).unwrap_err();
        assert!(error.contains("mesh.file"));
    }

    #[test]
    fn run_and_restart_advance_the_tutorial_case() {
        let case = std::env::temp_dir().join(format!("wasl-cli-run-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&case);
        copy_tutorial(&case);

        let run = command(&["run".into(), case.to_string_lossy().into_owned()]).unwrap();
        assert!(run.contains("time=0.080000000000000002"));
        let first = case.join("results/run-000001");
        assert!(first.join("history.csv").is_file());
        assert!(first.join("solution.pvd").is_file());
        let checkpoint = only_checkpoint(&first);

        let solver = fs::read_to_string(case.join("solver.toml"))
            .unwrap()
            .replace("end_time = 0.08", "end_time = 0.12");
        fs::write(case.join("solver.toml"), solver).unwrap();
        let restart = command(&[
            "restart".into(),
            case.to_string_lossy().into_owned(),
            checkpoint.to_string_lossy().into_owned(),
        ])
        .unwrap();
        assert!(restart.contains("time=0.12"));
        let second = case.join("results/run-000002");
        assert!(second.join("restart.txt").is_file());
        assert!(only_checkpoint(&second).is_file());
        let fresh_mass = masses(&first)[0];
        for value in masses(&first).into_iter().chain(masses(&second)) {
            assert!((value - fresh_mass).abs() < 1e-12);
        }
        fs::remove_dir_all(case).unwrap();
    }

    #[test]
    fn rotating_cone_tutorial_runs_and_conserves_mass() {
        let case = std::env::temp_dir().join(format!("wasl-cli-cone-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&case);
        copy_rotating_cone(&case);

        let run = command(&["run".into(), case.to_string_lossy().into_owned()]).unwrap();
        assert!(run.contains("time=0.02"));
        let output = case.join("results/run-000001");
        assert!(output.join("solution.pvd").is_file());
        assert!(only_checkpoint(&output).is_file());
        let mass = masses(&output);
        assert!(mass.iter().all(|value| (value - mass[0]).abs() < 1e-12));
        fs::remove_dir_all(case).unwrap();
    }

    #[test]
    fn euler_case_runs_through_the_cli() {
        let case = std::env::temp_dir().join(format!("wasl-cli-euler-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&case);
        copy_euler(&case);
        let run = command(&["run".into(), case.to_string_lossy().into_owned()]).unwrap();
        assert!(run.contains("time=0.080000000000000002"));
        let checkpoint = only_checkpoint(&case.join("results/run-000001"));
        assert!(wasl_bridge::inspect_checkpoint(&checkpoint)
            .unwrap()
            .contains("variables=4"));
        fs::remove_dir_all(case).unwrap();
    }
}
