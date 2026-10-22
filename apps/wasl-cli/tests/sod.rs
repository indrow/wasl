use std::path::Path;
use std::process::Command;

#[test]
fn sod_shock_tube_tutorial_converges_to_the_exact_riemann_solution() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new("python3")
        .arg(root.join("scripts/sod_shock_tube.py"))
        .arg(env!("CARGO_BIN_EXE_wasl-cli"))
        .arg(root.join("cases/tutorials/euler/sod-shock-tube"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("Shock tube verified"));
}
