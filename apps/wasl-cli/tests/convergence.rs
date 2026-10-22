use std::path::Path;
use std::process::Command;

#[test]
fn smooth_periodic_tutorial_meets_p2_convergence_threshold() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new("python3")
        .arg(root.join("scripts/convergence.py"))
        .arg(env!("CARGO_BIN_EXE_wasl-cli"))
        .arg(root.join("cases/tutorials/advection/smooth-periodic"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("P2 convergence verified"));
}
