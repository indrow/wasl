use std::path::Path;
use std::process::Command;

#[test]
fn line_inflow_tutorial_balances_mass_and_transports_the_inflow_state() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new("python3")
        .arg(root.join("scripts/inflow_advection.py"))
        .arg(env!("CARGO_BIN_EXE_wasl-cli"))
        .arg(root.join("cases/tutorials/advection/line-inflow"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("Inflow advection verified"));
}
