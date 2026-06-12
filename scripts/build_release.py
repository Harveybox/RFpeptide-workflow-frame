#!/usr/bin/env python3
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERSION = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
PACKAGE_NAME = f"RFpeptide-workflow-frame-v{VERSION}"
DIST_DIR = ROOT / "dist"
ZIP_PATH = DIST_DIR / f"{PACKAGE_NAME}.zip"

INCLUDE_FILES = [
    ".gitattributes",
    ".gitignore",
    "CHANGELOG.md",
    "README.md",
    "VERSION",
    "cluster_ops.py",
    "cluster_profile.example.json",
    "credential_store.py",
    "environment.yml",
    "generate_workflow.py",
    "launch_gui.ps1",
    "pdb_preprocess.py",
    "requirements.txt",
    "workflow_gui.py",
    "configs/DUMMY_TARGET.json",
    "examples/dummy_target_clean.pdb",
    "examples/scripts/PyRosetta_fullScoring_v4_debug.py",
    "examples/scripts/afcyc_predict_batch.py",
    "examples/scripts/merge_afcyc_csvs.py",
    "examples/scripts/merge_pyrosetta_csvs.py",
    "examples/scripts/rmsd_from_afcyc.py",
    "scripts/build_release.py",
]


def main() -> None:
    DIST_DIR.mkdir(parents=True, exist_ok=True)
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for relative_name in INCLUDE_FILES:
            source = ROOT / relative_name
            if not source.is_file():
                raise FileNotFoundError(source)
            archive.write(source, Path(PACKAGE_NAME) / relative_name)

    size_mb = ZIP_PATH.stat().st_size / (1024 * 1024)
    print(f"Built {ZIP_PATH} ({size_mb:.2f} MB)")

    latest_path = DIST_DIR / "RFpeptide-workflow-frame-latest.zip"
    if latest_path.exists():
        latest_path.unlink()
    shutil.copy2(ZIP_PATH, latest_path)
    print(f"Updated {latest_path}")


if __name__ == "__main__":
    main()
