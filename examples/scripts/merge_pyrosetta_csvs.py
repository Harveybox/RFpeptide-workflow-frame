#!/usr/bin/env python3
"""Dummy PyRosetta CSV merge entrypoint for release/demo workflow generation."""
from pathlib import Path


def main():
    Path("pyrosetta_scores_merged.csv").write_text("description,dG\nDUMMY,0.0\n", encoding="utf-8")
    print("Dummy PyRosetta merge script ran. Replace this file path with the real merge script for production.")


if __name__ == "__main__":
    main()
