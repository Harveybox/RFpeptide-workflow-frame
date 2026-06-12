#!/usr/bin/env python3
"""Dummy PyRosetta scoring entrypoint for release/demo workflow generation."""
from pathlib import Path


def main():
    Path("pyrosetta_scores.csv").write_text("description,dG\nDUMMY,0.0\n", encoding="utf-8")
    print("Dummy PyRosetta script ran. Replace this file path with the real PyRosetta script for production.")


if __name__ == "__main__":
    main()
