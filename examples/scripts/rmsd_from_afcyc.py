#!/usr/bin/env python3
"""Dummy RMSD post-processing entrypoint for release/demo workflow generation."""
from pathlib import Path


def main():
    Path("dummy_rmsd_results.csv").write_text("candidate,rmsd\nDUMMY,0.0\n", encoding="utf-8")
    print("Dummy RMSD script ran. Replace this file path with the real RMSD script for production.")


if __name__ == "__main__":
    main()
