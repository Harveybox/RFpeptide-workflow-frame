#!/usr/bin/env python3
"""Dummy AfCycDesign CSV merge entrypoint for release/demo workflow generation."""
from pathlib import Path


def main():
    Path("results_merged.csv").write_text("candidate,pae_interaction,rmsd\nDUMMY,0.0,0.0\n", encoding="utf-8")
    print("Dummy AfCycDesign merge script ran. Replace this file path with the real merge script for production.")


if __name__ == "__main__":
    main()
