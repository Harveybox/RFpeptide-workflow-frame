#!/usr/bin/env python3
"""Dummy AfCycDesign entrypoint for release/demo workflow generation."""
from pathlib import Path


def main():
    Path("dummy_afcyc_results.csv").write_text("candidate,pae_interaction\nDUMMY,0.0\n", encoding="utf-8")
    print("Dummy AfCycDesign script ran. Replace this file path with the real AfCycDesign script for production.")


if __name__ == "__main__":
    main()
