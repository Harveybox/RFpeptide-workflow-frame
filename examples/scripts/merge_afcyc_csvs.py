#!/usr/bin/env python3
"""Dummy AfCycDesign CSV merge entrypoint for release/demo workflow generation."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard_root", default="")
    parser.add_argument("--out_csv", default="results_merged_DUMMY_TARGET_pilot0.csv")
    args = parser.parse_args()
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_csv).write_text("candidate,pae_interaction,rmsd\nDUMMY,0.0,0.0\n", encoding="utf-8")
    print("Dummy AfCycDesign merge script ran. Replace this file path with the real merge script for production.")


if __name__ == "__main__":
    main()
