# Roadmap

## GUI Refinement Backlog

- Add an embedded PDB/PyMOL-style viewer so downloaded structures can be inspected without launching an external PyMOL window.
- Rework the Cluster Dashboard into an automatic workflow orchestrator: on refresh, detect the current stage of a target/pilot/main run from job status, logs, and expected outputs; provide one-click full workflow submission; after each stage finishes, automatically validate job/output health and submit the next stage.

## Completed

- Added independent `Data Processing` filters with separate target/pilot/stage/shard/search state.
- Moved project-level settings that directly control RFDiffusion into the `RFDiffusion` page.
- Added score-item selection controls for `AfCycDesign` and `PyRosetta`, including metric direction and default metric presets.
- Added hit screening from AfCycDesign and PyRosetta merged CSVs, with RFpeptides-inspired thresholds and structure retrieval.
