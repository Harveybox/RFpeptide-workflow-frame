#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
from string import Template


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPT_DIR / "configs" / "PGLYRP1.json"


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sh_quote(value) -> str:
    text = str(value)
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render(template: str, values: dict) -> str:
    return Template(template).safe_substitute(values).rstrip() + "\n"


def write_file(path: Path, content: str, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists; pass --force to overwrite")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")
    if path.suffix == ".sh":
        try:
            os.chmod(path, 0o755)
        except OSError:
            pass


def common_values(config: dict) -> dict:
    target = config["target"]
    pilot = config.get("pilot", "pilot0")
    scratch_root = f"/scratch/{config['scratch_date']}/{config['cluster_user']}/{config['project_dir_name']}"
    home_project_dir = config["home_project_dir"].rstrip("/")
    return {
        "target": target,
        "target_lower": target.lower(),
        "pilot": pilot,
        "scratch_root": scratch_root,
        "home_project_dir": home_project_dir,
        "input_pdb": config["input_pdb"],
        "target_chain": config.get("target_chain", "A"),
        "binder_chain": config.get("binder_chain", "B"),
        "contigs": config["contigs"],
        "micromamba": "$HOME/bin/micromamba",
    }


RFDIFFUSION_TEMPLATE = r"""#!/bin/bash
set -euo pipefail

QUEUE=${rf_queue}
MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME=${rf_env}
RFDIFFUSION_DIR=${rf_dir}

WORKDIR=${scratch_root}
INPUT_PDB=${input_pdb}
OUTPUT_PREFIX="$WORKDIR/example_outputs/${pilot}/diffused_binder_cyclic_${target}_${pilot}"

N_SHARDS=${rf_n_shards}
DESIGNS_PER_SHARD=${rf_designs_per_shard}

mkdir -p "$WORKDIR"
mkdir -p "$(dirname "$OUTPUT_PREFIX")"
cd "$WORKDIR"

for i in $(seq 0 $((N_SHARDS-1))); do
    shard=$(printf "%02d" "$i")
    shard_prefix="${OUTPUT_PREFIX}_shard${shard}"

    bsub <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_${pilot}_rfd_${shard}
#BSUB -q ${QUEUE}
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -gpu "num=1/host"
#BSUB -o ${WORKDIR}/${target_lower}_${pilot}_rfd_${shard}.%J.out
#BSUB -e ${WORKDIR}/${target_lower}_${pilot}_rfd_${shard}.%J.err

set -euo pipefail
date

module load cuda/11.8

MICROMAMBA="${MICROMAMBA}"
ENV_NAME="${ENV_NAME}"
RFDIFFUSION_DIR="${RFDIFFUSION_DIR}"
INPUT_PDB="${INPUT_PDB}"
SHARD_PREFIX="${shard_prefix}"

mkdir -p "\$(dirname "\$SHARD_PREFIX")"
export HYDRA_FULL_ERROR=1

"\$MICROMAMBA" run -n "\$ENV_NAME" python "\$RFDIFFUSION_DIR/scripts/run_inference.py" \
  --config-name base \
  inference.output_prefix="\$SHARD_PREFIX" \
  inference.num_designs=${DESIGNS_PER_SHARD} \
  'contigmap.contigs=[${contigs}]' \
  inference.input_pdb="\$INPUT_PDB" \
  inference.cyclic=True \
  diffuser.T=${rf_diffuser_T} \
  inference.cyc_chains='a'

date
EOF

done

echo "Submitted ${N_SHARDS} ${target} RFdiffusion jobs for ${pilot}, ${DESIGNS_PER_SHARD} designs each."
"""


PROTEINMPNN_TEMPLATE = r"""#!/bin/bash
set -euo pipefail

QUEUE="${QUEUE:-${mpnn_queue}}"
N_SHARDS="${N_SHARDS:-${mpnn_n_shards}}"

MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME=${mpnn_env}
SCRIPT=${mpnn_script}
PYTHON_BIN="$($MICROMAMBA run -n "$ENV_NAME" which python)"

INPUT_DIR=${scratch_root}/example_outputs/${pilot}
OUT_DIR=${scratch_root}/mpnn_relax${mpnn_relax_cycles}_out/${pilot}

RUNLIST_DIR="$OUT_DIR/runlists"
CHECKPOINT_DIR="$OUT_DIR/checkpoints"
LOG_DIR="$OUT_DIR/logs"
JOB_TMP_ROOT="$OUT_DIR/job_workdirs"

mkdir -p "$OUT_DIR" "$RUNLIST_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$JOB_TMP_ROOT"
ALL_TAGS="$RUNLIST_DIR/all_tags.txt"

python3 - <<'PY' "$INPUT_DIR" "$ALL_TAGS"
import sys, os, glob
input_dir, out_file = sys.argv[1], sys.argv[2]
tags = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(input_dir, "*.pdb")))
if not tags:
    raise SystemExit(f"No pdb files found in {input_dir}")
with open(out_file, "w") as handle:
    for tag in tags:
        handle.write(tag + "\n")
print(f"Wrote {len(tags)} tags to {out_file}")
PY

find "$RUNLIST_DIR" -maxdepth 1 -type f -name 'runlist_*.txt' -delete

python3 - <<'PY' "$ALL_TAGS" "$N_SHARDS" "$RUNLIST_DIR"
import sys, os, math
all_tags_file, n_shards, runlist_dir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
with open(all_tags_file) as handle:
    tags = [line.strip() for line in handle if line.strip()]
chunk = math.ceil(len(tags) / n_shards)
for i in range(n_shards):
    subset = tags[i * chunk:(i + 1) * chunk]
    if not subset:
        continue
    path = os.path.join(runlist_dir, f"runlist_{i:02d}.txt")
    with open(path, "w") as handle:
        for tag in subset:
            handle.write(tag + "\n")
    print(f"{path}: {len(subset)} tags")
PY

for RUNLIST in "$RUNLIST_DIR"/runlist_*.txt; do
    shard=$(basename "$RUNLIST" .txt)
    CHECKPOINT="$CHECKPOINT_DIR/${shard}.check.point"
    JOB_TMP="$JOB_TMP_ROOT/${shard}"

    bsub <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_mpnn_${shard}
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -o ${LOG_DIR}/${shard}.%J.out
#BSUB -e ${LOG_DIR}/${shard}.%J.err

set -euo pipefail
date

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

rm -f "$CHECKPOINT"
mkdir -p "$JOB_TMP"
cd "$JOB_TMP"
rm -f temp.pdb

"$PYTHON_BIN" "$SCRIPT" \
  -pdbdir "$INPUT_DIR" \
  -outpdbdir "$OUT_DIR" \
  -runlist "$RUNLIST" \
  -checkpoint_name "$CHECKPOINT" \
  -relax_cycles ${mpnn_relax_cycles} \
  -seqs_per_struct ${mpnn_seqs_per_struct} \
  -output_intermediates

date
EOF
done

echo "Submitted ${N_SHARDS} ${target} ProteinMPNN shards to ${QUEUE}."
"""


AFCYC_TEMPLATE = r"""#!/bin/bash
set -euo pipefail

GPU_QUEUE="${GPU_QUEUE:-${afcyc_gpu_queue}}"
GPU_NCPU="${GPU_NCPU:-1}"
GPU_PTILE="${GPU_PTILE:-1}"
GPU_REQ='num=1/host'

CPU_QUEUE="${CPU_QUEUE:-${afcyc_cpu_queue}}"
CPU_NCPU="${CPU_NCPU:-1}"
CPU_PTILE="${CPU_PTILE:-1}"

N_SHARDS="${N_SHARDS:-${afcyc_n_shards}}"
MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME=${afcyc_env}
AFCYC_SCRIPT=${afcyc_script}
RMSD_SCRIPT=${afcyc_rmsd_script}
MERGE_SCRIPT=${afcyc_merge_script}

INPUT_DIR=${scratch_root}/mpnn_relax${mpnn_relax_cycles}_out/${pilot}
OUT_BASE=${scratch_root}/afcyc_out/${pilot}
AF_PARAMS_DIR=${afcyc_params}

TARGET_CHAIN=${target_chain}
BINDER_CHAIN=${binder_chain}

RUNLIST_DIR="$OUT_BASE/runlists"
SHARD_INPUT_ROOT="$OUT_BASE/shard_inputs"
SHARD_RESULT_ROOT="$OUT_BASE/shard_results"
LOG_DIR="$OUT_BASE/logs"
MERGED_DIR="$OUT_BASE/merged"

mkdir -p "$RUNLIST_DIR" "$SHARD_INPUT_ROOT" "$SHARD_RESULT_ROOT" "$LOG_DIR" "$MERGED_DIR"
ALL_TAGS="$RUNLIST_DIR/all_tags.txt"

python3 - <<'PY' "$INPUT_DIR" "$ALL_TAGS"
import sys, os, glob
input_dir, out_file = sys.argv[1], sys.argv[2]
pdbs = sorted(glob.glob(os.path.join(input_dir, "*.pdb")))
if not pdbs:
    raise SystemExit(f"No pdb files found in {input_dir}")
with open(out_file, "w") as handle:
    for pdb in pdbs:
        handle.write(os.path.basename(pdb) + "\n")
print(f"Wrote {len(pdbs)} pdb names to {out_file}")
PY

python3 - <<'PY' "$ALL_TAGS" "$RUNLIST_DIR" "$N_SHARDS"
import sys, os, math
all_tags_file, runlist_dir, n_shards = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(all_tags_file) as handle:
    names = [line.strip() for line in handle if line.strip()]
chunk = math.ceil(len(names) / n_shards)
for filename in os.listdir(runlist_dir):
    if filename.startswith("runlist_") and filename.endswith(".txt"):
        os.remove(os.path.join(runlist_dir, filename))
for i in range(n_shards):
    subset = names[i * chunk:(i + 1) * chunk]
    if not subset:
        continue
    path = os.path.join(runlist_dir, f"runlist_{i:02d}.txt")
    with open(path, "w") as handle:
        for name in subset:
            handle.write(name + "\n")
    print(f"{path}: {len(subset)} files")
PY

for RUNLIST in "$RUNLIST_DIR"/runlist_*.txt; do
    shard=$(basename "$RUNLIST" .txt)
    SHARD_INPUT_DIR="$SHARD_INPUT_ROOT/$shard"
    mkdir -p "$SHARD_INPUT_DIR"
    find "$SHARD_INPUT_DIR" -maxdepth 1 -type l -name "*.pdb" -delete
    while IFS= read -r pdb_name; do
        ln -sf "$INPUT_DIR/$pdb_name" "$SHARD_INPUT_DIR/$pdb_name"
    done < "$RUNLIST"
done

AFCYC_JOB_IDS=()
RMSD_JOB_IDS=()

for RUNLIST in "$RUNLIST_DIR"/runlist_*.txt; do
    shard=$(basename "$RUNLIST" .txt)
    SHARD_INPUT_DIR="$SHARD_INPUT_ROOT/$shard"
    SHARD_OUT_DIR="$SHARD_RESULT_ROOT/$shard"
    mkdir -p "$SHARD_OUT_DIR"

    af_submit=$(bsub <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_afcyc_${shard}
#BSUB -q ${GPU_QUEUE}
#BSUB -n ${GPU_NCPU}
#BSUB -R "span[ptile=${GPU_PTILE}]"
#BSUB -gpu "${GPU_REQ}"
#BSUB -o ${LOG_DIR}/${shard}.afcyc.%J.out
#BSUB -e ${LOG_DIR}/${shard}.afcyc.%J.err
set -euo pipefail
date
"${MICROMAMBA}" run -n "${ENV_NAME}" python "${AFCYC_SCRIPT}" \
  --input_dir "${SHARD_INPUT_DIR}" \
  --out_dir "${SHARD_OUT_DIR}" \
  --data_dir "${AF_PARAMS_DIR}" \
  --target_chain "${TARGET_CHAIN}" \
  --binder_chain "${BINDER_CHAIN}" \
  --offset_type 2 \
  --num_recycles ${afcyc_num_recycles} \
  --num_models ${afcyc_num_models}
date
EOF
)
    echo "$af_submit"
    af_job_id=$(echo "$af_submit" | sed -n 's/Job <\([0-9]\+\)>.*/\1/p')
    [[ -n "$af_job_id" ]] && AFCYC_JOB_IDS+=("$af_job_id")

    rmsd_submit=$(bsub -w "done(${af_job_id})" <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_rmsd_${shard}
#BSUB -q ${CPU_QUEUE}
#BSUB -n ${CPU_NCPU}
#BSUB -R "span[ptile=${CPU_PTILE}]"
#BSUB -o ${LOG_DIR}/${shard}.rmsd.%J.out
#BSUB -e ${LOG_DIR}/${shard}.rmsd.%J.err
set -euo pipefail
date
"${MICROMAMBA}" run -n "${ENV_NAME}" python "${RMSD_SCRIPT}" \
  --design_dir "${SHARD_INPUT_DIR}" \
  --pred_dir "${SHARD_OUT_DIR}/pred_pdbs" \
  --out_csv "${SHARD_OUT_DIR}/rmsd_results.csv" \
  --target_chain "${TARGET_CHAIN}" \
  --binder_chain "${BINDER_CHAIN}"
date
EOF
)
    echo "$rmsd_submit"
    rmsd_job_id=$(echo "$rmsd_submit" | sed -n 's/Job <\([0-9]\+\)>.*/\1/p')
    [[ -n "$rmsd_job_id" ]] && RMSD_JOB_IDS+=("$rmsd_job_id")
done

if [[ ${#RMSD_JOB_IDS[@]} -gt 0 ]]; then
    DEP=$(printf "done(%s) && " "${RMSD_JOB_IDS[@]}")
    DEP=${DEP% && }
    MERGED_CSV="$MERGED_DIR/results_merged.csv"

    bsub -w "$DEP" <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_afcyc_merge
#BSUB -q ${CPU_QUEUE}
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -o ${LOG_DIR}/merge.%J.out
#BSUB -e ${LOG_DIR}/merge.%J.err
set -euo pipefail
date
"${MICROMAMBA}" run -n "${ENV_NAME}" python "${MERGE_SCRIPT}" \
  --shard_root "${SHARD_RESULT_ROOT}" \
  --out_csv "${MERGED_CSV}"
date
EOF
fi

echo "Submitted ${#AFCYC_JOB_IDS[@]} ${target} AfCyc shard jobs."
echo "Submitted ${#RMSD_JOB_IDS[@]} ${target} RMSD shard jobs."
echo "Merged csv: $MERGED_DIR/results_merged.csv"
"""


PYROSETTA_TEMPLATE = r"""#!/bin/bash
set -euo pipefail

QUEUE="${QUEUE:-${pyro_queue}}"
N_SHARDS="${N_SHARDS:-${pyro_n_shards}}"

MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME=${pyro_env}
PYTHON_BIN="$($MICROMAMBA run -n "$ENV_NAME" which python)"
SCRIPT=${pyro_script}
MERGE_SCRIPT=${pyro_merge_script}

INPUT_DIR=${scratch_root}/mpnn_relax${mpnn_relax_cycles}_out/${pilot}
OUT_BASE=${scratch_root}/pyrosetta_scores/${pilot}

TARGET_CHAIN=${target_chain}
BINDER_CHAIN=${binder_chain}

NCPU=${pyro_ncpu}
PTILE=${pyro_ptile}

RUNLIST_DIR="$OUT_BASE/runlists"
SHARD_INPUT_ROOT="$OUT_BASE/shard_inputs"
SHARD_RESULT_ROOT="$OUT_BASE/shard_results"
LOG_DIR="$OUT_BASE/logs"
MERGED_DIR="$OUT_BASE/merged"

mkdir -p "$RUNLIST_DIR" "$SHARD_INPUT_ROOT" "$SHARD_RESULT_ROOT" "$LOG_DIR" "$MERGED_DIR"
ALL_TAGS="$RUNLIST_DIR/all_pdbs.txt"

python3 - <<'PY' "$INPUT_DIR" "$ALL_TAGS"
import sys, os, glob
input_dir, out_file = sys.argv[1], sys.argv[2]
pdbs = sorted(glob.glob(os.path.join(input_dir, "*.pdb")))
if not pdbs:
    raise SystemExit(f"No pdb files found in {input_dir}")
with open(out_file, "w") as handle:
    for pdb in pdbs:
        handle.write(os.path.basename(pdb) + "\n")
print(f"Wrote {len(pdbs)} pdb names to {out_file}")
PY

python3 - <<'PY' "$ALL_TAGS" "$RUNLIST_DIR" "$N_SHARDS"
import sys, os, math
all_tags_file, runlist_dir, n_shards = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(all_tags_file) as handle:
    names = [line.strip() for line in handle if line.strip()]
chunk = math.ceil(len(names) / n_shards)
for filename in os.listdir(runlist_dir):
    if filename.startswith("runlist_") and filename.endswith(".txt"):
        os.remove(os.path.join(runlist_dir, filename))
for i in range(n_shards):
    subset = names[i * chunk:(i + 1) * chunk]
    if not subset:
        continue
    path = os.path.join(runlist_dir, f"runlist_{i:02d}.txt")
    with open(path, "w") as handle:
        for name in subset:
            handle.write(name + "\n")
    print(f"{path}: {len(subset)} files")
PY

for RUNLIST in "$RUNLIST_DIR"/runlist_*.txt; do
    shard=$(basename "$RUNLIST" .txt)
    SHARD_INPUT_DIR="$SHARD_INPUT_ROOT/$shard"
    mkdir -p "$SHARD_INPUT_DIR"
    find "$SHARD_INPUT_DIR" -maxdepth 1 -type l -name "*.pdb" -delete
    while IFS= read -r pdb_name; do
        ln -sf "$INPUT_DIR/$pdb_name" "$SHARD_INPUT_DIR/$pdb_name"
    done < "$RUNLIST"
done

JOB_IDS=()

for RUNLIST in "$RUNLIST_DIR"/runlist_*.txt; do
    shard=$(basename "$RUNLIST" .txt)
    SHARD_INPUT_DIR="$SHARD_INPUT_ROOT/$shard"
    SHARD_OUT_DIR="$SHARD_RESULT_ROOT/$shard"
    SHARD_OUT_CSV="$SHARD_OUT_DIR/pyrosetta_scores.csv"
    mkdir -p "$SHARD_OUT_DIR"

    job_submit_output=$(bsub <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_pyro_${shard}
#BSUB -q ${QUEUE}
#BSUB -n ${NCPU}
#BSUB -R "span[ptile=${PTILE}]"
#BSUB -o ${LOG_DIR}/${shard}.%J.out
#BSUB -e ${LOG_DIR}/${shard}.%J.err

set -euo pipefail
date

export OMP_NUM_THREADS=${NCPU}
export OPENBLAS_NUM_THREADS=${NCPU}
export MKL_NUM_THREADS=${NCPU}
export NUMEXPR_NUM_THREADS=${NCPU}

"$PYTHON_BIN" "$SCRIPT" \
  --input_dir "$SHARD_INPUT_DIR" \
  --out_csv "$SHARD_OUT_CSV" \
  --target_chain "$TARGET_CHAIN" \
  --binder_chain "$BINDER_CHAIN" \
${pyro_flags}

date
EOF
)

    echo "$job_submit_output"
    job_id=$(echo "$job_submit_output" | sed -n 's/Job <\([0-9]\+\)>.*/\1/p')
    [[ -n "$job_id" ]] && JOB_IDS+=("$job_id")
done

if [[ ${#JOB_IDS[@]} -gt 0 ]]; then
    DEP=$(printf "done(%s) && " "${JOB_IDS[@]}")
    DEP=${DEP% && }
    MERGED_CSV="$MERGED_DIR/pyrosetta_scores_merged.csv"

    bsub -w "$DEP" <<EOF
#!/bin/bash
#BSUB -J ${target_lower}_pyro_merge
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -o ${LOG_DIR}/merge.%J.out
#BSUB -e ${LOG_DIR}/merge.%J.err

set -euo pipefail
date

"$PYTHON_BIN" "$MERGE_SCRIPT" \
  --shard_root "$SHARD_RESULT_ROOT" \
  --out_csv "$MERGED_CSV"

date
EOF
fi

echo "Submitted ${#JOB_IDS[@]} ${target} PyRosetta shard jobs."
echo "Merged csv: $MERGED_DIR/pyrosetta_scores_merged.csv"
"""


def pyro_flags(config: dict) -> str:
    settings = config["pyrosetta"]
    flags = []
    if settings.get("pack_input", True):
        flags.append("  --pack_input \\")
    if settings.get("pack_separated", True):
        flags.append("  --pack_separated \\")
    if settings.get("packstat", True):
        flags.append("  --packstat")
    if not flags:
        return "  "
    flags[-1] = flags[-1].rstrip(" \\")
    return "\n".join(flags)


def build_values(config: dict) -> dict:
    values = common_values(config)
    rf = config["rfdiffusion"]
    mpnn = config["proteinmpnn"]
    afcyc = config["afcyc"]
    pyro = config["pyrosetta"]
    values.update({
        "rf_queue": sh_quote(rf.get("queue", "8v100-32-sc")),
        "rf_env": sh_quote(rf.get("env_name", "SE3nv")),
        "rf_dir": sh_quote(rf.get("rfdiffusion_dir", "$HOME/RFdiffusion")),
        "rf_n_shards": rf.get("n_shards", 10),
        "rf_designs_per_shard": rf.get("designs_per_shard", 20),
        "rf_diffuser_T": rf.get("diffuser_T", 50),
        "mpnn_queue": mpnn.get("queue", "33"),
        "mpnn_env": sh_quote(mpnn.get("env_name", "proteinmpnn_binder_design")),
        "mpnn_script": sh_quote(mpnn.get("script", "$HOME/dl_binder_design/mpnn_fr/dl_interface_design.py")),
        "mpnn_n_shards": mpnn.get("n_shards", 20),
        "mpnn_relax_cycles": mpnn.get("relax_cycles", 4),
        "mpnn_seqs_per_struct": mpnn.get("seqs_per_struct", 1),
        "afcyc_gpu_queue": afcyc.get("gpu_queue", "8v100-32-sc"),
        "afcyc_cpu_queue": afcyc.get("cpu_queue", "33"),
        "afcyc_env": sh_quote(afcyc.get("env_name", "afcycdesign")),
        "afcyc_script": sh_quote(afcyc.get("afcyc_script", f"{config['home_project_dir'].rstrip('/')}/3.AfCycDesign/v3/afcyc_predict_batch.py")),
        "afcyc_rmsd_script": sh_quote(afcyc.get("rmsd_script", f"{config['home_project_dir'].rstrip('/')}/3.AfCycDesign/v3/rmsd_from_afcyc.py")),
        "afcyc_merge_script": sh_quote(afcyc.get("merge_script", f"{config['home_project_dir'].rstrip('/')}/3.AfCycDesign/v3/merge_afcyc_csvs.py")),
        "afcyc_n_shards": afcyc.get("n_shards", 15),
        "afcyc_num_recycles": afcyc.get("num_recycles", 3),
        "afcyc_num_models": afcyc.get("num_models", 1),
        "afcyc_params": sh_quote(afcyc.get("af_params_dir", "$HOME/dl_binder_design/af2_initial_guess/model_weights/params")),
        "pyro_queue": pyro.get("queue", "33"),
        "pyro_env": sh_quote(pyro.get("env_name", "PyRosettaScore")),
        "pyro_script": sh_quote(pyro.get("script", f"{config['home_project_dir'].rstrip('/')}/4.PyRosetta/PyRosetta_fullScoring_v4_debug.py")),
        "pyro_merge_script": sh_quote(pyro.get("merge_script", "$HOME/Exercise_Phase2_design/merge_pyrosetta_csvs.py")),
        "pyro_n_shards": pyro.get("n_shards", 50),
        "pyro_ncpu": pyro.get("ncpu", 2),
        "pyro_ptile": pyro.get("ptile", 2),
        "pyro_flags": pyro_flags(config),
    })
    for key in ("scratch_root", "input_pdb", "target_chain", "binder_chain", "home_project_dir"):
        values[key] = sh_quote(values[key])
    return values


def generate(config: dict, output_root: Path, force: bool) -> list[Path]:
    values = build_values(config)
    target_root = output_root / f"{config['target']}_workflow"
    files = {
        target_root / "1.RFDiffusion" / f"submit_{config['target']}_rfdiffusion_{config.get('pilot', 'pilot0')}.sh": RFDIFFUSION_TEMPLATE,
        target_root / "2.ProteinMPNN" / f"submit_{config['target']}_mpnn_relax{config['proteinmpnn'].get('relax_cycles', 4)}_shards.sh": PROTEINMPNN_TEMPLATE,
        target_root / "3.AfCycDesign" / "submit_afcyc_shards_integrated.sh": AFCYC_TEMPLATE,
        target_root / "4.PyRosetta" / "submit_pyrosetta_shards.sh": PYROSETTA_TEMPLATE,
    }
    written = []
    for path, template in files.items():
        write_file(path, render(template, values), force)
        written.append(path)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate cluster submit scripts for the binder-design workflow.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Target JSON config.")
    parser.add_argument("--out", type=Path, default=SCRIPT_DIR / "generated_examples", help="Output directory.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing generated files.")
    args = parser.parse_args()

    config = load_config(args.config)
    written = generate(config, args.out, args.force)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
