#!/bin/bash
set -euo pipefail

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUEUE="${QUEUE:-33}"
N_SHARDS="${N_SHARDS:-50}"
WAIT_FOR_STABLE_INPUT="${WAIT_FOR_STABLE_INPUT:-1}"
INPUT_STABILITY_SECONDS="${INPUT_STABILITY_SECONDS:-30}"

MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME="PyRosettaScore"
PYTHON_BIN="$($MICROMAMBA run -n "$ENV_NAME" which python)"
SCRIPT="$WORKFLOW_DIR/scripts/PyRosetta_fullScoring_v4_debug.py"
MERGE_SCRIPT="$WORKFLOW_DIR/scripts/merge_pyrosetta_csvs.py"

INPUT_DIR="/scratch/2026-06-09/bme-yaozm/PGLYRP1_test"/mpnn_relax4_out/pilot0
OUT_BASE="/scratch/2026-06-09/bme-yaozm/PGLYRP1_test"/pyrosetta_scores/pilot0

TARGET_CHAIN="A"
BINDER_CHAIN="B"

NCPU=2
PTILE=2

RUNLIST_DIR="$OUT_BASE/runlists"
SHARD_INPUT_ROOT="$OUT_BASE/shard_inputs"
SHARD_RESULT_ROOT="$OUT_BASE/shard_results"
LOG_DIR="$OUT_BASE/logs"
MERGED_DIR="$OUT_BASE/merged"

mkdir -p "$RUNLIST_DIR" "$SHARD_INPUT_ROOT" "$SHARD_RESULT_ROOT" "$LOG_DIR" "$MERGED_DIR"
ALL_TAGS="$RUNLIST_DIR/all_pdbs.txt"

count_input_pdbs() {
    find "$INPUT_DIR" -maxdepth 1 -type f -name "*.pdb" | wc -l
}

INPUT_COUNT_1=$(count_input_pdbs)
echo "PyRosetta input pdb count before stability check: ${INPUT_COUNT_1}"
if [[ "$INPUT_COUNT_1" -eq 0 ]]; then
    echo "ERROR: no input pdb files found in $INPUT_DIR" >&2
    exit 2
fi

if [[ "$WAIT_FOR_STABLE_INPUT" == "1" ]]; then
    sleep "$INPUT_STABILITY_SECONDS"
    INPUT_COUNT_2=$(count_input_pdbs)
    echo "PyRosetta input pdb count after ${INPUT_STABILITY_SECONDS}s: ${INPUT_COUNT_2}"
    if [[ "$INPUT_COUNT_1" -ne "$INPUT_COUNT_2" ]]; then
        echo "ERROR: input pdb count changed from ${INPUT_COUNT_1} to ${INPUT_COUNT_2}; upstream output is still changing. Wait for upstream jobs to finish, then resubmit PyRosetta." >&2
        exit 3
    fi
fi

if [[ "$INPUT_COUNT_1" -lt "$N_SHARDS" ]]; then
    echo "WARNING: input pdb count (${INPUT_COUNT_1}) is lower than requested shards (${N_SHARDS}); only ${INPUT_COUNT_1} non-empty shard jobs can be submitted." >&2
fi

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
print(f"Created {len([name for name in os.listdir(runlist_dir) if name.startswith('runlist_') and name.endswith('.txt')])} non-empty runlists from {len(names)} pdbs; requested {n_shards} shards")
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
#BSUB -J pglyrp1_pyro_${shard}
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
  --pack_input \
  --pack_separated \
  --packstat

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
#BSUB -J pglyrp1_pyro_merge
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

if [[ ${#JOB_IDS[@]} -ne "$N_SHARDS" ]]; then
    echo "WARNING: submitted ${#JOB_IDS[@]} PyRosetta shard jobs, requested ${N_SHARDS}. Check input pdb count and runlists under $RUNLIST_DIR." >&2
fi
echo "Submitted ${#JOB_IDS[@]} PGLYRP1 PyRosetta shard jobs."
echo "Merged csv: $MERGED_DIR/pyrosetta_scores_merged.csv"
