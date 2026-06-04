#!/bin/bash
set -euo pipefail

GPU_QUEUE="${GPU_QUEUE:-8v100-32-sc}"
GPU_NCPU="${GPU_NCPU:-1}"
GPU_PTILE="${GPU_PTILE:-1}"
GPU_REQ='num=1/host'

CPU_QUEUE="${CPU_QUEUE:-33}"
CPU_NCPU="${CPU_NCPU:-1}"
CPU_PTILE="${CPU_PTILE:-1}"

N_SHARDS="${N_SHARDS:-15}"
MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME="afcycdesign"
AFCYC_SCRIPT="$HOME/Test3_PGLYRP1/3.AfCycDesign/v3/afcyc_predict_batch.py"
RMSD_SCRIPT="$HOME/Test3_PGLYRP1/3.AfCycDesign/v3/rmsd_from_afcyc.py"
MERGE_SCRIPT="$HOME/Test3_PGLYRP1/3.AfCycDesign/v3/merge_afcyc_csvs.py"

INPUT_DIR="/scratch/2026-05-24/bme-yaozm/PGLYRP1_test"/mpnn_relax4_out/pilot0
OUT_BASE="/scratch/2026-05-24/bme-yaozm/PGLYRP1_test"/afcyc_out/pilot0
AF_PARAMS_DIR="$HOME/dl_binder_design/af2_initial_guess/model_weights/params"

TARGET_CHAIN="A"
BINDER_CHAIN="B"

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
#BSUB -J pglyrp1_afcyc_${shard}
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
  --num_recycles 3 \
  --num_models 1
date
EOF
)
    echo "$af_submit"
    af_job_id=$(echo "$af_submit" | sed -n 's/Job <\([0-9]\+\)>.*/\1/p')
    [[ -n "$af_job_id" ]] && AFCYC_JOB_IDS+=("$af_job_id")

    rmsd_submit=$(bsub -w "done(${af_job_id})" <<EOF
#!/bin/bash
#BSUB -J pglyrp1_rmsd_${shard}
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
#BSUB -J pglyrp1_afcyc_merge
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

echo "Submitted ${#AFCYC_JOB_IDS[@]} PGLYRP1 AfCyc shard jobs."
echo "Submitted ${#RMSD_JOB_IDS[@]} PGLYRP1 RMSD shard jobs."
echo "Merged csv: $MERGED_DIR/results_merged.csv"
