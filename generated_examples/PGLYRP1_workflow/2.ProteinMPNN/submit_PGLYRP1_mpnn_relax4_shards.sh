#!/bin/bash
set -euo pipefail

QUEUE="${QUEUE:-33}"
N_SHARDS="${N_SHARDS:-20}"

MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME="proteinmpnn_binder_design"
SCRIPT="$HOME/dl_binder_design/mpnn_fr/dl_interface_design.py"
PYTHON_BIN="$($MICROMAMBA run -n "$ENV_NAME" which python)"

INPUT_DIR="/scratch/2026-05-24/bme-yaozm/PGLYRP1_test"/example_outputs/pilot0
OUT_DIR="/scratch/2026-05-24/bme-yaozm/PGLYRP1_test"/mpnn_relax4_out/pilot0

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
#BSUB -J pglyrp1_mpnn_${shard}
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
  -relax_cycles 4 \
  -seqs_per_struct 1 \
  -output_intermediates

date
EOF
done

echo "Submitted ${N_SHARDS} PGLYRP1 ProteinMPNN shards to ${QUEUE}."
