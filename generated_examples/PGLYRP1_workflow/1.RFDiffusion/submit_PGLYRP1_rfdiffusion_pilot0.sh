#!/bin/bash
set -euo pipefail

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUEUE="hgx-aais-didier"
GPU_NCPU=8
GPU_REQ="num=1:aff=no"
GPU_SPAN="span[ptile=8]"
MICROMAMBA="$HOME/bin/micromamba"
ENV_NAME="SE3nv"
RFDIFFUSION_DIR="$HOME/RFdiffusion"

WORKDIR="/scratch/2026-06-09/bme-yaozm/PGLYRP1_test"
INPUT_PDB="$WORKFLOW_DIR/inputs/1YCK_clean.pdb"
OUTPUT_PREFIX="$WORKDIR/example_outputs/pilot0/diffused_binder_cyclic_PGLYRP1_pilot0"

N_SHARDS=20
DESIGNS_PER_SHARD=20

mkdir -p "$WORKDIR"
mkdir -p "$(dirname "$OUTPUT_PREFIX")"
cd "$WORKDIR"

echo "RFDiffusion GPU resource request: queue=${QUEUE}, ncpu=${GPU_NCPU}, span=${GPU_SPAN}, gpu=${GPU_REQ}"

for i in $(seq 0 $((N_SHARDS-1))); do
    shard=$(printf "%02d" "$i")
    shard_prefix="${OUTPUT_PREFIX}_shard${shard}"

    bsub <<EOF
#!/bin/bash
#BSUB -J pglyrp1_pilot0_rfd_${shard}
#BSUB -q ${QUEUE}
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -gpu "num=1:aff=no"
#BSUB -o ${WORKDIR}/pglyrp1_pilot0_rfd_${shard}.%J.out
#BSUB -e ${WORKDIR}/pglyrp1_pilot0_rfd_${shard}.%J.err

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
  'contigmap.contigs=[A9-175/0 8-16]' \
  inference.input_pdb="\$INPUT_PDB" \
  inference.cyclic=True \
  diffuser.T=50 \
  inference.cyc_chains='a'

date
EOF

done

echo "Submitted ${N_SHARDS} PGLYRP1 RFdiffusion jobs for pilot0, ${DESIGNS_PER_SHARD} designs each."
