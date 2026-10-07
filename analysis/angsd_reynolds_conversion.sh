#!/bin/bash
#SBATCH -J fst_reynolds -A dewoody -p cpu -t 1-00:00:00
#SBATCH --cpus-per-task=16 --mem=128G --array=0-2
#SBATCH -o logs/%x_%A_%a.out -e logs/%x_%A_%a.err
set -euo pipefail
OUT="${CLUSTER_SCRATCH}/GROUSE/nexus/fst/species_LEPC_GRPC_STGR"
PAIRS=("LEPC GRPC" "LEPC STGR" "GRPC STGR")       # same order as the main script
read -r A B <<< "${PAIRS[$SLURM_ARRAY_TASK_ID]}"
P="${OUT}/${A}_${B}"
ml biocontainers angsd/0.940
export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""

[[ $(wc -l < "${P}.2dsfs.ml") -eq 1 ]] || { echo "ERROR: ${P}.2dsfs.ml has >1 line; sum lines first"; exit 1; }
realSFS fst index "${OUT}/${A}.saf.idx" "${OUT}/${B}.saf.idx" -sfs "${P}.2dsfs.ml" \
    -fold 1 -whichFst 0 -P 16 -fstout "${P}.reynolds"
realSFS fst stats "${P}.reynolds.fst.idx" > "${P}.reynolds.fst.global.txt"
realSFS fst stats2 "${P}.reynolds.fst.idx" -win 50000 -step 10000 -type 2 > "${P}.reynolds.fst.windows_50000_10000.txt"
echo "$A x $B  Reynolds unweighted/weighted: $(cat "${P}.reynolds.fst.global.txt")"
