#!/bin/bash
# =============================================================================
# nexus_ngsadmix.sh -- NGSadmix on the genotype-likelihood (beagle) file made
# by nexus_pca_admix.sh, K = K_MIN..K_MAX with REPS independent runs per K.
# Modelled on S. Mathur's admixture.sh (DeWoody lab): same NGSadmix settings
# (-minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9), rewritten for SLURM
# as one array task per (K, replicate) so runs go in parallel.
#
#   bash nexus_ngsadmix.sh submit    # (optional thinning) + array + summary
#   bash nexus_ngsadmix.sh summary   # rebuild the likelihood table only
#
# OUTPUT (${OUT}):
#   r<rep>/nexus_k<K>.{qopt,fopt.gz,log}   NGSadmix output per run
#   loglik_by_run.tsv                      K, rep, log-likelihood
#   evanno.tsv                             mean/SD log-likelihood and delta K
#   best/nexus_k<K>.qopt                   highest-likelihood run per K
#   samples.txt                            sample IDs in .qopt row order
#
# Runs are skipped if their .qopt and .log already exist, so resubmitting
# only reruns what failed.
#
# NOT USED IN THE USFWS REPORT (Objective 2 reports the PCAngsd PCA only).
# =============================================================================
#SBATCH --job-name=nexus_ngsadmix
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail

# =============================================================================
# USER SETTINGS
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/nexus"
PCA_DIR="${PROJECT_DIR}/pca_admix_4.66x"            # from nexus_pca_admix.sh
BEAGLE_IN="${PCA_DIR}/autosomes.beagle.gz"
SAMPLES_IN="${PCA_DIR}/samples.txt"                 # beagle column order
OUT="${PROJECT_DIR}/ngsadmix_4.66x"
K_MIN=1
K_MAX=10
REPS=10
MINMAF=0.05
MAXITER=50000
TOL=1e-9
TOLLIKE50=1e-9
THIN_BP=10000        # 0 = use every SNP; e.g. 10000 keeps >= 10 kb between SNPs
                     #     (less LD among loci, much less memory)
CPUS=20
MAX_PARALLEL=10      # simultaneous NGSadmix runs
# NGSadmix is NOT in the RCAC biocontainers or the angsd module. Build it once
# (commands below) into TOOLS_DIR, or point NGSADMIX at an existing binary.
TOOLS_DIR="${PROJECT_DIR}/tools/NGSadmix"
NGSADMIX="${NGSADMIX:-${TOOLS_DIR}/NGSadmix}"
#   one-time build (login node):
#     mkdir -p ${TOOLS_DIR} && cd ${TOOLS_DIR}
#     wget http://popgen.dk/software/download/NGSadmix/ngsadmix32.cpp
#     ml gcc; g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
#     ./NGSadmix            # prints usage if the build worked
# =============================================================================

mkdir -p "$OUT"/best logs
NK=$((K_MAX - K_MIN + 1)); NTASK=$((NK * REPS))
BEAGLE="$BEAGLE_IN"; [[ "$THIN_BP" -gt 0 ]] && BEAGLE="${OUT}/autosomes.thin${THIN_BP}.beagle.gz"

load_tools() {
    module --force purge 2>/dev/null || true
    ml gcc 2>/dev/null || true     # runtime libstdc++ for the self-built binary
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
}

STAGE="${1:-}"

# -----------------------------------------------------------------------------
# submit (login node)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "submit" ]]; then
    [[ -s "$BEAGLE_IN" ]] || { echo "ERROR: $BEAGLE_IN not found (run nexus_pca_admix.sh first)" >&2; exit 1; }
    # fail now, not in 100 array tasks, if the binary is missing
    [[ -x "$NGSADMIX" ]] || command -v "$NGSADMIX" >/dev/null 2>&1 \
        || { echo "ERROR: NGSadmix not found at $NGSADMIX -- build it first (see USER SETTINGS)" >&2; exit 1; }
    cp -f "$SAMPLES_IN" "${OUT}/samples.txt"
    N=$(wc -l < "${OUT}/samples.txt")

    # optional distance thinning; beagle marker = <contig>_<pos>, and contig
    # names contain underscores, so the position is the text after the LAST one
    if [[ "$THIN_BP" -gt 0 && ! -s "$BEAGLE" ]]; then
        echo "  thinning to >= ${THIN_BP} bp between SNPs ..."
        zcat "$BEAGLE_IN" | awk -v d="$THIN_BP" 'NR==1 {print; next}
            { m=$1; p=m; sub(/.*_/,"",p); c=substr(m,1,length(m)-length(p)-1)
              if (c!=lc || p-lp>=d) {print; lc=c; lp=p} }' | gzip > "${BEAGLE}.tmp"
        mv -f "${BEAGLE}.tmp" "$BEAGLE"
    fi
    M=$(( $(zcat "$BEAGLE" | wc -l) - 1 ))
    NCOL=$( (zcat "$BEAGLE" 2>/dev/null | awk 'NR==1{print (NF-3)/3; exit}') || true )
    [[ "$NCOL" == "$N" ]] || { echo "ERROR: beagle has $NCOL samples, samples.txt has $N" >&2; exit 1; }

    # NGSadmix holds the likelihoods as doubles: ~ 3 x N x M x 8 bytes, plus headroom
    MEM_GB=$(awk -v n="$N" -v m="$M" 'BEGIN{g=3*n*m*8*1.4/1e9 + 4; g=(g<8)?8:int(g+1); print g}')
    echo "  samples: $N   SNPs: $M   K: ${K_MIN}..${K_MAX} x ${REPS} reps = ${NTASK} runs"
    echo "  memory per run: ${MEM_GB}G   (if this exceeds a node, set THIN_BP, e.g. 10000)"
    SELF="$(readlink -f "$0")"
    J1=$(sbatch --parsable --array=0-$((NTASK-1))%${MAX_PARALLEL} \
         --cpus-per-task="$CPUS" --mem="${MEM_GB}G" "$SELF" run)
    J2=$(sbatch --parsable --dependency=afterany:"$J1" --job-name=nexus_ngsadmix_sum \
         --array=0 --cpus-per-task=1 --mem=2G -t 0-00:30:00 "$SELF" summary)
    echo "  submitted: runs $J1 ; summary $J2"
    exit 0
fi

# -----------------------------------------------------------------------------
# run (array task -> K, replicate)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "run" ]]; then
    T=$SLURM_ARRAY_TASK_ID
    K=$(( K_MIN + T / REPS )); R=$(( 1 + T % REPS ))
    D="${OUT}/r${R}"; P="${D}/nexus_k${K}"; mkdir -p "$D"
    if [[ -s "${P}.qopt" && -s "${P}.log" ]]; then echo "K=$K rep=$R done -- skipping"; exit 0; fi
    load_tools
    [[ -x "$NGSADMIX" ]] || command -v "$NGSADMIX" >/dev/null 2>&1 \
        || { echo "ERROR: NGSadmix not found; set NGSADMIX to its full path" >&2; exit 1; }
    SEED=$(( 1000 * K + R ))       # distinct, reproducible seed per run
    echo ">>> K=$K rep=$R seed=$SEED  $(date)"
    "$NGSADMIX" -likes "$BEAGLE" -K "$K" -P "${SLURM_CPUS_PER_TASK:-$CPUS}" \
        -minMaf "$MINMAF" -maxiter "$MAXITER" -tol "$TOL" -tolLike50 "$TOLLIKE50" \
        -seed "$SEED" -outfiles "${P}.tmp"
    for e in qopt fopt.gz filter log; do
        [[ -e "${P}.tmp.${e}" ]] && mv -f "${P}.tmp.${e}" "${P}.${e}"
    done
    echo ">>> K=$K rep=$R: $(grep -o 'best like=[^ ]*' "${P}.log" | tail -1)  $(date)"
    exit 0
fi

# -----------------------------------------------------------------------------
# summary: log-likelihoods, Evanno delta K, best run per K
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "summary" ]]; then
    TAB="${OUT}/loglik_by_run.tsv"
    printf "K\trep\tloglik\n" > "$TAB"
    missing=()
    for K in $(seq "$K_MIN" "$K_MAX"); do
        for R in $(seq 1 "$REPS"); do
            L="${OUT}/r${R}/nexus_k${K}.log"
            LL=$(grep -o 'best like=[^ ]*' "$L" 2>/dev/null | tail -1 | cut -d= -f2 || true)
            if [[ -z "$LL" ]]; then missing+=("K${K}r${R}"); continue; fi
            printf "%s\t%s\t%s\n" "$K" "$R" "$LL" >> "$TAB"
        done
    done
    # Evanno et al. (2005): L'(K)=mean L(K)-mean L(K-1); dK=|L'(K+1)-L'(K)|/sd L(K)
    awk -F'\t' 'NR>1 {n[$1]++; s[$1]+=$3; ss[$1]+=$3*$3; if(!($1 in mx)||$3>mx[$1]) mx[$1]=$3}
        END {
          for (k in n) {m[k]=s[k]/n[k]; v=(n[k]>1)?(ss[k]-n[k]*m[k]^2)/(n[k]-1):0; sd[k]=(v>0)?sqrt(v):0}
          print "K\tn_runs\tmean_loglik\tsd_loglik\tmax_loglik\tdeltaK"
          for (k=1; k<=1000; k++) if (k in n) {
            dk="NA"
            if ((k-1) in n && (k+1) in n && sd[k]>0) dk=sprintf("%.4f", ((m[k+1]-m[k])-(m[k]-m[k-1]))/sd[k])
            if (dk!="NA" && dk<0) dk=sprintf("%.4f",-dk)
            printf "%d\t%d\t%.2f\t%.4f\t%.2f\t%s\n", k, n[k], m[k], sd[k], mx[k], dk }
        }' "$TAB" > "${OUT}/evanno.tsv"
    # copy the highest-likelihood run for each K
    awk -F'\t' 'NR>1 {if(!($1 in b)||$3>b[$1]){b[$1]=$3; r[$1]=$2}} END{for(k in r) print k"\t"r[k]}' "$TAB" |
    while read -r K R; do
        cp -f "${OUT}/r${R}/nexus_k${K}.qopt" "${OUT}/best/nexus_k${K}.qopt"
    done
    column -t "${OUT}/evanno.tsv" 2>/dev/null || cat "${OUT}/evanno.tsv"
    echo "Best run per K copied to ${OUT}/best/ (rows follow ${OUT}/samples.txt)"
    (( ${#missing[@]} )) && echo "WARNING: ${#missing[@]} run(s) missing: ${missing[*]:0:40}" >&2
    exit 0
fi

echo "Usage: bash $0 submit | summary" >&2
exit 1
