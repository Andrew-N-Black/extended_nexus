#!/bin/bash
# =============================================================================
# nexus_pca_admix.sh -- genotype likelihoods (ANGSD beagle), PCA and
# admixture (PCAngsd), and plots, for the depth-harmonized nexus panel.
#
# Stages (chained with SLURM dependencies by "submit"):
#   beagle  array over autosomal region chunks: angsd -doGlf 2 on ALL
#           samples together (one beagle file per chunk)
#   pcangsd one job: concatenate chunk beagles -> PCAngsd covariance +
#           admixture for K = 2..K_MAX -> plot_pca_admix.R
#
# All 506 birds are included (the two STGR x GRPC hybrids are plotted as
# their own group, which is the point of an admixture plot). Autosomes only
# (Z scaffolds excluded). Filters follow the report Methods: -GL 1,
# -minMapQ 30, -minQ 30, -skipTriallelic 1, SNP P < 1e-6, MAF >= 0.01,
# site kept if >= MIN_IND_FRAC of samples have data.
#
# USAGE (login node, from the folder holding this script and plot_pca_admix.R):
#   bash nexus_pca_admix.sh check            # sample count, chunks; no jobs
#   bash nexus_pca_admix.sh submit           # beagle array + pcangsd/plot job
#   bash nexus_pca_admix.sh submit pcangsd   # rerun PCAngsd + plots only
#   Rscript plot_pca_admix.R <OUT_DIR>       # replot anywhere (e.g., laptop)
#
# OUTPUT (${OUT}):
#   samples.txt                     sample IDs in beagle column order
#   autosomes.beagle.gz             all chunks, one header
#   pcangsd/nexus.cov               covariance matrix (PCA)
#   pcangsd/nexus_K<k>.admix.*.Q    admixture proportions per K
#   plots/pca.pdf|png, plots/admixture.pdf|png
# =============================================================================
#SBATCH --job-name=nexus_pca
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail

# =============================================================================
# USER SETTINGS
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/nexus"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
CRAMLIST="${PROJECT_DIR}/crams_4.66x/final_cramlist_4.66x.txt"
POPMAP="${PROJECT_DIR}/popmap_species.txt"         # ID<TAB>SPECIES (header OK)
OUT="${PROJECT_DIR}/pca_admix_4.66x"
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
N_CHUNKS=100
MIN_IND_FRAC=0.80        # report Methods: sites present in >= 80% of samples
MINMAF=0.01
SNP_PVAL=1e-6
MINMAPQ=30
MINQ=30
K_MAX=4                  # admixture K = 2..K_MAX (PCAngsd: K = eigenvectors + 1)
BEAGLE_RES=(8 32G 2-00:00:00)
PCANGSD_RES=(32 250G 3-00:00:00)
MAX_PARALLEL=50
PLOT_R="${SLURM_SUBMIT_DIR:-$PWD}/plot_pca_admix.R"
# =============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUT"/{chunks,beagle,pcangsd,plots} logs
cram_id() { local b; b=$(basename "$1"); echo "${b%%.*}"; }

xalt_fix() {
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
}

build_chunks() {   # autosomal contigs split into N_CHUNKS length-balanced lists
    [[ -s "${OUT}/chunks/chunk_000.rf" ]] && return 0
    [[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: ${REF_FASTA}.fai not found" >&2; exit 1; }
    local c
    for c in ${Z_SCAFFOLDS//,/ }; do
        awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF_FASTA}.fai" \
            || { echo "ERROR: Z scaffold $c not in .fai" >&2; exit 1; }
    done
    awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
        !($1 in Z) {print $1"\t"$2}' "${REF_FASTA}.fai" | sort -k2,2nr > "${OUT}/autosomes.len"
    awk -v n="$N_CHUNKS" -v d="${OUT}/chunks" '
        { b=0; for(i=1;i<n;i++) if (t[i] < t[b]) b=i
          t[b]+=$2; f=sprintf("%s/chunk_%03d.rf", d, b); print $1":" >> f; close(f) }' "${OUT}/autosomes.len"
}
n_chunks_real() { ls "${OUT}"/chunks/chunk_*.rf | wc -l; }

prepare() {
    [[ -s "$CRAMLIST" ]] || { echo "ERROR: $CRAMLIST not found" >&2; exit 1; }
    grep . "$CRAMLIST" | tr -d '\r' > "${OUT}/bamlist.txt"
    while read -r c; do cram_id "$c"; done < "${OUT}/bamlist.txt" > "${OUT}/samples.txt"
    while read -r c; do
        [[ -s "$c" ]] || { echo "ERROR: CRAM missing: $c" >&2; exit 1; }
    done < "${OUT}/bamlist.txt"
    build_chunks
}

STAGE="${1:-}"

# -----------------------------------------------------------------------------
# check / submit (login node)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "check" || "$STAGE" == "submit" ]]; then
    prepare
    N=$(wc -l < "${OUT}/bamlist.txt"); NC=$(n_chunks_real)
    MININD=$(awk -v n="$N" -v f="$MIN_IND_FRAC" 'BEGIN{printf "%d", n*f+0.5}')
    echo "  samples: $N   minInd: $MININD   chunks: $NC   K: 2..$K_MAX"
    echo "  output : $OUT"
    [[ -s "$PLOT_R" ]] || echo "  WARNING: $PLOT_R not found; plots will be skipped (run the R script later)."
    [[ "$STAGE" == "check" ]] && exit 0
    FROM="${2:-beagle}"; SELF="$(readlink -f "$0")"; DEP=""
    if [[ "$FROM" == beagle ]]; then
        DEP=$(sbatch --parsable --job-name=nexus_pca_beagle --array=0-$((NC-1))%${MAX_PARALLEL} \
              --cpus-per-task="${BEAGLE_RES[0]}" --mem="${BEAGLE_RES[1]}" -t "${BEAGLE_RES[2]}" "$SELF" beagle)
        echo "  beagle : job $DEP (array 0-$((NC-1)))"
    fi
    J=$(sbatch --parsable ${DEP:+--dependency=afterok:$DEP} --job-name=nexus_pcangsd --array=0 \
        --cpus-per-task="${PCANGSD_RES[0]}" --mem="${PCANGSD_RES[1]}" -t "${PCANGSD_RES[2]}" "$SELF" pcangsd)
    echo "  pcangsd: job $J"
    exit 0
fi

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "Run with: bash $0 submit   (or: bash $0 check)" >&2; exit 1; }

# -----------------------------------------------------------------------------
# beagle (task = chunk)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "beagle" ]]; then
    CH=$(printf "%03d" "$SLURM_ARRAY_TASK_ID"); B="${OUT}/beagle/chunk_${CH}"
    if [[ -s "${B}.beagle.gz" ]]; then echo "${B}.beagle.gz exists -- skipping"; exit 0; fi
    N=$(wc -l < "${OUT}/bamlist.txt")
    MININD=$(awk -v n="$N" -v f="$MIN_IND_FRAC" 'BEGIN{printf "%d", n*f+0.5}')
    module --force purge 2>/dev/null || true
    ml biocontainers angsd/0.940; xalt_fix
    echo ">>> chunk $CH: $(wc -l < "${OUT}/chunks/chunk_${CH}.rf") contigs, N=$N, minInd=$MININD  $(date)"
    angsd -bam "${OUT}/bamlist.txt" -ref "$REF_FASTA" -rf "${OUT}/chunks/chunk_${CH}.rf" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf "$MINMAF" -SNP_pval "$SNP_PVAL" \
        -skipTriallelic 1 -minMapQ "$MINMAPQ" -minQ "$MINQ" \
        -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -minInd "$MININD" -P "$THREADS" -out "${B}.tmp"
    mv -f "${B}.tmp.mafs.gz" "${B}.mafs.gz"; mv -f "${B}.tmp.arg" "${B}.arg"
    mv -f "${B}.tmp.beagle.gz" "${B}.beagle.gz"
    echo ">>> chunk $CH: $(( $(zcat "${B}.beagle.gz" | wc -l) - 1 )) SNPs  $(date)"
    exit 0
fi

# -----------------------------------------------------------------------------
# pcangsd: merge beagles, PCA + admixture, plots
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "pcangsd" ]]; then
    NC=$(n_chunks_real)
    mapfile -t PARTS < <(ls "${OUT}"/beagle/chunk_*.beagle.gz 2>/dev/null | sort)
    (( ${#PARTS[@]} == NC )) || { echo "ERROR: ${#PARTS[@]} of $NC chunk beagles present" >&2; exit 1; }
    BEAGLE="${OUT}/autosomes.beagle.gz"
    if [[ ! -s "$BEAGLE" ]]; then
        echo ">>> concatenating $NC chunk beagles  $(date)"
        { zcat "${PARTS[0]}" | head -n 1
          for f in "${PARTS[@]}"; do zcat "$f" | tail -n +2; done; } | gzip > "${BEAGLE}.tmp"
        mv -f "${BEAGLE}.tmp" "$BEAGLE"
    fi
    NSNP=$(( $(zcat "$BEAGLE" | wc -l) - 1 )); echo ">>> $NSNP autosomal SNPs"
    # beagle columns must match samples.txt (3 GL columns per sample)
    NCOL=$(zcat "$BEAGLE" | head -n 1 | awk '{print (NF-3)/3}')
    [[ "$NCOL" == "$(wc -l < "${OUT}/samples.txt")" ]] \
        || { echo "ERROR: beagle has $NCOL samples, samples.txt has $(wc -l < "${OUT}/samples.txt")" >&2; exit 1; }

    module --force purge 2>/dev/null || true
    ml biocontainers pcangsd; xalt_fix
    # PCAngsd >= 1.0 uses --maf; some older builds used --minMaf
    if pcangsd --help 2>&1 | grep -q -- '--minMaf'; then MAFOPT="--minMaf"; else MAFOPT="--maf"; fi
    PC="${OUT}/pcangsd/nexus"
    if [[ ! -s "${PC}.cov" ]]; then
        pcangsd -b "$BEAGLE" -o "$PC" -t "$THREADS" $MAFOPT "$MINMAF"
    fi
    for K in $(seq 2 "$K_MAX"); do
        compgen -G "${PC}_K${K}.admix*.Q" >/dev/null && { echo "K=$K done -- skipping"; continue; }
        echo ">>> admixture K=$K  $(date)"
        pcangsd -b "$BEAGLE" -o "${PC}_K${K}" -t "$THREADS" $MAFOPT "$MINMAF" -e $((K-1)) --admix
    done

    # sample metadata for plotting: ID, species (from popmap), in beagle order
    awk -F'\t' 'NR==FNR {gsub(/\r/,""); k=$1; gsub(/[ \t]+/,"",k); sub(/^normal_/,"",k);
                         s=$0; sub(/^[^\t]*\t+/,"",s); gsub(/^[ \t]+|[ \t]+$/,"",s); sp[k]=s; next}
                {print $1"\t"(($1 in sp)?sp[$1]:"unassigned")}' "$POPMAP" "${OUT}/samples.txt" \
        | awk 'BEGIN{print "sample_id\tspecies"} {print}' > "${OUT}/samples_species.tsv"

    module --force purge 2>/dev/null || true
    ml r 2>/dev/null || ml R 2>/dev/null || true
    if command -v Rscript >/dev/null && [[ -s "$PLOT_R" ]]; then
        Rscript "$PLOT_R" "$OUT"
    else
        echo "Plots skipped: copy ${OUT}/{samples_species.tsv,pcangsd/} and run 'Rscript plot_pca_admix.R <dir>'." >&2
    fi
    echo ">>> done $(date)"
    exit 0
fi

echo "Unknown stage '$STAGE' (use: check | submit [beagle|pcangsd])" >&2
exit 1
