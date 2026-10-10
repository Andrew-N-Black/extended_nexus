#!/bin/bash
# =============================================================================
# nexus_heterozygosity.sh -- per-sample autosomal heterozygosity (ANGSD +
# realSFS) for the depth-harmonized nexus panel, ending in ONE table:
#
#   ${HET_DIR}/heterozygosity_summary.tsv
#       sample_id  species  heterozygosity  sites  (sites = SFS[0]+SFS[1])
#
# One array task per CRAM, then a summary job (submitted with a dependency)
# that collects every sample into the table above and lists any sample that
# did not finish.
#
# Heterozygosity = SFS[1] / (SFS[0] + SFS[1]) from the folded single-sample
# SFS (the standard ANGSD single-sample estimate). Autosomes only: the two
# Z scaffolds are excluded, because hemizygous females would otherwise show
# artificially low heterozygosity.
#
# USAGE (login node, from the folder holding this script):
#   bash nexus_heterozygosity.sh submit     # array + summary job
#   bash nexus_heterozygosity.sh summary    # rebuild the table only
#
# Reruns are safe: samples with a finished _est.ml are skipped.
# =============================================================================
#SBATCH --job-name=nexus_het
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
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
HET_DIR="${PROJECT_DIR}/heterozygosity_4.66x"
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
MINMAPQ=30
MINQ=30            # base quality (report Methods: -minQ 30)
MIN_DEPTH=3        # per-site minimum depth (-setMinDepth)
CLEANUP_SAF=1      # 1 = delete the per-sample .saf files once the SFS exists
MAX_PARALLEL=50
# =============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-4}
AUTOSOME_RF="${HET_DIR}/autosomes.rf"
mkdir -p "$HET_DIR"/per_sample logs
cram_id() { local b; b=$(basename "$1"); echo "${b%%.*}"; }     # F10.cram -> F10

load_angsd() {
    module --force purge 2>/dev/null || true
    ml biocontainers
    ml angsd/0.940
    # xalt injects LD_PRELOAD into the container (GLIBC_2.33/2.34 errors);
    # blank it inside the container -- a host-side unset alone does not hold
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
}

build_rf() {   # autosomes-only region file, built once and atomically
    [[ -s "$AUTOSOME_RF" ]] && return 0
    [[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: ${REF_FASTA}.fai not found" >&2; exit 1; }
    local c tmp
    for c in ${Z_SCAFFOLDS//,/ }; do
        awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF_FASTA}.fai" \
            || { echo "ERROR: Z scaffold '$c' not in ${REF_FASTA}.fai" >&2; exit 1; }
    done
    tmp=$(mktemp "${AUTOSOME_RF}.XXXXXX")
    awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
        !($1 in Z) {print $1":"}' "${REF_FASTA}.fai" > "$tmp"
    mv -n "$tmp" "$AUTOSOME_RF" 2>/dev/null || true
    rm -f "$tmp"
}

STAGE="${1:-}"

# -----------------------------------------------------------------------------
# submit (login node)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "submit" ]]; then
    [[ -s "$CRAMLIST" ]] || { echo "ERROR: $CRAMLIST not found" >&2; exit 1; }
    build_rf
    N=$(grep -c . "$CRAMLIST")
    echo "  $N CRAMs; autosomes: $(wc -l < "$AUTOSOME_RF") contigs (Z excluded); output: $HET_DIR"
    SELF="$(readlink -f "$0")"
    J1=$(sbatch --parsable --array=0-$((N-1))%${MAX_PARALLEL} "$SELF" sample)
    J2=$(sbatch --parsable --dependency=afterany:"$J1" --job-name=nexus_het_sum \
         --array=0 --cpus-per-task=1 --mem=2G -t 0-00:30:00 "$SELF" summary)
    echo "  submitted: per-sample array $J1 ; summary $J2"
    exit 0
fi

# -----------------------------------------------------------------------------
# sample (array task = line of CRAMLIST)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "sample" ]]; then
    CRAM=$(grep . "$CRAMLIST" | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
    [[ -n "$CRAM" ]] || { echo "ERROR: no CRAM at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE=$(cram_id "$CRAM")
    P="${HET_DIR}/per_sample/${SAMPLE}"
    echo ">>> task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE}  ${CRAM}  $(date)"
    if [[ -s "${P}_est.ml" && -s "${P}_heterozygosity.txt" ]]; then
        echo "already done -- skipping"; exit 0
    fi
    build_rf
    load_angsd

    # 1) site allele frequency likelihoods, autosomes only
    if [[ ! -s "${P}.saf.idx" ]]; then
        angsd -i "$CRAM" -ref "$REF_FASTA" -anc "$REF_FASTA" -rf "$AUTOSOME_RF" \
            -dosaf 1 -GL 1 -minMapQ "$MINMAPQ" -minQ "$MINQ" \
            -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
            -doCounts 1 -setMinDepth "$MIN_DEPTH" \
            -P "$THREADS" -out "${P}.tmp"
        for e in saf.gz saf.pos.gz arg; do mv -f "${P}.tmp.${e}" "${P}.${e}"; done
        mv -f "${P}.tmp.saf.idx" "${P}.saf.idx"          # idx last = complete
    fi

    # 2) folded single-sample SFS
    realSFS "${P}.saf.idx" -P "$THREADS" -fold 1 > "${P}_est.ml.tmp"
    mv -f "${P}_est.ml.tmp" "${P}_est.ml"

    # 3) heterozygosity = SFS[1] / (SFS[0] + SFS[1]); if realSFS printed
    #    more than one line (blocks of sites), the lines are summed
    awk '{a+=$1; b+=$2} END {if (a+b>0) printf "%.8f\t%.0f\n", b/(a+b), a+b; else print "NA\tNA"}' \
        "${P}_est.ml" | awk -v s="$SAMPLE" '{print s"\t"$0}' > "${P}_heterozygosity.txt"
    echo ">>> ${SAMPLE}: H = $(cut -f2 "${P}_heterozygosity.txt")  $(date)"

    [[ "$CLEANUP_SAF" == 1 ]] && rm -f "${P}.saf.gz" "${P}.saf.pos.gz"
    exit 0
fi

# -----------------------------------------------------------------------------
# summary: one table, in CRAM-list order, with species
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "summary" ]]; then
    OUTF="${HET_DIR}/heterozygosity_summary.tsv"
    printf "sample_id\tspecies\theterozygosity\tsites\n" > "${OUTF}.tmp"
    missing=()
    while read -r cram; do
        [[ -z "$cram" ]] && continue
        S=$(cram_id "$cram"); F="${HET_DIR}/per_sample/${S}_heterozygosity.txt"
        if [[ ! -s "$F" ]]; then missing+=("$S"); continue; fi
        SP=$(awk -F'\t' -v id="$S" '{gsub(/\r/,""); k=$1; gsub(/[ \t]+/,"",k); sub(/^normal_/,"",k)}
                 k==id {s=$0; sub(/^[^\t]*\t+/,"",s); gsub(/^[ \t]+|[ \t]+$/,"",s); print s; exit}' "$POPMAP" 2>/dev/null || true)
        read -r _ H SITES < "$F"
        printf "%s\t%s\t%s\t%s\n" "$S" "${SP:-unassigned}" "$H" "$SITES" >> "${OUTF}.tmp"
    done < "$CRAMLIST"
    mv -f "${OUTF}.tmp" "$OUTF"
    echo "Wrote $OUTF ($(($(wc -l < "$OUTF")-1)) samples)"
    awk -F'\t' 'NR>1 && $3!="NA" {n[$2]++; s[$2]+=$3} END {for (k in n) printf "  %-12s n=%-4d mean H = %.6f\n", k, n[k], s[k]/n[k]}' "$OUTF"
    if (( ${#missing[@]} )); then
        echo "WARNING: ${#missing[@]} sample(s) not finished: ${missing[*]:0:40}" >&2
        echo "         rerun 'bash $0 submit' -- finished samples are skipped." >&2
    fi
    exit 0
fi

echo "Usage: bash $0 submit | summary" >&2
exit 1
