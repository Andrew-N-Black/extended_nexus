#!/bin/bash
# =============================================================================
# nexus_downsample_depth_mapping.sh
#
# For EVERY CRAM in CRAM_DIR (one SLURM array task per sample):
#   1. measure genome-wide mean depth            (samtools coverage)
#   2. if depth > TARGET_DEPTH x (1 + SKIP_TOL): subsample reads to
#      TARGET_DEPTH with a fixed seed and write a new CRAM + index;
#      otherwise link the original unchanged
#   3. on the FINAL CRAM: per-contig depth (samtools coverage) and
#      read/mapping counts (samtools flagstat)
# Then one summary job builds a per-sample table, a QC flag column, the new
# CRAM list, and plots by species (plot_depth_mapping.R).
#
#   bash nexus_downsample_depth_mapping.sh submit    # array + summary job
#   bash nexus_downsample_depth_mapping.sh summary   # rebuild table/plots only
#
# Originals are never modified. Rerunning "submit" skips finished samples.
#
# Why samtools and not mosdepth: the only mosdepth on Negishi (0.3.3) cannot
# decode these CRAM 3.1 files ("Slice decode failure"). samtools coverage
# applies the same default read filters (unmapped, secondary, QC-fail and
# duplicate reads excluded) and gives the same mean depth.
#
# Depth targeting uses the GENOME-WIDE mean (all contigs), as in the original
# 4.66x subsampling; the table also reports autosomal and Z means.
# =============================================================================
#SBATCH --job-name=nexus_ds
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail

# =============================================================================
# USER SETTINGS
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/nexus"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
CRAM_DIR="${PROJECT_DIR}/crams"                   # every *.cram here is processed
POPMAP="${PROJECT_DIR}/popmap_species.txt"        # ID<TAB>SPECIES (header OK)
OUT="${PROJECT_DIR}/crams_4.66x"                  # new CRAMs / links + QC
TARGET_DEPTH=4.66
SKIP_TOL=0.10      # leave a sample unchanged if depth <= target x 1.10
QC_TOL=0.15        # flag final depth more than 15% from target
SEED=42
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
PLOT_R="${SLURM_SUBMIT_DIR:-$PWD}/plot_depth_mapping.R"
MAX_PARALLEL=60
# =============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-4}
QC="${OUT}/qc"
mkdir -p "$OUT" "$QC"/per_sample logs
INLIST="${OUT}/input_cramlist.txt"
cram_id() { local b; b=$(basename "$1"); echo "${b%%.*}"; }

load_tools() {
    module --force purge 2>/dev/null || true
    ml biocontainers samtools
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
}

# per-contig coverage in a simple table: chrom length bases mean
coverage_table() {   # coverage_table <cram> <out>
    samtools coverage --input-fmt-option reference="$REF_FASTA" "$1" \
      | awk 'BEGIN{OFS="\t"; print "chrom","length","bases","mean"}
             /^#/ {next}
             {len=$3-$2+1; b=$7*len; L+=len; B+=b; print $1,len,sprintf("%.0f",b),$7}
             END {print "total",L,sprintf("%.0f",B),(L?B/L:0)}' > "$2.tmp"
    mv -f "$2.tmp" "$2"
}
total_depth() { awk '$1=="total"{print $4}' "$1"; }

STAGE="${1:-}"

# -----------------------------------------------------------------------------
# submit (login node)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "submit" ]]; then
    ls -1 "$CRAM_DIR"/*.cram | sort -V > "${INLIST}.tmp"
    if [[ -s "$INLIST" ]] && ! cmp -s "$INLIST" "${INLIST}.tmp"; then
        echo "NOTE: the set of CRAMs in $CRAM_DIR changed since the last run; list updated." >&2
    fi
    mv -f "${INLIST}.tmp" "$INLIST"
    N=$(wc -l < "$INLIST")
    [[ -s "$PLOT_R" ]] || { echo "ERROR: $PLOT_R not found (keep it next to this script)" >&2; exit 1; }
    SELF="$(readlink -f "$0")"
    J1=$(sbatch --parsable --array=0-$((N-1))%${MAX_PARALLEL} "$SELF" sample)
    J2=$(sbatch --parsable --dependency=afterany:"$J1" --job-name=nexus_ds_sum \
         --array=0 --cpus-per-task=1 --mem=4G -t 0-01:00:00 "$SELF" summary)
    echo "Submitted $N per-sample tasks: $J1 ; summary: $J2"
    echo "Output: $OUT"
    exit 0
fi

# -----------------------------------------------------------------------------
# sample (array task = line of input_cramlist.txt)
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "sample" ]]; then
    IN=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$INLIST")
    ID=$(cram_id "$IN")
    P="${QC}/per_sample/${ID}"
    FINAL="${OUT}/${ID}.cram"
    echo ">>> $ID  $IN  $(date)"
    if [[ -s "${P}.stats.tsv" && -s "${P}.final.coverage.txt" && -s "${P}.flagstat.tsv" ]]; then
        echo "already done -- skipping"; exit 0
    fi
    load_tools

    # 1) depth of the input CRAM
    [[ -s "${P}.input.coverage.txt" ]] || coverage_table "$IN" "${P}.input.coverage.txt"
    D_IN=$(total_depth "${P}.input.coverage.txt")
    FRAC=$(awk -v t="$TARGET_DEPTH" -v d="$D_IN" 'BEGIN{f=(d>0)?t/d:1; if(f>1)f=1; printf "%.4f", f}')
    DO_DS=$(awk -v d="$D_IN" -v t="$TARGET_DEPTH" -v k="$SKIP_TOL" 'BEGIN{print (d > t*(1+k)) ? 1 : 0}')
    echo "    input depth ${D_IN}x -> fraction ${FRAC} (subsample: ${DO_DS})"

    # 2) subsample, or link the original
    rm -f "$FINAL" "${FINAL}.crai"
    if [[ "$DO_DS" == 1 ]]; then
        # -s SEED.FRACTION: same seed for every sample; reads kept by name
        # hash, so both mates of a pair are kept or dropped together
        SEEDFRAC="${SEED}.${FRAC#0.}"
        samtools view -@ "$THREADS" -T "$REF_FASTA" -s "$SEEDFRAC" -C \
            -o "${FINAL}.tmp.cram" "$IN"
        mv -f "${FINAL}.tmp.cram" "$FINAL"
        samtools index -@ "$THREADS" "$FINAL"
        STATUS=subsampled
    else
        ln -s "$(readlink -f "$IN")" "$FINAL"
        if   [[ -s "${IN}.crai" ]];          then ln -s "$(readlink -f "${IN}.crai")" "${FINAL}.crai"
        elif [[ -s "${IN%.cram}.crai" ]];    then ln -s "$(readlink -f "${IN%.cram}.crai")" "${FINAL}.crai"
        else samtools index "$FINAL" "${FINAL}.crai"; fi
        SEEDFRAC=NA; FRAC=1.0000
        STATUS=$(awk -v d="$D_IN" -v t="$TARGET_DEPTH" 'BEGIN{print (d < t) ? "kept_below_target" : "kept_near_target"}')
    fi

    # 3) QC on the final CRAM
    if [[ "$STATUS" == subsampled ]]; then
        coverage_table "$FINAL" "${P}.final.coverage.txt"
    else
        cp -f "${P}.input.coverage.txt" "${P}.final.coverage.txt"   # unchanged file
    fi
    samtools flagstat -@ "$THREADS" -O tsv --input-fmt-option reference="$REF_FASTA" "$FINAL" > "${P}.flagstat.tmp"
    mv -f "${P}.flagstat.tmp" "${P}.flagstat.tsv"
    D_OUT=$(total_depth "${P}.final.coverage.txt")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$ID" "$D_IN" "$FRAC" "$SEEDFRAC" "$STATUS" "$D_OUT" > "${P}.stats.tsv"
    echo ">>> $ID: ${D_IN}x -> ${D_OUT}x ($STATUS)  $(date)"
    exit 0
fi

# -----------------------------------------------------------------------------
# summary
# -----------------------------------------------------------------------------
if [[ "$STAGE" == "summary" ]]; then
    TABLE="${QC}/depth_mapping_summary.tsv"
    NEWLIST="${OUT}/final_cramlist_4.66x.txt"
    printf "sample\tspecies\tdepth_before\tsubsample_fraction\tseed_fraction\tstatus\tmean_depth_genome\tmean_depth_autosomal\tmean_depth_Z\tpct_from_target\tdepth_flag\ttotal_reads\tpct_mapped\tpct_both_mates_mapped\tpct_properly_paired\n" > "${TABLE}.tmp"
    : > "${NEWLIST}.tmp"
    missing=()
    while read -r cram; do
        ID=$(cram_id "$cram"); P="${QC}/per_sample/${ID}"
        if [[ ! -s "${P}.stats.tsv" || ! -s "${P}.flagstat.tsv" ]]; then missing+=("$ID"); continue; fi
        SP=$(awk -F'\t' -v id="$ID" '{gsub(/\r/,""); sub(/^normal_/,"",$1)} $1==id {print $2; exit}' "$POPMAP")
        [[ -z "$SP" ]] && SP="unassigned"
        read -r _ D_IN FRAC SEEDFRAC STATUS D_OUT < "${P}.stats.tsv"
        DEPTH=$(awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
            NR==1 {next} $1=="total" {next}
            ($1 in Z) {zl+=$2; zb+=$3; next} {al+=$2; ab+=$3}
            END {printf "%.3f\t%.3f", (al?ab/al:0), (zl?zb/zl:0)}' "${P}.final.coverage.txt")
        QCF=$(awk -v d="$D_OUT" -v t="$TARGET_DEPTH" -v k="$QC_TOL" 'BEGIN{p=100*(d-t)/t;
            f=(p<-100*k)?"LOW":((p>100*k)?"HIGH":"OK"); printf "%.2f\t%s", p, f}')
        MAP=$(awk -F'\t' '
            $3=="primary"                     {tot=$1}
            $3=="primary mapped"              {m=$1}
            $3=="paired in sequencing"        {pr=$1}
            $3=="properly paired"             {pp=$1}
            $3=="with itself and mate mapped" {bm=$1}
            END {printf "%d\t%.3f\t%.3f\t%.3f", tot, (tot?100*m/tot:0), (pr?100*bm/pr:0), (pr?100*pp/pr:0)}' "${P}.flagstat.tsv")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ID" "$SP" "$D_IN" "$FRAC" "$SEEDFRAC" "$STATUS" "$D_OUT" "$DEPTH" "$QCF" "$MAP" >> "${TABLE}.tmp"
        echo "${OUT}/${ID}.cram" >> "${NEWLIST}.tmp"
    done < "$INLIST"
    mv -f "${TABLE}.tmp" "$TABLE"; mv -f "${NEWLIST}.tmp" "$NEWLIST"

    echo "Wrote $TABLE ($(($(wc -l < "$TABLE")-1)) samples) and $NEWLIST"
    awk -F'\t' 'NR>1 {s[$6]++; f[$11]++} END {
        printf "  status: "; for (k in s) printf "%s=%d  ", k, s[k]; print "";
        printf "  depth vs %s target: ", "'"$TARGET_DEPTH"'"; for (k in f) printf "%s=%d  ", k, f[k]; print ""}' "$TABLE"
    awk -F'\t' 'NR>1 && $11!="OK" {printf "    %s %s  %.2fx (%s%%)  %s\n", $11, $1, $7, $10, $6}' "$TABLE" | head -40
    if (( ${#missing[@]} )); then
        echo "WARNING: ${#missing[@]} sample(s) not finished, left out of the table AND the new list: ${missing[*]:0:30}" >&2
        echo "         rerun 'bash $0 submit' -- finished samples are skipped." >&2
    fi

    module --force purge 2>/dev/null || true
    ml r 2>/dev/null || ml R 2>/dev/null || true
    if command -v Rscript >/dev/null; then
        Rscript "$PLOT_R" "$TABLE" "${QC}/depth_mapping_by_species"
    else
        echo "Rscript not found; plot locally: Rscript plot_depth_mapping.R depth_mapping_summary.tsv depth_mapping_by_species" >&2
    fi
    exit 0
fi

echo "Usage: bash $0 submit | summary" >&2
exit 1
