#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: ROHan only (nexus project, all samples)
#
# Standalone Step 7 -- factored out of 08b_roh_resume.sh so ROHan no longer
# has to run after (or depend on) the ANGSD/bcftools ROH steps (3-6) in that
# script. Only requires 06_downsample_and_finalize.sh (final_cramlist.txt)
# to have completed.
#
# NOTE: 08b_roh_resume.sh's own Step 7 (--rohmu 2e-5, no --size, i.e. the
# 1Mb default window) is now superseded by this script and should be
# skipped/removed there -- otherwise ROHan ends up run twice per sample
# under two different, inconsistent parameter sets.
#
# Uses the --rohmu and --size values, and the post-hoc length-class
# parsing (100kb-1Mb / >1Mb, comparable to bcftools roh's own classes),
# established for the old_vs_new ROHan comparison (run_rohan.sh):
#   - ROHMU=1e-5 (not nexus's earlier 2e-5)
#   - ROHAN_WINDOW_SIZE=50000 (not ROHan's 1Mb default -- see CONFIG note)
#   - Step 2 parsing below is copied over unchanged, only with nexus's
#     REF_FASTA/ROHAN_OUT_DIR substituted in.
#
# Requires install_rohan.sh to have been run on the login node first.
#
# USAGE:
#   sbatch 08c_rohan.sh
# =============================================================================
#SBATCH --job-name=nexus_rohan
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/nexus"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist_4.66x.txt"

# ROHan binary + its GSL build are shared from the old_vs_new project
# rather than rebuilt here -- built once, reused across projects (see
# install_rohan.sh).
ROHAN_BIN="/scratch/gautschi/blackan/GROUSE/old_vs_new/tools/ROHan/bin/rohan"
GSL_PREFIX="/scratch/gautschi/blackan/GROUSE/old_vs_new/tools/gsl"

ROHAN_OUT_DIR="${PROJECT_DIR}/results/rohan"
ROHAN_THREADS=16

# ROHMU: ROHan's expected within-ROH heterozygosity rate parameter. Set to
# the value used for the old_vs_new comparison, not nexus's earlier 2e-5.
ROHMU=1e-5

# ROHAN_WINDOW_SIZE: ROHan's HMM operates on the genome tiled into fixed
# windows (--size). At the 1Mb default it cannot resolve anything in the
# 100kb-1Mb class at all -- confirmed empirically in the old_vs_new
# comparison: ROH segment lengths came back as exact multiples of
# 1,000,000 across every sample, with near-zero total ROH, while bcftools
# roh found substantial 100kb-1Mb signal in the same individuals. 50kb
# gives several windows of resolution within that bin (enough to
# distinguish it from the >1Mb class) without pushing per-window
# segregating-site counts so low that local theta estimates become
# unreliably noisy.
ROHAN_WINDOW_SIZE=50000

mkdir -p "$ROHAN_OUT_DIR" logs

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found. Run 06_downsample_and_finalize.sh first." >&2
    exit 1
fi
if [ ! -x "$ROHAN_BIN" ]; then
    echo "ERROR: ROHan binary not found at $ROHAN_BIN" >&2
    echo "  Run install_rohan.sh on the login node first." >&2
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found." >&2
    exit 1
fi

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")
echo ">>> 08c_rohan.sh (ROHan only, rohmu=${ROHMU}, window=${ROHAN_WINDOW_SIZE})"
echo ">>> N samples : ${N_SAMPLES}"
echo ">>> Start time: $(date)"

module --force purge
module load gcc/14.1.0
module load biocontainers
module load samtools

# rohan was linked against a custom-built GSL (no 'gsl' module exists on
# Gautschi), so it needs this to find libgsl.so at runtime, not just at
# build time -- see install_rohan.sh.
export LD_LIBRARY_PATH="$GSL_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# =============================================================================
# STEP 1: Run ROHan per sample
# =============================================================================
echo ""
echo ">>> Step 1: ROHan per-sample analysis"

while IFS= read -r CRAM; do
    SAMPLE=$(basename "$CRAM" | sed -E 's/\.cram$//')

    # Window size baked into the output prefix -- keeps this run distinct
    # from nexus's earlier 1Mb-window/2e-5 ROHan output for the same
    # samples (08b_roh_resume.sh), which used unsuffixed filenames.
    OUT_PREFIX="${ROHAN_OUT_DIR}/${SAMPLE}.win${ROHAN_WINDOW_SIZE}"
    if [[ -f "${OUT_PREFIX}.hEst.gz" ]]; then
        echo "  ${SAMPLE}: ${OUT_PREFIX}.hEst.gz already exists -- skipping."
        continue
    fi
    echo "=== Processing $SAMPLE ==="

    # ROHan's most reliably supported input format is BAM. Rather than rely
    # on ROHan's own CRAM/reference handling (undocumented in what we could
    # verify), convert to a temporary indexed BAM first -- slower, but
    # removes any ambiguity about reference resolution.
    TMP_BAM="${ROHAN_OUT_DIR}/${SAMPLE}.tmp.bam"
    samtools view -@ "$ROHAN_THREADS" -b -T "$REF_FASTA" -o "$TMP_BAM" "$CRAM"
    samtools index "$TMP_BAM"

    "$ROHAN_BIN" \
        -t "$ROHAN_THREADS" \
        --rohmu "$ROHMU" \
        --size "$ROHAN_WINDOW_SIZE" \
        -o "$OUT_PREFIX" \
        "$REF_FASTA" "$TMP_BAM"

    rm -f "$TMP_BAM" "${TMP_BAM}.bai"
    echo "  Done: ${OUT_PREFIX}.*"
done < "$FINAL_CRAMLIST"

echo ""
echo ">>> Step 1 complete. Per-sample outputs in ${ROHAN_OUT_DIR}/"
echo "Each sample produces (filenames per ROHan's own convention):"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.hEst.gz          -- genome-wide heterozygosity estimate + CI, outside ROH"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.mid.hmmp.gz      -- per-window HMM posterior probabilities"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz   -- called ROH segments (chrom, begin, end, length)"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.summary.txt      -- genome-wide summary stats"

# =============================================================================
# STEP 2: Parse per-segment ROH output into bcftools-roh-comparable bins
#
# hmmrohl format (confirmed from actual output):
#   #ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
# One row per called ROH segment; a sample with zero ROH has no data rows
# (header only).
# =============================================================================
echo ""
echo ">>> Step 2: Parsing ROH segments into 100kb-1Mb / >1Mb bins"

# Genome length denominator for fROH, summed from the reference .fai --
# matches the convention rohparser.py uses on the ANGSD/bcftools roh side
# (08b_roh_resume.sh), so the two methods' fROH values are computed the
# same way.
GENOME_LEN=$(awk '{sum+=$2} END{print sum}' "${REF_FASTA}.fai")
echo "  Genome length (from .fai): ${GENOME_LEN} bp"

OUT_TSV="${ROHAN_OUT_DIR}/rohan_froh_summary.win${ROHAN_WINDOW_SIZE}.tsv"
echo -e "sample\tfROH_100kb-1Mb_rohan\tfROH_1Mb_rohan\tfROH_total_rohan\tn_segments_short\tn_segments_long" > "$OUT_TSV"

N_PARSED=0
for HMMROHL in "$ROHAN_OUT_DIR"/*.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz; do
    [[ -e "$HMMROHL" ]] || { echo "ERROR: no *.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz files found in ${ROHAN_OUT_DIR}" >&2; exit 1; }

    SAMPLE=$(basename "$HMMROHL" | sed -E "s/\.win${ROHAN_WINDOW_SIZE}\.mid\.hmmrohl\.gz$//")

    zcat "$HMMROHL" 2>/dev/null | awk -v genome="$GENOME_LEN" -v sample="$SAMPLE" '
        NR > 1 {
            len = $5
            if (len >= 100000 && len <= 1000000) { short += len; n_short++ }
            else if (len > 1000000) { long += len; n_long++ }
        }
        END {
            froh_short = short / genome
            froh_long  = long / genome
            froh_total = (short + long) / genome
            printf "%s\t%.9f\t%.9f\t%.9f\t%d\t%d\n", sample, froh_short, froh_long, froh_total, n_short+0, n_long+0
        }' >> "$OUT_TSV"

    N_PARSED=$((N_PARSED + 1))
done

echo "  Parsed ${N_PARSED} samples"
echo "  Written: ${OUT_TSV}"
echo ""
column -t "$OUT_TSV" 2>/dev/null || cat "$OUT_TSV"

echo ""
echo ">>> All done."
echo ">>> End time: $(date)"
