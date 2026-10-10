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
#   bash nexus_pca_admix.sh check            # sample count, chunks, ref cache; no jobs
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
#
# CRAM REFERENCE: htslib looks reference sequences up by MD5 and, failing
# that, writes to ~/.cache/hts-ref and downloads from EBI; with 50 tasks x
# 506 CRAMs that races, fails, and silently drops slices. This script builds
# a COMPLETE local MD5 cache from REF_FASTA once (ref/hts-cache-full) and
# points REF_PATH/REF_CACHE at it, so no task touches the network or $HOME.
# Each chunk's ANGSD log is checked; any decode error fails the task and
# its output is not kept.
# =============================================================================
#SBATCH --job-name=nexus_pca
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
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
PCANGSD_RES=(32 200G 3-00:00:00)   # must fit one node: check  sinfo -p cpu -o "%c %m %l"
MAX_PARALLEL=50
PCA_THIN_BP=auto   # auto = thin only if the SNP set will not fit PCANGSD_RES memory;
                   # 0 = never thin; or a distance in bp (e.g. 5000)
HTS_CACHE="${PROJECT_DIR}/ref/hts-cache-full"   # complete MD5 cache built from REF_FASTA
PLOT_R="${SLURM_SUBMIT_DIR:-$PWD}/plot_pca_admix.R"
# =============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUT"/{chunks,beagle,pcangsd,plots} logs
cram_id() { local b; b=$(basename "$1"); echo "${b%%.*}"; }

xalt_fix() {
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
}

# Complete local MD5 reference cache: every sequence in REF_FASTA stored as
# <cache>/<md5[0:2]>/<md5[2:4]>/<md5[4:]> (uppercase, no newlines), the layout
# htslib expects. Built atomically, once; DONE marker lists the sequence count.
build_ref_cache() {
    [[ -s "${HTS_CACHE}/DONE" ]] && return 0
    echo "  building MD5 reference cache in $HTS_CACHE (one time, a few minutes) ..."
    rm -rf "${HTS_CACHE}.tmp"; mkdir -p "${HTS_CACHE}.tmp"
    python3 - "$REF_FASTA" "${HTS_CACHE}.tmp" <<'PY'
import sys, os, hashlib
fa, root = sys.argv[1], sys.argv[2]
n = 0
def put(name, chunks):
    global n
    if name is None: return
    seq = b"".join(chunks).upper()
    m = hashlib.md5(seq).hexdigest()
    d = os.path.join(root, m[:2], m[2:4]); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, m[4:]), "wb") as f: f.write(seq)
    n += 1
name, chunks = None, []
with open(fa, "rb") as f:
    for line in f:
        if line.startswith(b">"):
            put(name, chunks); name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
put(name, chunks)
print(n)
PY
    [[ $? -eq 0 ]] || { echo "ERROR: cache build failed" >&2; exit 1; }
    local nseq; nseq=$(find "${HTS_CACHE}.tmp" -type f | wc -l)
    echo "$nseq" > "${HTS_CACHE}.tmp/DONE"
    rm -rf "$HTS_CACHE"; mv "${HTS_CACHE}.tmp" "$HTS_CACHE"
    echo "  cache: $(cat "${HTS_CACHE}/DONE") sequences"
}

# Confirm every M5 in a CRAM header is in the cache (needs samtools; header
# reading does not need the reference). Run on the first CRAM in the list.
check_ref_cache() {
    local c m miss=0 n=0
    c=$(head -n 1 "${OUT}/bamlist.txt")
    module --force purge 2>/dev/null || true
    ml biocontainers samtools 2>/dev/null || true; xalt_fix
    command -v samtools >/dev/null || { echo "  (samtools not found; skipping M5 check)"; return 0; }
    while read -r m; do
        n=$((n+1)); [[ -s "${HTS_CACHE}/${m:0:2}/${m:2:2}/${m:4}" ]] || miss=$((miss+1))
    done < <(samtools view -H "$c" | grep '^@SQ' | grep -o 'M5:[0-9a-f]*' | cut -d: -f2)
    if (( n == 0 )); then echo "  WARNING: no M5 tags in $c header"; return 0; fi
    (( miss == 0 )) || { echo "ERROR: $miss of $n CRAM reference MD5s not in $HTS_CACHE -- wrong REF_FASTA?" >&2; exit 1; }
    echo "  cache check: all $n reference MD5s in $(basename "$c") found"
}

use_ref_cache() {   # htslib reads ONLY the local cache: no $HOME cache, no EBI
    export REF_PATH="${HTS_CACHE}/%2s/%2s/%s" REF_CACHE="${HTS_CACHE}/%2s/%2s/%s"
    export SINGULARITYENV_REF_PATH="$REF_PATH" SINGULARITYENV_REF_CACHE="$REF_CACHE"
    export APPTAINERENV_REF_PATH="$REF_PATH" APPTAINERENV_REF_CACHE="$REF_CACHE"
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
    build_ref_cache
    check_ref_cache
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
    ml biocontainers angsd/0.940; xalt_fix; use_ref_cache
    [[ -s "${HTS_CACHE}/DONE" ]] || { echo "ERROR: $HTS_CACHE missing; run: bash $0 check" >&2; exit 1; }
    rm -f "${B}".tmp.*
    echo ">>> chunk $CH: $(wc -l < "${OUT}/chunks/chunk_${CH}.rf") contigs, N=$N, minInd=$MININD  $(date)"
    angsd -bam "${OUT}/bamlist.txt" -ref "$REF_FASTA" -rf "${OUT}/chunks/chunk_${CH}.rf" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf "$MINMAF" -SNP_pval "$SNP_PVAL" \
        -skipTriallelic 1 -minMapQ "$MINMAPQ" -minQ "$MINQ" \
        -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -minInd "$MININD" -P "$THREADS" -out "${B}.tmp" 2> "${B}.angsd.log" || true
    tail -n 5 "${B}.angsd.log" >&2
    if grep -qE 'cram_decode_slice|cram_next_slice|cram_get_ref|Unable to fetch reference|Failed to populate reference|find_file_url' "${B}.angsd.log"; then
        echo "ERROR: CRAM reference/decode errors in ${B}.angsd.log -- output discarded" >&2
        grep -m 5 -E 'cram_|reference' "${B}.angsd.log" >&2; rm -f "${B}".tmp.*; exit 1
    fi
    grep -q 'ALL done' "${B}.angsd.log" || { echo "ERROR: ANGSD did not finish (see ${B}.angsd.log)" >&2; rm -f "${B}".tmp.*; exit 1; }
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
        # header from the first chunk, then every chunk without its header.
        # (header read with awk 'NR==1{print;exit}' || true: "zcat | head" ends
        #  with SIGPIPE, which set -o pipefail would treat as a fatal error)
        GZ="gzip"; command -v pigz >/dev/null && GZ="pigz -p ${THREADS}"
        { zcat "${PARTS[0]}" 2>/dev/null | awk 'NR==1{print; exit}' || true
          for f in "${PARTS[@]}"; do zcat "$f" | tail -n +2; done; } | $GZ > "${BEAGLE}.tmp"
        mv -f "${BEAGLE}.tmp" "$BEAGLE"
        echo ">>> concatenated  $(date)"
    fi
    NSNP=$(( $(zcat "$BEAGLE" | wc -l) - 1 )); echo ">>> $NSNP autosomal SNPs"
    # beagle columns must match samples.txt (3 GL columns per sample)
    NCOL=$( (zcat "$BEAGLE" 2>/dev/null | awk 'NR==1{print (NF-3)/3; exit}') || true )
    [[ "$NCOL" == "$(wc -l < "${OUT}/samples.txt")" ]] \
        || { echo "ERROR: beagle has $NCOL samples, samples.txt has $(wc -l < "${OUT}/samples.txt")" >&2; exit 1; }

    # PCAngsd holds ~12 bytes per SNP per sample (float32 likelihoods + working
    # matrices). If that exceeds ~70% of the job memory, thin by distance:
    # d = autosome length / SNPs that fit, which guarantees the kept set fits.
    NS=$(wc -l < "${OUT}/samples.txt")
    MEMGB=$(echo "${PCANGSD_RES[1]}" | tr -dc '0-9')
    FIT=$(awk -v g="$MEMGB" -v n="$NS" 'BEGIN{printf "%d", 0.70*g*1e9/(12*n)}')
    echo ">>> ~$(awk -v m="$NSNP" -v n="$NS" 'BEGIN{printf "%.0f", 12*m*n/1e9}') GB needed for all SNPs; ${MEMGB} GB job fits ~${FIT} SNPs"
    THIN="$PCA_THIN_BP"
    if [[ "$THIN" == auto ]]; then
        if (( NSNP > FIT )); then
            THIN=$(awk -v L="$(awk '{s+=$2} END{print s}' "${OUT}/autosomes.len")" -v f="$FIT" 'BEGIN{d=L/f; print (d==int(d))?d:int(d)+1}')
        else THIN=0; fi
    fi
    if (( THIN > 0 )); then
        TB="${OUT}/autosomes.thin${THIN}bp.beagle.gz"
        if [[ ! -s "$TB" ]]; then
            echo ">>> thinning to >= ${THIN} bp between SNPs  $(date)"
            # marker = <contig>_<pos>; contig names contain '_', so pos is after the LAST one
            zcat "$BEAGLE" | awk -v d="$THIN" 'NR==1 {print; next}
                { m=$1; p=m; sub(/.*_/,"",p); c=substr(m,1,length(m)-length(p)-1)
                  if (c!=lc || p-lp>=d) {print; lc=c; lp=p} }' | gzip > "${TB}.tmp"
            mv -f "${TB}.tmp" "$TB"
        fi
        BEAGLE="$TB"; NSNP=$(( $(zcat "$BEAGLE" | wc -l) - 1 ))
        echo ">>> PCAngsd input: $NSNP SNPs after thinning (${THIN} bp)" | tee "${OUT}/pcangsd/thinning.txt"
    fi

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
