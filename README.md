# extended_nexus

Range-wide *Tympanuchus* resequencing panel: the 433-bird panel of
[Black et al. 2024, *PNAS Nexus*](https://academic.oup.com/pnasnexus/article/3/8/pgae298/7720645)
plus 73 newly sequenced birds (15 LEPC, 27 GRPC, 31 STGR; five more LEPC
libraries replaced lower-depth data for birds already in the panel), for
**506 birds**: 426 Lesser Prairie-Chicken (LEPC), 49 Greater Prairie-Chicken
(GRPC), 29 Sharp-tailed Grouse (STGR) and 2 putative STGR x GRPC hybrids.
Every alignment is harmonized to the original panel's mean depth (4.66x).

This repository is the code for **Objective 2** of the USFWS report
*Grouse genomics, October 2026* ("Expanded range-wide resequencing
panel"). Sequence data: NCBI BioProjects PRJNA1513026 (embargoed) and
PRJNA986511.

- Reference: LEPC `pur_lepc_1.0` (GCF_026119805.1) for all three species.
  Z-linked scaffolds `NW_026294758.1` and `NW_026294813.1` are excluded from
  every analysis; f_ROH denominator = autosomal length (~920 Mb).
  (Black et al. 2024 included Z and used the full genome length, so their
  f_ROH values are not directly comparable.)
- Cluster: Purdue RCAC (Gautschi; depth step tested on Negishi), SLURM
  account `fnrdewoody`. Paths are under `${CLUSTER_SCRATCH}/GROUSE/nexus`.

## Pipeline (report section -> scripts)

| Step | Report (Objective 2) | Script | Key settings |
|---|---|---|---|
| 1. Mapping | Panel and depth harmonization | `nf-core/` (nf-core/sarek v3.8.1, `step: mapping`) | fastp (min length 75), BWA-MEM, duplicate marking |
| 2. Depth harmonization + QC | same; Fig 7 | `processing/nexus_downsample_depth_mapping.sh`, `R/plot_depth_mapping.R` | subsample to 4.66x genome-wide mean (fixed seed); samples below target kept at native depth; samtools coverage/flagstat |
| 3. Heterozygosity | Diversity, inbreeding...; Fig 8 | `analysis/nexus_heterozygosity.sh`, `R/plot_heterozygosity.R` | ANGSD 0.940 `-dosaf 1 -GL 1 -minQ 30 -minMapQ 30 -setMinDepth 3 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1`, autosomes; folded SFS. Shapiro-Wilk, Kruskal-Wallis, BH-adjusted pairwise Wilcoxon, rank-biserial r |
| 4. ROH / f_ROH | same; Fig 9 | `analysis/run_bcftools_roh.sh`, `analysis/roh_parse_autosomal.sh` | ANGSD `-GL 1 -doBcf 1 -doPost 1 -minQ 30 -SNP_pval 1e-6` -> `bcftools roh`; quality >= 30; Z excluded; classes 100 kb-1 Mb and > 1 Mb (ROH < 100 kb not counted) |
| 5. f_ROH statistics | same (species, allopatric vs sympatric, depth check) | `R/fROH_statistics.R`, `R/plot_roh.R` | Kruskal-Wallis + BH pairwise Wilcoxon (hybrids excluded); allopatric vs sympatric within species, BH across 9 tests; Spearman f_ROH vs depth; LEPC `lm(rank(f_ROH) ~ GROUP + depth)` |
| 6. PCA | same; Fig 10 | `analysis/nexus_pca_admix.sh`, `R/plot_pca.R` (`R/plot_pca_admix.R` for the automatic plots) | ANGSD beagle GLs on all 506 birds (autosomal chunks, local MD5 reference cache) -> PCAngsd |
| 7. Admixture | same; Fig 11 | `analysis/nexus_ngsadmix.sh`, `R/plot_ngsadmix.R` | NGSadmix on the same beagle file thinned to 1 SNP / 10 kb; K = 1-10, 10 runs per K (`-minMaf 0.05 -maxiter 50000 -tol 1e-9`); log-likelihood convergence + Evanno delta K; K = 2-4 plotted. NGSadmix is built from `ngsadmix32.cpp` (see script header) |
| 8. F_ST | same; Fig 12 | `analysis/angsd_fst_species.sh`, `R/plot_fst_sliding_window.R` | SNPs discovered jointly in the 504 non-hybrid birds; folded 2D-SFS prior; Hudson (`-whichFst 1`) and Reynolds (`-whichFst 0`); 100-kb windows, 20-kb step; windows with < 500 SNPs dropped in the plot |

Kept for reference, not used in the report: the PCAngsd admixture output of
`nexus_pca_admix.sh`, and `analysis/run_ROHan.sh` (ROHan cross-check).

## Layout

```
nf-core/       sarek config, params, sample sheet
processing/    depth harmonization + depth/mapping QC
analysis/      heterozygosity, ROH, PCA/admixture, F_ST (SLURM)
R/             statistics and report figures (run locally)
```

## Running

Submit from the folder holding each script. Multi-stage scripts chain
their own SLURM dependencies:

```bash
bash processing/nexus_downsample_depth_mapping.sh submit
bash analysis/nexus_heterozygosity.sh submit
sbatch analysis/run_bcftools_roh.sh && sbatch analysis/roh_parse_autosomal.sh
bash analysis/nexus_pca_admix.sh submit
bash analysis/nexus_ngsadmix.sh submit        # after nexus_pca_admix.sh beagle stage
bash analysis/angsd_fst_species.sh submit
```

Notes for RCAC: `xalt` injects `LD_PRELOAD` into biocontainers (scripts
blank it via `SINGULARITYENV_/APPTAINERENV_LD_PRELOAD=""`); CRAM decoding
uses a local MD5 reference cache (`REF_PATH`/`REF_CACHE`) built from the
reference, so no EBI lookups are needed.

The R scripts expect `sub` / metadata objects built from
`nexus_metadata.xlsx` (columns `ID`, `GROUP`, `common`, f_ROH classes,
depth). Two hybrid rows in that workbook are column-shifted; scripts coerce
numeric columns with `as.numeric()`.
