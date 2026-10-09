#!/usr/bin/env Rscript
# Plot mean depth and paired mapping rate by species from
# depth_mapping_summary.tsv (written by nexus_depth_mapping.sh summary).
#
#   Rscript plot_depth_mapping.R depth_mapping_summary.tsv depth_mapping_by_species
#
# Writes <prefix>.pdf, <prefix>.png and <prefix>_by_species.tsv (n, mean, SD,
# median, range per species). Needs ggplot2; patchwork

library(ggplot2)
library(patchwork)
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "depth_mapping_summary.tsv"
prefix <- if (length(args) >= 2) args[2] else "depth_mapping_by_species"

# Read robustly: a stray extra tab in a row (e.g., an empty species field
# from a double tab in the popmap) would otherwise shift that row's columns
# or wrap it onto a new row. Empty fields after the sample ID are dropped
# when a row has more fields than the header.
lines <- readLines(infile, warn = FALSE)
lines <- sub("\r$", "", lines[nzchar(lines)])
hdr   <- strsplit(lines[1], "\t", fixed = FALSE)[[1]]
hdr   <- trimws(hdr[nzchar(trimws(hdr))])     # a header name is never empty: drop stray tabs
rows  <- lapply(strsplit(lines[-1], "\t"), function(f) {
  if (length(f) > length(hdr)) f <- c(f[1], f[-1][nzchar(f[-1])])
  length(f) <- length(hdr); f
})
d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
# empty header names (a double tab in the header line): drop the column if it
# holds no data, otherwise give it a placeholder name
empty <- which(!nzchar(trimws(hdr)))
for (i in empty) hdr[i] <- paste0("unnamed_", i)
names(d) <- hdr
drop <- names(d)[startsWith(names(d), "unnamed_") &
                 vapply(d, function(z) all(is.na(z) | !nzchar(trimws(z))), logical(1))]
if (length(drop)) d <- d[, setdiff(names(d), drop), drop = FALSE]
hdr <- names(d)
for (v in setdiff(names(d), c("sample", "sample_id", "species", "status", "depth_flag", "seed_fraction")))
  suppressWarnings(d[[v]] <- as.numeric(d[[v]]))
cat("Columns:", paste(hdr, collapse = ", "), "\n")
cat("Rows:", nrow(d), "\n")
d$species <- trimws(d$species)
d$species[grepl("/", d$species)] <- "STGR x GRPC"   # putative hybrids
lev <- c("LEPC", "GRPC", "STGR", "STGR x GRPC", "unassigned")
d$species <- factor(d$species, levels = c(intersect(lev, unique(d$species)),
                                          setdiff(unique(d$species), lev)))
cols <- c(LEPC = "#A52A2A", GRPC = "#DAA520", STGR = "#000000",
          "STGR x GRPC" = "#BEBEBE", unassigned = "#888888")
n_lab <- function(x) paste0(levels(x), "\n(n=", table(x), ")")

panel <- function(y, ylab, hline = NULL) {
  if (!y %in% names(d)) { message("Column ", y, " not in table -- panel skipped"); return(NULL) }
  p <- ggplot(d, aes(species, .data[[y]], colour = species)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey30", fill = NA) +
    geom_jitter(width = 0.18, height = 0, size = 1.3, alpha = 0.65) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_discrete(labels = n_lab(d$species)) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11)
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2, colour = "grey50")
  p
}

target <- 4.66
DPI    <- 800
# Final (post-subsampling) data only: mean depth and properly paired reads
depth_col <- if ("mean_depth_genome" %in% names(d)) "mean_depth_genome" else "mean_depth_autosomal"
pA <- panel(depth_col, "Mean depth (X)", hline = target)
pB <- panel("pct_properly_paired", "Properly paired reads (%)")
pl <- Filter(Negate(is.null), list(pA, pB))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- wrap_plots(pl, ncol = 2) + plot_annotation(tag_levels = "A")
  ggsave(paste0(prefix, ".pdf"), fig, width = 9, height = 3.8)
  ggsave(paste0(prefix, ".png"), fig, width = 9, height = 3.8, dpi = DPI)
} else {
  pdf(paste0(prefix, ".pdf"), width = 4.5, height = 3.8)
  for (pp in pl) print(pp); invisible(dev.off())
  ggsave(paste0(prefix, "_depth.png"), pA, width = 4.5, height = 3.8, dpi = DPI)
  ggsave(paste0(prefix, "_properly_paired.png"), pB, width = 4.5, height = 3.8, dpi = DPI)
}

# per-species summary table
vars <- c("mean_depth_genome", "mean_depth_autosomal", "pct_both_mates_mapped", "pct_properly_paired", "pct_mapped")
vars <- intersect(vars, names(d))
s <- do.call(rbind, lapply(split(d, d$species, drop = TRUE), function(g) {
  data.frame(species = as.character(g$species[1]), n = nrow(g),
             do.call(cbind, lapply(vars, function(v) {
               x <- g[[v]]
               setNames(data.frame(round(mean(x), 3), round(sd(x), 3), round(median(x), 3),
                                   round(min(x), 3), round(max(x), 3)),
                        paste0(v, c("_mean", "_sd", "_median", "_min", "_max")))
             })))
}))
write.table(s, paste0(prefix, "_by_species.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(s[, c("species", "n", "mean_depth_autosomal_mean", "mean_depth_autosomal_sd",
            "pct_both_mates_mapped_mean", "pct_properly_paired_mean")], row.names = FALSE)
cat("Wrote", paste0(prefix, c(".pdf", ".png", "_by_species.tsv")), sep = "\n  ")
