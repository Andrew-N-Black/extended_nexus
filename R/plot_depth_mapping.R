#!/usr/bin/env Rscript
# Plot mean depth and paired mapping rate by species from
# depth_mapping_summary.tsv (written by nexus_depth_mapping.sh summary).
#
#   Rscript plot_depth_mapping.R depth_mapping_summary.tsv depth_mapping_by_species
#
# Writes <prefix>.pdf, <prefix>.png and <prefix>_by_species.tsv (n, mean, SD,
# median, range per species). Needs ggplot2; patchwork is used if installed.

suppressPackageStartupMessages(library(ggplot2))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "depth_mapping_summary.tsv"
prefix <- if (length(args) >= 2) args[2] else "depth_mapping_by_species"

d <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)
d$species <- trimws(d$species)
d$species[grepl("/", d$species)] <- "STGR x GRPC"   # putative hybrids
lev <- c("LEPC", "GRPC", "STGR", "STGR x GRPC", "unassigned")
d$species <- factor(d$species, levels = c(intersect(lev, unique(d$species)),
                                          setdiff(unique(d$species), lev)))
cols <- c(LEPC = "#2a78b5", GRPC = "#c4762b", STGR = "#3a9a5b",
          "STGR x GRPC" = "#8a6bbe", unassigned = "#888888")
n_lab <- function(x) paste0(levels(x), "\n(n=", table(x), ")")

panel <- function(y, ylab, hline = NULL) {
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
if ("mean_depth_genome" %in% names(d) && "depth_before" %in% names(d)) {
  p1 <- panel("mean_depth_genome", "Mean depth after\nsubsampling (X)", hline = target)
} else {
  p1 <- panel("mean_depth_autosomal", "Mean autosomal depth (X)")
}
p0 <- if ("depth_before" %in% names(d)) panel("depth_before", "Mean depth before\nsubsampling (X)", hline = target) else NULL
p2 <- panel("pct_both_mates_mapped", "Paired mapping rate (%)\n(both mates mapped)")
p3 <- panel("pct_properly_paired", "Properly paired (%)")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- if (is.null(p0)) p1 + p2 + p3 else (p0 | p1) / (p2 | p3)
  fig <- fig + plot_annotation(tag_levels = "A")
  w <- if (is.null(p0)) 11 else 9; h <- if (is.null(p0)) 3.8 else 7
  ggsave(paste0(prefix, ".pdf"), fig, width = w, height = h)
  ggsave(paste0(prefix, ".png"), fig, width = w, height = h, dpi = 300)
} else {
  pdf(paste0(prefix, ".pdf"), width = 4.2, height = 3.8)
  if (!is.null(p0)) print(p0); print(p1); print(p2); print(p3); invisible(dev.off())
  for (k in list(list(p1, "depth"), list(p2, "paired_mapping"), list(p3, "properly_paired")))
    ggsave(paste0(prefix, "_", k[[2]], ".png"), k[[1]], width = 4.2, height = 3.8, dpi = 300)
}

# per-species summary table
vars <- c(intersect(c("depth_before", "mean_depth_genome"), names(d)), "mean_depth_autosomal", "pct_both_mates_mapped", "pct_properly_paired", "pct_mapped")
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
