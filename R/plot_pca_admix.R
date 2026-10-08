#!/usr/bin/env Rscript
# PCA and admixture plots from PCAngsd output (nexus_pca_admix.sh).
#
#   Rscript plot_pca_admix.R <OUT_DIR>
#
# Reads <OUT_DIR>/samples_species.tsv and <OUT_DIR>/pcangsd/nexus.cov plus
# nexus_K<k>.admix*.Q; writes <OUT_DIR>/plots/{pca,admixture}.{pdf,png} and
# pca_scores.tsv (PC1-PC10 per sample, for joining with other metadata).
# Needs ggplot2; patchwork (optional) combines panels.

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args) >= 1) args[1] else "."
pdir <- file.path(dir, "plots"); dir.create(pdir, showWarnings = FALSE)

meta <- read.delim(file.path(dir, "samples_species.tsv"), stringsAsFactors = FALSE)
meta$species <- trimws(meta$species)
meta$species[grepl("/", meta$species)] <- "STGR x GRPC"     # putative hybrids
lev  <- c("LEPC", "GRPC", "STGR", "STGR x GRPC", "unassigned")
meta$species <- factor(meta$species, levels = c(intersect(lev, unique(meta$species)),
                                                setdiff(unique(meta$species), lev)))
cols <- c(LEPC = "#2a78b5", GRPC = "#c4762b", STGR = "#3a9a5b",
          "STGR x GRPC" = "#8a6bbe", unassigned = "#888888")
shps <- c(LEPC = 16, GRPC = 17, STGR = 15, "STGR x GRPC" = 8, unassigned = 4)

# ---------------- PCA ----------------
C <- as.matrix(read.table(file.path(dir, "pcangsd", "nexus.cov")))
stopifnot(nrow(C) == nrow(meta))
e <- eigen(C, symmetric = TRUE)
pve <- 100 * e$values / sum(e$values)
npc <- min(10, ncol(e$vectors))
sc  <- data.frame(meta, e$vectors[, 1:npc])
names(sc)[-(1:2)] <- paste0("PC", 1:npc)
write.table(sc, file.path(pdir, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca_panel <- function(a, b) {
  ggplot(sc[order(sc$species == "STGR x GRPC"), ],          # hybrids drawn on top
         aes(.data[[paste0("PC", a)]], .data[[paste0("PC", b)]], colour = species, shape = species)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    scale_shape_manual(values = shps, name = NULL, drop = TRUE) +
    labs(x = sprintf("PC%d (%.1f%%)", a, pve[a]), y = sprintf("PC%d (%.1f%%)", b, pve[b])) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom")
}
p12 <- pca_panel(1, 2); p13 <- pca_panel(1, 3)
have_pw <- requireNamespace("patchwork", quietly = TRUE)
if (have_pw) {
  library(patchwork)
  fig <- (p12 | p13) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
         theme(legend.position = "bottom")
  ggsave(file.path(pdir, "pca.pdf"), fig, width = 10, height = 5)
  ggsave(file.path(pdir, "pca.png"), fig, width = 10, height = 5, dpi = 300)
} else {
  ggsave(file.path(pdir, "pca.pdf"), p12, width = 6, height = 5)
  ggsave(file.path(pdir, "pca.png"), p12, width = 6, height = 5, dpi = 300)
}
cat(sprintf("PCA: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%\n", pve[1], pve[2], pve[3]))

# ---------------- admixture ----------------
qfiles <- list.files(file.path(dir, "pcangsd"), pattern = "^nexus_K[0-9]+\\.admix.*\\.Q$", full.names = TRUE)
if (length(qfiles) == 0) { cat("No admixture .Q files found; skipping admixture plot\n"); quit(save = "no") }
Ks <- as.integer(sub("^nexus_K([0-9]+)\\..*$", "\\1", basename(qfiles)))
qfiles <- qfiles[order(Ks)]; Ks <- sort(Ks)

# one sample order for all K: by species, then by ancestry at the largest K
Qmax <- as.matrix(read.table(qfiles[length(qfiles)]))
stopifnot(nrow(Qmax) == nrow(meta))
dom  <- max.col(Qmax, ties.method = "first")
ord  <- order(meta$species, dom, -Qmax[cbind(seq_len(nrow(Qmax)), dom)])
pos  <- integer(nrow(meta)); pos[ord] <- seq_along(ord)

long <- do.call(rbind, lapply(seq_along(qfiles), function(i) {
  Q <- as.matrix(read.table(qfiles[i]))
  data.frame(x = rep(pos, ncol(Q)), K = paste0("K = ", Ks[i]),
             cluster = factor(rep(seq_len(ncol(Q)), each = nrow(Q))), q = as.vector(Q))
}))
long$K <- factor(long$K, levels = paste0("K = ", Ks))
anc_cols <- c("#2a78b5", "#c4762b", "#3a9a5b", "#8a6bbe", "#d1495b", "#edae49", "#00798c", "#66a182")

# species blocks along the x axis
sp_ord <- meta$species[ord]
brk <- cumsum(table(sp_ord)); brk <- brk[brk > 0]
mid <- brk - table(sp_ord)[names(brk)] / 2

p_adm <- ggplot(long, aes(x, q, fill = cluster)) +
  geom_col(width = 1) +
  geom_vline(xintercept = brk[-length(brk)] + 0.5, colour = "white", linewidth = 0.6) +
  facet_grid(K ~ .) +
  scale_fill_manual(values = anc_cols, guide = "none") +
  scale_x_continuous(breaks = as.numeric(mid), labels = names(mid), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.3, "lines"))
h <- 1.2 + 1.3 * length(Ks)
ggsave(file.path(pdir, "admixture.pdf"), p_adm, width = 11, height = h)
ggsave(file.path(pdir, "admixture.png"), p_adm, width = 11, height = h, dpi = 300)

# per-species mean ancestry at each K (quick check for the text)
for (i in seq_along(qfiles)) {
  Q <- as.matrix(read.table(qfiles[i]))
  m <- aggregate(Q, by = list(species = meta$species), FUN = mean)
  cat(sprintf("\nK = %d mean ancestry by species:\n", Ks[i])); print(m, digits = 3, row.names = FALSE)
}
cat("\nWrote", file.path(pdir, c("pca.pdf", "pca.png", "admixture.pdf", "admixture.png", "pca_scores.tsv")), sep = "\n  ")
