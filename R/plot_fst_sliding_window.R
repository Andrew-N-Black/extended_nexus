# =============================================================================
# Sliding-window FST (ANGSD realSFS fst stats2) for the three species pairs,
# nexus panel. One panel per pair, stacked in one column. Scaffolds are laid
# end to end (longest first) and unlabeled; point colors alternate between
# scaffolds using the two species' colors of each pair.
# Dashed line = genome-wide weighted FST for that pair.
# =============================================================================
library(ggplot2)

FST_DIR <- "/path/to/nexus/fst/species_LEPC_GRPC_STGR"   # <- edit (copy of the cluster OUT dir)
WIN_FILE_SUFFIX <- ".fst.windows_50000_10000.txt"
MIN_SITES <- 10          # drop windows with fewer SNPs than this (noisy FST)

species_col <- c("Tympanuchus cupido"         = "goldenrod",
                 "Tympanuchus pallidicinctus" = "brown",
                 "Tympanuchus phasianellus"   = "black")
code2sp <- c(GRPC = "Tympanuchus cupido",
             LEPC = "Tympanuchus pallidicinctus",
             STGR = "Tympanuchus phasianellus")
pairs <- list(c("LEPC", "GRPC"), c("LEPC", "STGR"), c("GRPC", "STGR"))   # file prefixes A_B

short_name <- function(sp) sub("^Tympanuchus ", "T. ", sp)

# ---- read windows ----------------------------------------------------------
# stats2 writes a 4-name header (region chr midPos Nsites) above 5-column rows,
# so the header is skipped and columns are named here.
read_windows <- function(A, B) {
  f <- file.path(FST_DIR, paste0(A, "_", B, WIN_FILE_SUFFIX))
  if (!file.exists(f)) stop("Missing: ", f)
  d <- read.table(f, header = FALSE, skip = 1, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("region", "chr", "midPos", "Nsites", "fst"))
  d$A <- A; d$B <- B
  d
}
read_global <- function(A, B) {
  f <- file.path(FST_DIR, paste0(A, "_", B, ".fst.global.txt"))
  if (!file.exists(f)) return(NA_real_)
  as.numeric(strsplit(trimws(tail(readLines(f), 1)), "[[:space:]]+")[[1]][2])  # weighted
}

win <- do.call(rbind, lapply(pairs, function(p) read_windows(p[1], p[2])))
win <- win[win$Nsites >= MIN_SITES & is.finite(win$fst), ]

# ---- scaffold order and cumulative positions --------------------------------
# Use the autosome lengths written by the FST script if present; otherwise the
# largest window midpoint stands in for scaffold length.
len_file <- file.path(FST_DIR, "autosomes.len")
if (file.exists(len_file)) {
  scaf <- read.table(len_file, col.names = c("chr", "len"), stringsAsFactors = FALSE)
} else {
  scaf <- aggregate(midPos ~ chr, data = win, FUN = max); names(scaf)[2] <- "len"
}
scaf <- scaf[scaf$chr %in% win$chr, ]
scaf <- scaf[order(-scaf$len), ]
scaf$offset <- c(0, cumsum(as.numeric(scaf$len))[-nrow(scaf)])
scaf$parity <- seq_len(nrow(scaf)) %% 2          # same alternation in every panel

win <- merge(win, scaf[, c("chr", "offset", "parity")], by = "chr")
win$pos <- win$midPos + win$offset

# ---- per-pair labels and alternating colors ---------------------------------
win$spA <- code2sp[win$A]; win$spB <- code2sp[win$B]
win$pair <- paste0(short_name(win$spA), " vs. ", short_name(win$spB))
win$col  <- ifelse(win$parity == 1, species_col[win$spA], species_col[win$spB])
pair_levels <- sapply(pairs, function(p) paste0(short_name(code2sp[p[1]]), " vs. ", short_name(code2sp[p[2]])))
win$pair <- factor(win$pair, levels = pair_levels)

glob <- data.frame(pair = factor(pair_levels, levels = pair_levels),
                   fst  = sapply(pairs, function(p) read_global(p[1], p[2])))
glob$label <- sprintf("weighted F[ST] == %.3f", glob$fst)

# ---- plot -------------------------------------------------------------------
p <- ggplot(win, aes(x = pos, y = fst)) +
  geom_point(aes(color = col), size = 0.4, alpha = 0.7) +
  scale_color_identity() +
  geom_hline(data = glob[is.finite(glob$fst), ], aes(yintercept = fst),
             linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_text(data = glob[is.finite(glob$fst), ], aes(label = label), parse = TRUE,
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~ pair, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0.005, 0)) +
  labs(x = "Autosomal scaffolds (ordered by length)", y = expression(F[ST]~"(50 kb windows)")) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(color = "black", fill = "grey90", linewidth = 1),
        strip.text = element_text(size = 13, face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 14))
print(p)
ggsave("fst_windows_species.png", p, width = 12, height = 8, dpi = 300)
