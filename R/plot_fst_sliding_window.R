# =============================================================================
# Sliding-window FST (ANGSD realSFS fst stats2) for the three species pairs,
# nexus panel. One panel per pair, stacked in one column. Scaffolds are laid
# end to end (longest first) and unlabeled; point colors alternate between
# scaffolds using the two species' colors of each pair.
# Dashed line = genome-wide weighted FST for that pair.
# =============================================================================
library(ggplot2)

FST_DIR <- "/path/to/nexus/fst/species_LEPC_GRPC_STGR"   # <- edit (copy of the cluster OUT dir)
WIN  <- 100000           # window size used in realSFS fst stats2 (50000 or 100000)
STEP <- 20000            # step used (10000 for the 50 kb run)
EST  <- ""               # "" = Hudson files; ".reynolds" = Reynolds files
WIN_FILE_SUFFIX <- sprintf("%s.fst.windows_%d_%d.txt", EST, WIN, STEP)
MIN_SITES <- 10          # drop windows with fewer SNPs than this (noisy FST)

# Genome-wide weighted FST per pair, if the .fst.global.txt files are not in
# FST_DIR: paste the second number from `cat <A>_<B>.fst.global.txt` here.
# NA = read from file.
# Values from realSFS fst stats (weighted, Hudson / -whichFst 1); set to NA to
# read the .fst.global.txt files from FST_DIR instead.
GLOBAL_FST <- c(LEPC_GRPC = 0.071541, LEPC_STGR = 0.418855, GRPC_STGR = 0.414708)

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
# realSFS fst stats prints "unweighted weighted" (two numbers). Older builds
# print "FST.Unweight[nObs:N]:x Fst.Weight:y" instead; both are handled by
# taking the last number on the last non-empty line.
read_global <- function(A, B) {
  key <- paste0(A, "_", B)
  if (!is.na(GLOBAL_FST[key])) return(unname(GLOBAL_FST[key]))
  f <- file.path(FST_DIR, paste0(key, EST, ".fst.global.txt"))
  if (!file.exists(f)) { warning("No global FST file: ", f); return(NA_real_) }
  l <- readLines(f, warn = FALSE); l <- l[nzchar(trimws(l))]
  if (!length(l)) { warning("Empty global FST file: ", f); return(NA_real_) }
  nums <- regmatches(tail(l, 1), gregexpr("-?[0-9.]+(e-?[0-9]+)?", tail(l, 1)))[[1]]
  as.numeric(tail(nums, 1))
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
# Fallback: SNP-weighted mean of window FST, clearly labelled as approximate
approx_fst <- tapply(win$fst * win$Nsites, win$pair, sum) / tapply(win$Nsites, win$pair, sum)
glob$approx <- is.na(glob$fst)
glob$fst[glob$approx] <- approx_fst[as.character(glob$pair[glob$approx])]
glob$label <- ifelse(glob$approx,
                     sprintf("approx.~weighted~F[ST] == %.3f", glob$fst),
                     sprintf("weighted~F[ST] == %.3f", glob$fst))
print(glob)

# ---- plot -------------------------------------------------------------------
p <- ggplot(win, aes(x = pos, y = fst)) +
  geom_point(aes(color = col), size = 0.4, alpha = 0.7) +
  scale_color_identity() +
  geom_hline(data = glob, aes(yintercept = fst),
             linetype = "dashed", linewidth = 0.7, color = "red3") +
  geom_label(data = glob, aes(label = label), parse = TRUE, label.size = 0, fill = "white",
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~ pair, ncol = 1, scales = "fixed") +   # shared 0-1 axis so LEPC-GRPC and STGR pairs compare
  scale_x_continuous(expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(min(0, min(win$fst)), 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Autosomal scaffolds (ordered by length)", y = bquote(F[ST]~"("*.(WIN/1000)~"kb windows)")) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(color = "black", fill = "grey90", linewidth = 1),
        strip.text = element_text(size = 13, face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 14))
print(p)
ggsave(sprintf("fst_windows_species_%dkb%s.png", WIN / 1000, EST), p, width = 12, height = 8, dpi = 300)
