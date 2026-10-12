# =============================================================================
# plot_evaladmix.R -- evalAdmix residual-correlation heatmaps for the NGSadmix
# runs (USFWS report Objective 2; choosing K for Figure 11).
#
# Inputs (copy from ${OUT} on the cluster; see analysis/nexus_ngsadmix.sh):
#   evaladmix/nexus_k<K>.corres.txt   one per K, rows/cols in beagle sample order
#   nexus_metadata.xlsx               SAME row order as samples.txt (as for
#                                     plot_ngsadmix.R)
#
# Outputs (in <base>/evaladmix/plots/; one page per K in every PDF):
#   evaladmix_overview.pdf            all birds, species > sympatry > latitude
#                                     (the admixture barplot order)
#   evaladmix_by_<VAR>.pdf            all birds, species > VAR > latitude
#   evaladmix_LEPC_by_<VAR>.pdf       LEPC only, DPS > VAR > latitude
#   evaladmix_means_by_<VAR>.csv      mean residual within each species x VAR
#                                     cell, per K (all off-diagonal pairs)
#   evaladmix_LEPC_DPS.csv            LEPC: mean residual within Northern DPS,
#                                     within Southern DPS and between them, per K
#   evaladmix_top_pairs.csv           the highest residual pairs at the largest K
#                                     (duplicates / close relatives)
#   VAR = STATE, YEAR, SEX, plus HABITAT and GROUP (sympatry) for LEPC.
#
# Annotation columns are matched by NAME, not position. STATE comes from
# STATE_NAME (clean state names) when present, otherwise from STATE with any
# "County, " prefix dropped. SEX is reduced to Male / Female from entries such
# as "M; ASY; Mass(g)...". Missing values are kept and drawn grey ("NA").
#
# Reading the plots: a good fit leaves residual correlations near 0
# everywhere. A red block on the diagonal (positive correlation among a set of
# birds) with blue against the rest means those birds share ancestry the model
# does not explain -- structure K has not captured. Single dark-red cells are
# pairs of close relatives or duplicate samples. A group mean near 0 can hide
# such blocks (positive within sub-blocks, negative between them), so read the
# DPS table and the heatmaps, not only the means (Garcia-Erill & Albrechtsen
# 2020, Mol Ecol Resour 20:936-949).
# =============================================================================
library(readxl)

base    <- "/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files"
evaldir <- file.path(base, "evaladmix")      # where the .corres.txt files are
outdir  <- file.path(evaldir, "plots")
Ks      <- 1:10                              # K values to plot (missing files skipped)
zlim    <- 0.1                               # colour scale capped at +/- zlim
vars_all  <- c("STATE", "YEAR", "SEX")       # sorted / labelled for all birds
vars_lepc <- c("STATE", "HABITAT", "GROUP", "YEAR", "SEX")   # LEPC-only views
north_dps <- c("Kansas", "Nebraska", "Colorado")   # LEPC Northern DPS
south_dps <- c("New Mexico", "Texas")              # LEPC Southern DPS
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- metadata ----------------------------------------------------------------
md <- as.data.frame(read_xlsx(file.path(base, "nexus_metadata.xlsx")))
col <- function(name) { if (!name %in% names(md)) stop("metadata has no column ", name); md[[name]] }
blank_na <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA", "nan", "NaN")] <- NA; x }

species <- blank_na(col("common"))
sympat  <- blank_na(col("GROUP"))
lat     <- suppressWarnings(as.numeric(col("Lat")))
state   <- if ("STATE_NAME" %in% names(md)) blank_na(md$STATE_NAME) else
             blank_na(sub(".*,\\s*", "", col("STATE")))
sex_raw <- toupper(blank_na(col("SEX")))
sex     <- ifelse(grepl("^(M|MALE)\\b", sex_raw), "Male",
           ifelse(grepl("^(F|FEMALE)\\b", sex_raw), "Female", NA))
year    <- blank_na(col("YEAR"))
habitat <- blank_na(col("HABITAT")); habitat[species != "LEPC"] <- NA
dps     <- ifelse(species == "LEPC" & state %in% north_dps, "Northern DPS",
           ifelse(species == "LEPC" & state %in% south_dps, "Southern DPS", NA))
ids     <- if ("ID" %in% names(md)) as.character(md$ID) else as.character(seq_len(nrow(md)))
grp     <- ifelse(is.na(sympat), species, paste(species, sympat))
annot   <- list(STATE = state, YEAR = year, SEX = sex, HABITAT = habitat,
                GROUP = sympat, DPS = dps, species = species, group = grp)

# ---- helpers -----------------------------------------------------------------
read_corres <- function(f) {
  m <- as.matrix(read.table(f, na.strings = c("NA", "nan", "NaN", "-nan")))
  dimnames(m) <- NULL
  diag(m) <- NA
  m
}
pal <- colorRampPalette(c("#2B5C8A", "white", "#B5506B"))(101)
brk <- seq(-zlim, zlim, length.out = 102)
level_cols <- function(lv, shift = 0) {   # shift: start further along the palette
  base_cols <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948",
                 "#B07AA1", "#FF9DA7", "#9C755F", "#86BCB6", "#D37295", "#8CD17D",
                 "#499894", "#F1CE63", "#79706E")
  k <- shift %% length(base_cols)
  if (k > 0) base_cols <- c(base_cols[-seq_len(k)], base_cols[seq_len(k)])
  setNames(rep_len(base_cols, length(lv)), lv)
}
# order: each key in turn, NA last; numeric-looking keys sort numerically
ord_by <- function(keys) {
  k <- lapply(keys, function(x) {
    if (is.numeric(x)) return(x)
    xn <- suppressWarnings(as.numeric(x))
    if (all(is.na(xn) == is.na(x))) xn else x
  })
  do.call(order, c(k, list(na.last = TRUE)))
}

# draw one heatmap page: m = matrix in plotting order, strips = named list of
# label vectors (same order), bounds = vector whose changes get a dark line
plot_page <- function(m, strips, bounds, title, minlab = 3) {
  n  <- nrow(m)
  sw <- max(n * 0.025, 2)                          # strip width (matrix units)
  S  <- length(strips)
  cm <- matrix(pal[findInterval(pmax(pmin(m, zlim), -zlim), brk, all.inside = TRUE)], n)
  cm[is.na(m)] <- "grey90"
  par(mar = c(1, 1, 3, 1))
  plot.new()
  plot.window(xlim = c(-S * sw - n * 0.22, n + 0.5 + n * 0.30),
              ylim = c(-n * 0.02, n + 0.5 + S * sw + n * 0.12), asp = 1)
  rasterImage(as.raster(cm), 0.5, 0.5, n + 0.5, n + 0.5, interpolate = FALSE)
  # dark lines where `bounds` changes, light lines where the first strip changes
  cut_at <- function(v) { v[is.na(v)] <- "NA"; which(head(v, -1) != tail(v, -1)) }
  b1 <- cut_at(strips[[1]])
  segments(b1 + 0.5, 0.5, b1 + 0.5, n + 0.5, col = "grey55", lwd = 0.3)
  segments(0.5, n - b1 + 0.5, n + 0.5, n - b1 + 0.5, col = "grey55", lwd = 0.3)
  b0 <- cut_at(bounds)
  segments(b0 + 0.5, 0.5, b0 + 0.5, n + 0.5, col = "grey15", lwd = 0.8)
  segments(0.5, n - b0 + 0.5, n + 0.5, n - b0 + 0.5, col = "grey15", lwd = 0.8)
  # annotation strips: top (columns) and left (rows)
  legend_y <- n + 0.5
  for (s in seq_len(S)) {
    v  <- strips[[s]]; vv <- ifelse(is.na(v), "NA", v)
    lv <- unique(vv[vv != "NA"]); cols <- c(level_cols(sort(lv), shift = 9 * (s - 1)), "NA" = "grey80")
    y0 <- n + 0.5 + (s - 1) * sw; x0 <- -s * sw + 0.5
    rect(seq_len(n) - 0.5, y0, seq_len(n) + 0.5, y0 + sw * 0.9, col = cols[vv], border = NA)
    rect(x0, n - seq_len(n) + 0.5, x0 + sw * 0.9, n - seq_len(n) + 1.5, col = cols[vv], border = NA)
    text(-S * sw - n * 0.13, y0 + sw * 0.45, names(strips)[s], adj = 1, cex = 0.55, font = 2)
    # level labels of the FIRST strip (the sorting variable), placed outside
    # all strips (top: rotated; left: horizontal); runs < minlab birds unlabelled
    if (s == 1) {
      r <- rle(vv); e <- cumsum(r$lengths); st <- e - r$lengths + 1
      keep <- r$lengths >= minlab
      text((st + e)[keep] / 2, n + 0.5 + S * sw + sw * 0.2, r$values[keep], srt = 90, adj = 0, cex = 0.45)
      text(-S * sw + 0.5 - sw * 0.2, n - (st + e)[keep] / 2 + 1, r$values[keep], adj = 1, cex = 0.45)
    }
    # legend per strip on the right
    lg <- legend(n + 0.5 + n * 0.02, legend_y, legend = names(cols), fill = cols,
                 border = NA, bty = "n", cex = 0.5, title = names(strips)[s],
                 title.adj = 0, y.intersp = 0.9, xpd = NA)
    legend_y <- lg$rect$top - lg$rect$h - n * 0.01
  }
  # colour key
  kx <- n + 0.5 + n * 0.22; ky <- seq(n * 0.05, n * 0.35, length.out = 102)
  rect(kx, head(ky, -1), kx + n * 0.025, tail(ky, -1), col = pal, border = NA)
  text(kx + n * 0.03, c(min(ky), mean(ky), max(ky)), c(-zlim, 0, zlim), adj = 0, cex = 0.5)
  text(kx, max(ky) + n * 0.02, "residual r", adj = 0, cex = 0.5)
  title(title, cex.main = 0.9)
}

# mean of off-diagonal pairs within each cell of `cell` (vector of labels)
cell_means <- function(m, cell) {
  cl <- ifelse(is.na(cell), "NA", cell)
  do.call(rbind, lapply(unique(cl), function(g) {
    i <- which(cl == g); v <- m[i, i][upper.tri(m[i, i])]
    data.frame(cell = g, n = length(i), mean_r = if (length(v)) mean(v, na.rm = TRUE) else NA)
  }))
}

# ---- read matrices -----------------------------------------------------------
mats <- list()
for (K in Ks) {
  f <- file.path(evaldir, paste0("nexus_k", K, ".corres.txt"))
  if (!file.exists(f)) { message("skip K=", K, " (no ", basename(f), ")"); next }
  m <- read_corres(f)
  if (nrow(m) != length(species) || ncol(m) != length(species))
    stop(basename(f), " is ", nrow(m), "x", ncol(m), " but metadata has ", length(species), " rows")
  mats[[as.character(K)]] <- m
}
if (!length(mats)) stop("no .corres.txt files found in ", evaldir)
Kfound <- as.integer(names(mats))

# ---- 1. overview (admixture barplot order) -----------------------------------
o <- ord_by(list(species, sympat, -lat))
pdf(file.path(outdir, "evaladmix_overview.pdf"), width = 9, height = 8.5)
for (K in Kfound) plot_page(mats[[as.character(K)]][o, o],
                            list(GROUP = sympat[o], species = species[o]), species[o],
                            paste0("evalAdmix residual correlation, all birds, K = ", K))
dev.off()

# ---- 2. all birds, sorted by each variable within species --------------------
for (v in vars_all) {
  x <- annot[[v]]; o <- ord_by(list(species, x, -lat))
  pdf(file.path(outdir, paste0("evaladmix_by_", v, ".pdf")), width = 9, height = 8.5)
  for (K in Kfound) plot_page(mats[[as.character(K)]][o, o],
                              setNames(list(x[o], species[o]), c(v, "species")), species[o],
                              paste0("All birds by species > ", v, ", K = ", K))
  dev.off()
  tab <- do.call(rbind, lapply(Kfound, function(K)
    cbind(K = K, cell_means(mats[[as.character(K)]], paste(species, ifelse(is.na(x), "NA", x), sep = " | ")))))
  write.csv(tab, file.path(outdir, paste0("evaladmix_means_by_", v, ".csv")), row.names = FALSE)
}

# ---- 3. LEPC only, sorted by DPS then each variable --------------------------
L <- which(species == "LEPC")
for (v in vars_lepc) {
  x <- annot[[v]][L]; o <- L[ord_by(list(dps[L], x, -lat[L]))]; xo <- annot[[v]][o]
  pdf(file.path(outdir, paste0("evaladmix_LEPC_by_", v, ".pdf")), width = 9, height = 8.5)
  for (K in Kfound) plot_page(mats[[as.character(K)]][o, o],
                              setNames(list(xo, dps[o]), c(v, "DPS")), dps[o],
                              paste0("LEPC by DPS > ", v, ", K = ", K))
  dev.off()
}

# ---- 4. LEPC DPS table -----------------------------------------------------------
N <- which(dps == "Northern DPS"); S <- which(dps == "Southern DPS")
dps_tab <- do.call(rbind, lapply(Kfound, function(K) {
  m <- mats[[as.character(K)]]
  data.frame(K = K,
             within_north = mean(m[N, N][upper.tri(m[N, N])], na.rm = TRUE),
             within_south = mean(m[S, S][upper.tri(m[S, S])], na.rm = TRUE),
             between      = mean(m[N, S], na.rm = TRUE))
}))
write.csv(dps_tab, file.path(outdir, "evaladmix_LEPC_DPS.csv"), row.names = FALSE)
cat("\nLEPC: mean residual r within and between DPS (n North =", length(N),
    ", n South =", length(S), ")\n"); print(dps_tab, digits = 3, row.names = FALSE)

# ---- 5. highest-residual pairs at the largest K ------------------------------
m <- mats[[as.character(max(Kfound))]]
up <- which(upper.tri(m), arr.ind = TRUE)
r <- m[up]; top <- head(order(r, decreasing = TRUE), 25)
pairs <- data.frame(K = max(Kfound), id1 = ids[up[top, 1]], id2 = ids[up[top, 2]],
                    species1 = species[up[top, 1]], species2 = species[up[top, 2]],
                    state1 = state[up[top, 1]], state2 = state[up[top, 2]],
                    year1 = year[up[top, 1]], year2 = year[up[top, 2]], r = round(r[top], 3))
write.csv(pairs, file.path(outdir, "evaladmix_top_pairs.csv"), row.names = FALSE)
cat("\nHighest residual pairs at K =", max(Kfound), "(r > 0.4 suggests duplicates / first-degree relatives):\n")
print(head(pairs, 10), row.names = FALSE)
message("\nWrote plots and tables to ", outdir)
