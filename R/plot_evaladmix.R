# =============================================================================
# plot_evaladmix.R -- evalAdmix residual-correlation heatmaps for the NGSadmix
# runs (USFWS report Objective 2; choosing K for Figure 11).
#
# Inputs (copy from ${OUT} on the cluster; see analysis/nexus_ngsadmix.sh):
#   evaladmix/nexus_k<K>.corres.txt   one per K, rows/cols in beagle sample order
#   nexus_metadata.xlsx               SAME row order as samples.txt (as for
#                                     plot_ngsadmix.R)
# Outputs (working directory):
#   evaladmix_residuals_K<..>.pdf     one heatmap per K, same sample order and
#                                     group labels as the admixture barplot
#   evaladmix_group_means.csv         mean residual within each group, per K
#
# Reading the plots: a good fit leaves residual correlations near 0
# everywhere. A red block on the diagonal (positive correlation among members
# of one group) means that group is not modelled well at that K -- typically
# K is too low for it. Blue between groups is expected where the model forces
# a shared ancestry component. Choose the smallest K beyond which the
# within-group blocks stop shrinking (Garcia-Erill & Albrechtsen 2020).
# =============================================================================
library(readxl)

base    <- "/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files"
evaldir <- file.path(base, "evaladmix")      # where the .corres.txt files are
Ks      <- 1:6                               # K values to plot (files that exist)
zlim    <- 0.1                               # colour scale capped at +/- zlim

metadata <- read_xlsx(file.path(base, "nexus_metadata.xlsx"))
species  <- as.character(metadata[[32]])
sympatry <- as.character(metadata[[14]])
lat      <- as.numeric(metadata[[11]])
grp <- ifelse(is.na(sympatry) | sympatry == "", species, paste(species, sympatry))
ord <- order(species, sympatry, -lat, na.last = TRUE)   # same order as plot_ngsadmix.R
g   <- grp[ord]

read_corres <- function(f) {
  m <- as.matrix(read.table(f, na.strings = c("NA", "nan", "NaN")))
  dimnames(m) <- NULL
  m
}

pal <- colorRampPalette(c("#2B5C8A", "white", "#B5506B"))(101)
brk <- seq(-zlim, zlim, length.out = 102)

# group boundaries and label positions along the ordered axis
runs   <- rle(g)
ends   <- cumsum(runs$lengths)
starts <- c(1, head(ends, -1) + 1)
mids   <- (starts + ends) / 2

means <- list()
files_found <- c()
for (K in Ks) {
  f <- file.path(evaldir, paste0("nexus_k", K, ".corres.txt"))
  if (!file.exists(f)) { message("skip K=", K, " (no ", basename(f), ")"); next }
  m <- read_corres(f)
  stopifnot(nrow(m) == length(grp), ncol(m) == length(grp))
  m <- m[ord, ord]
  diag(m) <- NA

  # mean residual within each group (off-diagonal pairs only)
  for (gg in unique(g)) {
    i <- which(g == gg)
    v <- m[i, i][upper.tri(m[i, i])]
    means[[length(means) + 1]] <- data.frame(K = K, group = gg, n = length(i),
                                             mean_r = if (length(v)) mean(v, na.rm = TRUE) else NA)
  }

  mm <- pmax(pmin(m, zlim), -zlim)            # cap for colouring
  out <- paste0("evaladmix_residuals_K", K, ".pdf")
  pdf(out, width = 7.5, height = 7)
  par(mar = c(9, 9, 3, 6))
  n <- nrow(mm)
  image(1:n, 1:n, mm[, n:1], col = pal, breaks = brk, axes = FALSE,
        xlab = "", ylab = "", useRaster = TRUE,
        main = paste0("evalAdmix residual correlation, K = ", K))
  abline(v = ends[-length(ends)] + 0.5, h = n - ends[-length(ends)] + 0.5,
         col = "grey30", lwd = 0.4)
  axis(1, at = mids, labels = runs$values, las = 2, cex.axis = 0.6, tick = FALSE)
  axis(2, at = n - mids + 1, labels = runs$values, las = 2, cex.axis = 0.6, tick = FALSE)
  # colour key
  par(new = TRUE, fig = c(0.9, 0.93, 0.35, 0.75), mar = c(0, 0, 0, 0))
  image(1, seq(-zlim, zlim, length.out = 101), t(seq(-zlim, zlim, length.out = 101)),
        col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.6)
  dev.off()
  files_found <- c(files_found, out)
}

if (length(means)) {
  tab <- do.call(rbind, means)
  write.csv(tab, "evaladmix_group_means.csv", row.names = FALSE)
  print(reshape(tab[, c("K", "group", "mean_r")], idvar = "group",
                timevar = "K", direction = "wide"), digits = 3)
}
message("Wrote: ", paste(files_found, collapse = ", "))
