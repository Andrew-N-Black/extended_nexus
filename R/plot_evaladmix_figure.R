# =============================================================================
# plot_evaladmix_figure.R -- report Figure S5 (USFWS report, Objective 2):
# evalAdmix fit of the best NGSadmix run per K.
#   (A) mean absolute residual correlation over all pairs of birds, K = 1-10
#   (B) LEPC: mean residual within the Northern DPS, within the Southern DPS
#       and between them, K = 1-10
#   (C, D) LEPC residual-correlation matrices at K = 3 and K = 4, birds ordered
#       by DPS, then ecoregion, then latitude (north to south)
# Inputs as for plot_evaladmix.R (evaladmix/nexus_k<K>.corres.txt and
# nexus_metadata.xlsx in samples.txt row order).
# =============================================================================
library(readxl)

base    <- "/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files"
evaldir <- file.path(base, "evaladmix")
outfile <- file.path(evaldir, "plots", "FigS_evaladmix")   # .png and .pdf
Ks      <- 1:10
Kheat   <- c(3, 4)
zlim    <- 0.1
north_dps <- c("Kansas", "Nebraska", "Colorado")
south_dps <- c("New Mexico", "Texas")
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

md <- as.data.frame(read_xlsx(file.path(base, "nexus_metadata.xlsx")))
blank_na <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA", "nan")] <- NA; x }
species <- blank_na(md$common)
state   <- blank_na(md$STATE_NAME)
habitat <- blank_na(md$HABITAT)
lat     <- suppressWarnings(as.numeric(md$Lat))
dps     <- ifelse(state %in% north_dps, "Northern", ifelse(state %in% south_dps, "Southern", NA))

read_corres <- function(K) {
  m <- as.matrix(read.table(file.path(evaldir, paste0("nexus_k", K, ".corres.txt")),
                            na.strings = c("NA", "nan", "NaN", "-nan")))
  dimnames(m) <- NULL; diag(m) <- NA
  stopifnot(nrow(m) == nrow(md))
  m
}
mats <- lapply(Ks, read_corres); names(mats) <- Ks

L <- which(species == "LEPC" & !is.na(dps))
N <- which(species == "LEPC" & dps == "Northern")
S <- which(species == "LEPC" & dps == "Southern")
upm <- function(x) mean(x[upper.tri(x)], na.rm = TRUE)
stats <- data.frame(
  K        = Ks,
  mean_abs = sapply(mats, function(m) mean(abs(m[upper.tri(m)]), na.rm = TRUE)),
  north    = sapply(mats, function(m) upm(m[N, N])),
  south    = sapply(mats, function(m) upm(m[S, S])),
  between  = sapply(mats, function(m) mean(m[N, S], na.rm = TRUE)))
write.csv(stats, paste0(outfile, "_data.csv"), row.names = FALSE)
print(stats, digits = 3)

# ---- colours -----------------------------------------------------------------
pal  <- colorRampPalette(c("#2B5C8A", "white", "#B5506B"))(101)
brk  <- seq(-zlim, zlim, length.out = 102)
dcol <- c(Northern = "#7A9E5A", Southern = "#E0A33B")
hab_lv  <- c("SHGPR", "MGPR", "SSBPR", "SSOPR")
hab_lab <- c(SHGPR = "Short-grass prairie/CRP", MGPR = "Mixed-grass prairie",
             SSBPR = "Sand sagebrush prairie", SSOPR = "Shinnery oak prairie")
hcol <- c(SHGPR = "#4E79A7", MGPR = "#76B7B2", SSBPR = "#B07AA1", SSOPR = "#9C755F", "NA" = "grey80")

o  <- L[order(dps[L], factor(habitat[L], levels = hab_lv), -lat[L], na.last = TRUE)]
n  <- length(o)
hv <- ifelse(is.na(habitat[o]), "NA", habitat[o])

heat <- function(K, lab) {
  m  <- mats[[as.character(K)]][o, o]
  cm <- matrix(pal[findInterval(pmax(pmin(m, zlim), -zlim), brk, all.inside = TRUE)], n)
  cm[is.na(m)] <- "grey90"
  sw <- n * 0.03
  par(mar = c(0.5, 0.5, 2, 0.5))
  plot.new(); plot.window(xlim = c(-2.2 * sw, n + 0.5), ylim = c(0.5, n + 0.5 + 2.2 * sw), asp = 1)
  rasterImage(as.raster(cm), 0.5, 0.5, n + 0.5, n + 0.5, interpolate = FALSE)
  for (s in 1:2) {
    v  <- if (s == 1) hcol[hv] else dcol[dps[o]]
    y0 <- n + 0.5 + (s - 1) * sw * 1.1; x0 <- 0.5 - s * sw * 1.1
    rect(seq_len(n) - 0.5, y0, seq_len(n) + 0.5, y0 + sw, col = v, border = NA)
    rect(x0, n - seq_len(n) + 0.5, x0 + sw, n - seq_len(n) + 1.5, col = v, border = NA)
  }
  b <- which(head(dps[o], -1) != tail(dps[o], -1))
  segments(b + 0.5, 0.5, b + 0.5, n + 0.5, lwd = 0.8)
  segments(0.5, n - b + 0.5, n + 0.5, n - b + 0.5, lwd = 0.8)
  mtext(lab, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  mtext(paste0("LEPC, K = ", K), side = 3, line = 0.6, cex = 0.7)
}

panel_lab <- function(x) {
  u <- par("usr"); mtext(x, side = 3, line = 0.6, at = u[1] - 0.14 * diff(u[1:2]), adj = 0, font = 2, cex = 0.8)
}

draw <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, byrow = TRUE), heights = c(1, 1.35, 0.25))
  # A
  par(mar = c(3.2, 3.8, 2, 1), mgp = c(2, 0.6, 0), cex.axis = 0.75, cex.lab = 0.8, las = 1)
  plot(stats$K, stats$mean_abs, type = "b", pch = 19, cex = 0.7, xaxt = "n",
       xlab = "K", ylab = "Mean |residual r|", ylim = c(0, max(stats$mean_abs) * 1.08))
  axis(1, at = Ks)
  panel_lab("A")
  # B
  yl <- range(stats[, c("north", "south", "between")]) * 1.1
  plot(stats$K, stats$north, type = "b", pch = 19, cex = 0.7, col = dcol["Northern"],
       xaxt = "n", ylim = yl, xlab = "K", ylab = "Mean residual r (LEPC)")
  lines(stats$K, stats$south, type = "b", pch = 19, cex = 0.7, col = dcol["Southern"])
  lines(stats$K, stats$between, type = "b", pch = 1, cex = 0.7, col = "grey30", lty = 2)
  abline(h = 0, col = "grey70", lwd = 0.5); axis(1, at = Ks)
  legend("topright", c("Within Northern DPS", "Within Southern DPS", "Between DPSs"),
         col = c(dcol, "grey30"), pch = c(19, 19, 1), lty = c(1, 1, 2), bty = "n", cex = 0.7)
  panel_lab("B")
  # C, D
  heat(Kheat[1], "C"); heat(Kheat[2], "D")
  # legend row
  par(mar = c(0, 0.5, 0, 0.5)); plot.new()
  legend(0, 0.95, c("Northern DPS", "Southern DPS"), fill = dcol, border = NA, bty = "n",
         cex = 0.7, title = "DPS (outer bar)", title.adj = 0)
  legend(0.22, 0.95, c(hab_lab[hab_lv], "Not recorded"), fill = hcol[c(hab_lv, "NA")], border = NA, bty = "n",
         cex = 0.7, ncol = 2, title = "Ecoregion (inner bar)", title.adj = 0)
  kx <- seq(0.78, 0.98, length.out = 102)
  rect(head(kx, -1), 0.45, tail(kx, -1), 0.65, col = pal, border = NA)
  text(c(0.78, 0.88, 0.98), 0.32, c(-zlim, 0, zlim), cex = 0.65)
  text(0.88, 0.82, "Residual correlation", cex = 0.7)
}

png(paste0(outfile, ".png"), width = 7, height = 8, units = "in", res = 400)
draw(); dev.off()
pdf(paste0(outfile, ".pdf"), width = 7, height = 8)
draw(); dev.off()
message("Wrote ", outfile, ".png/.pdf")
