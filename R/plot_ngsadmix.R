# =============================================================================
# plot_ngsadmix.R -- USFWS report Objective 2, Figure 11: NGSadmix ancestry
# proportions, K = 2-4, stacked into one figure (pophelper).
# Inputs: nexus_k{2,3,4}.qopt (best run per K from analysis/nexus_ngsadmix.sh,
# rows in beagle sample order) and nexus_metadata.xlsx in the SAME row order.
# =============================================================================
library(readxl)
library(pophelper)

base <- "/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files"

metadata <- read_xlsx(file.path(base, "nexus_metadata.xlsx"))

species  <- as.character(metadata[[32]])
sympatry <- as.character(metadata[[14]])
lat      <- as.numeric(metadata[[11]])     # Latitude (col 12 = longitude, unused)

# ---- one combined label row (e.g. "LEPC Sympatric") ----
# Two stacked rows of angled labels collide in pophelper, so merge them.
# Samples with no sympatry value (e.g. STGR / GRPC) keep the species name only.
grp <- ifelse(is.na(sympatry) | sympatry == "",
              species,
              paste(species, sympatry))

# ---- sample order: by species, then sympatry, then latitude N -> S ----
ord <- order(species, sympatry, -lat, na.last = TRUE)

# ---- read all K, align clusters across K, then apply the order ----
qfiles <- file.path(base, paste0("nexus_k", 2:4, ".qopt"))
qlist  <- readQ(files = qfiles)
qlist  <- alignK(qlist, type = "across")
qlist  <- lapply(qlist, function(q) {
    q2 <- q[ord, , drop = FALSE]
    rownames(q2) <- NULL
    q2
})

grplab <- data.frame(" " = grp[ord], check.names = FALSE,
                     stringsAsFactors = FALSE)

# ---- colour palette (muted, colourblind-safe) ----
cols <- c("#2B5C8A",   # deep blue
          "#E0A33B",   # ochre
          "#7A9E5A",   # sage green
          "#B5506B")   # muted rose

plotQ(qlist,
      imgoutput     = "join",
      returnplot    = TRUE,
      exportplot    = TRUE,
      clustercol    = cols,
      grplab        = grplab,
      ordergrp      = FALSE,        # order already set above (keeps latitude sort)
      sharedindlab  = FALSE,
      showlegend    = FALSE,
      showsp        = TRUE,
      splab         = paste0("K = ", 2:4),
      splabsize     = 6,
      splabcol      = "black",
      spbgcol       = "white",
      barbordersize   = 0,
      barbordercolour = NA,
      barsize         = 1,
      divtype   = 1,
      divcol    = "white",
      divsize   = 0.2,
      grplabsize   = 1.8,
      grplabangle  = 45,
      grplabjust   = 1,
      grplabheight = 7,
      linepos      = 0.95,
      grplabpos    = 0.85,
      linesize     = 0.4,
      pointsize    = 1.5,
      width  = 24,
      height = 2,
      dpi    = 600,
      outputfilename = "grouse_nexus_K2-4",
      imgtype  = "pdf",
      exportpath = getwd())
