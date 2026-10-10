# =============================================================================
# fROH summaries and tests for the nexus panel, by species ("common") and
# sampling context ("GROUP").  Input: `sub` with columns ID, GROUP, common and
# the three fROH columns (100 kb-1 Mb, >1 Mb, total), as used for the plot.
# USFWS report, Objective 2: f_ROH Results, Figure 9; depth check (section 5).
# Needs: dplyr, tidyr, rstatix   install.packages(c("dplyr","tidyr","rstatix"))
# =============================================================================
library(dplyr)
library(tidyr)
library(rstatix)

## ---- input ------------------------------------------------------------------
## `sub`: one row per bird with columns ID, GROUP (Allopatric/Sympatric),
## common (LEPC, GRPC, STGR, "STGR / GRPC") and the three fROH columns
## (100 kb-1 Mb, >1 Mb, total) from roh_parse_autosomal.sh, merged with the
## sample metadata (nexus_metadata.xlsx).
if (!is.data.frame(get0("sub")))
    stop("Create `sub` first (ID, GROUP, common, three fROH columns); see header.")
## ---------------------------------------------------------------------------

DEPTH_COL <- "depth"   # optional mean-depth column in `sub` (section 5); change to your name
HYB <- "STGR / GRPC"                     # exact label in your data (with spaces)

# long format with readable class names (matched on column name)
roh <- sub %>%
    pivot_longer(-any_of(c("ID", "GROUP", "common", DEPTH_COL)), names_to = "variable", values_to = "fROH") %>%
    mutate(class  = case_when(grepl("tot", variable, ignore.case = TRUE) ~ "Total",
                              grepl("100", variable)                     ~ "100 kb-1 Mb",
                              TRUE                                       ~ ">1 Mb"),
           class  = factor(class,  levels = c("100 kb-1 Mb", ">1 Mb", "Total")),
           common = factor(common, levels = c("LEPC", "GRPC", "STGR", HYB)),
           GROUP  = factor(GROUP,  levels = c("Allopatric", "Sympatric")))

summ <- function(d) summarise(d, n = n(),
                              median = median(fROH), IQR = IQR(fROH), mean = mean(fROH), sd = sd(fROH),
                              min = min(fROH), max = max(fROH),
                              n_nonzero = sum(fROH > 0), .groups = "drop")   # birds with any ROH in that class

# ---- 1. Descriptive summaries ----------------------------------------------
by_species        <- roh %>% group_by(class, common)        %>% summ()
by_group          <- roh %>% group_by(class, GROUP)         %>% summ()
by_species_group  <- roh %>% group_by(class, common, GROUP) %>% summ()
print(by_species, n = Inf); print(by_species_group, n = Inf)

# ---- 2. Species differences (hybrids excluded: n = 2) -----------------------
sp <- roh %>% filter(common != HYB) %>% droplevels()

kw_species <- sp %>% group_by(class) %>% kruskal_test(fROH ~ common)
pw_species <- sp %>% group_by(class) %>%
    pairwise_wilcox_test(fROH ~ common, p.adjust.method = "BH")
es_species <- sp %>% group_by(class) %>% wilcox_effsize(fROH ~ common)

# same, separately within each sampling context
kw_species_by_group <- sp %>% group_by(class, GROUP) %>% kruskal_test(fROH ~ common)
pw_species_by_group <- sp %>% group_by(class, GROUP) %>%
    pairwise_wilcox_test(fROH ~ common, p.adjust.method = "BH")

# ---- 3. Allopatric vs sympatric, within each species ------------------------
# (only species sampled in both contexts; hybrids are all sympatric)
both <- sp %>% group_by(common) %>% filter(n_distinct(GROUP) == 2) %>% ungroup()
wx_context <- both %>% group_by(class, common) %>%
    wilcox_test(fROH ~ GROUP) %>%
    adjust_pvalue(method = "BH") %>% add_significance()
es_context <- both %>% group_by(class, common) %>% wilcox_effsize(fROH ~ GROUP)

# and pooled across species (descriptive; confounded by species composition)
wx_context_all <- sp %>% group_by(class) %>% wilcox_test(fROH ~ GROUP)

# ---- 4. Print and save ------------------------------------------------------
for (nm in c("kw_species", "pw_species", "es_species", "kw_species_by_group",
             "pw_species_by_group", "wx_context", "es_context", "wx_context_all")) {
    cat("\n====", nm, "====\n"); print(as.data.frame(get(nm)))
}
write.csv(by_species_group, "fROH_summary_species_by_group.csv", row.names = FALSE)
write.csv(by_species,       "fROH_summary_species.csv",          row.names = FALSE)
write.csv(as.data.frame(pw_species), "fROH_pairwise_species.csv", row.names = FALSE)
write.csv(as.data.frame(wx_context), "fROH_allo_vs_sym.csv",      row.names = FALSE)

# ---- 5. Does f_ROH track sequencing depth? -----------------------------------
# Needs a numeric mean-depth column in `sub` (post-harmonization depth from
# nexus_downsample_depth_mapping.sh). Two hybrid rows in nexus_metadata.xlsx
# are column-shifted, so coerce to numeric (non-numeric -> NA, dropped).
if (DEPTH_COL %in% names(sub)) {
    d <- sub %>% filter(common != HYB) %>%
        mutate(depth = suppressWarnings(as.numeric(.data[[DEPTH_COL]])),
               fROH_tot = as.numeric(.data[[grep("tot", names(sub), ignore.case = TRUE, value = TRUE)[1]]])) %>%
        filter(is.finite(depth), is.finite(fROH_tot))
    cat("\n==== Spearman: total fROH vs depth ====\n")
    print(cor.test(d$fROH_tot, d$depth, method = "spearman", exact = FALSE))
    for (s in c("LEPC", "GRPC", "STGR")) {
        x <- filter(d, common == s)
        ct <- cor.test(x$fROH_tot, x$depth, method = "spearman", exact = FALSE)
        cat(sprintf("%s: rho = %.2f, P = %.2g, n = %d\n", s, ct$estimate, ct$p.value, nrow(x)))
    }
    # LEPC: ranked total fROH ~ sampling context + depth
    lepc <- filter(d, common == "LEPC") %>% mutate(r = rank(fROH_tot))
    print(summary(lm(r ~ GROUP + depth, data = lepc)))
} else message("Section 5 skipped: no '", DEPTH_COL, "' column in `sub`.")
