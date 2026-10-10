# =============================================================================
# fROH summaries and tests for the nexus panel, by species ("common") and
# sampling context ("GROUP").  Input: `sub` with columns ID, GROUP, common and
# the three fROH columns (100 kb-1 Mb, >1 Mb, total), as used for the plot.
# Needs: dplyr, tidyr, rstatix   install.packages(c("dplyr","tidyr","rstatix"))
# =============================================================================
library(dplyr)
library(tidyr)
library(rstatix)

## ---- test data (ignored when your own `sub` data frame exists) -------------
if (!is.data.frame(get0("sub"))) {
    set.seed(1)
    n <- c(LEPC = 426, GRPC = 49, STGR = 29, "STGR / GRPC" = 2)
    sub <- data.frame(ID = paste0("F", 1:506), common = rep(names(n), n))
    sub$GROUP <- ifelse(runif(506) < 0.32, "Sympatric", "Allopatric")
    sub$GROUP[sub$common == "STGR / GRPC"] <- "Sympatric"
    sub$fROH_100kb <- rbeta(506, 2, 40); sub$fROH_1Mb <- rbeta(506, 0.3, 60)
    sub$fROH_total <- sub$fROH_100kb + sub$fROH_1Mb
}
## ---------------------------------------------------------------------------

HYB <- "STGR / GRPC"                     # exact label in your data (with spaces)

# long format with readable class names (matched on column name)
roh <- sub %>%
    pivot_longer(-c(ID, GROUP, common), names_to = "variable", values_to = "fROH") %>%
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
