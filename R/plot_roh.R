library(ggplot2)
library(dplyr)
library(reshape2)

## ---- test data (delete this block; use your own `sub`) ----
if (!is.data.frame(get0("sub"))) {
  set.seed(1)
  n <- c(LEPC = 426, GRPC = 49, STGR = 29, "STGR / GRPC" = 2)
  sub <- data.frame(ID = paste0("F", 1:506),
                    common = rep(names(n), n))
  sub$GROUP <- ifelse(runif(506) < 0.32, "Sympatric", "Allopatric")
  sub$fROH_100kb <- rbeta(506, 2, 40)
  sub$fROH_1Mb   <- rbeta(506, 0.3, 60)
  sub$fROH_total <- sub$fROH_100kb + sub$fROH_1Mb
}
## ------------------------------------------------------------

melt_data <- melt(sub, id.vars = c("ID", "GROUP", "common"))

# readable facet labels: matched on the column name, whatever it is exactly
lab <- function(v) ifelse(grepl("tot", v, ignore.case = TRUE), "Total",
                   ifelse(grepl("100", v), "100 kb-1 Mb", "> 1 Mb"))
melt_data$class <- factor(lab(as.character(melt_data$variable)),
                          levels = c("100 kb-1 Mb", "> 1 Mb", "Total"))
melt_data$GROUP  <- factor(melt_data$GROUP, levels = c("Allopatric", "Sympatric"))
melt_data$common <- factor(melt_data$common, levels = c("LEPC", "GRPC", "STGR", "STGR / GRPC"))

# one order for every row: each bird keeps the same x position in all three
# panels, sorted by its TOTAL fROH within its sampling context
ord <- melt_data %>% filter(class == "Total") %>% arrange(GROUP, value) %>% pull(ID)
# (alternative: group birds by species first, then by total fROH)
# ord <- melt_data %>% filter(class == "Total") %>% arrange(GROUP, common, value) %>% pull(ID)
melt_data$ID <- factor(melt_data$ID, levels = unique(ord))

n_grp <- table(sub$GROUP)
grp_lab <- setNames(sprintf("%s (n = %d)", names(n_grp), n_grp), names(n_grp))

p <- ggplot(melt_data, aes(x = ID, y = value, fill = common)) +
  geom_col(width = 1) +                                  # no gaps between bars
  facet_grid(class ~ GROUP,
             scales = "free", space = "free_x",          # each context shows only its own birds
             labeller = labeller(GROUP = grp_lab)) +
  scale_fill_manual(NULL, values = c(LEPC = "brown", GRPC = "goldenrod",
                                     STGR = "black", "STGR / GRPC" = "grey60")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Individuals (n = 506)", y = expression(italic(f)[ROH])) +
  theme_classic(base_size = 14) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line        = element_blank(),
        strip.background = element_rect(colour = "black", fill = "grey90", linewidth = 0.6),
        strip.text       = element_text(size = 14),
        strip.text.y     = element_text(angle = 0, size = 12),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(size = 11),
        axis.title.x     = element_text(size = 14),
        panel.spacing    = unit(0.4, "lines"),
        legend.position  = "bottom")
print(p)
ggsave("nexus_fROH.png", p, width = 9, height = 6.5, dpi = 600)
ggsave("nexus_fROH.pdf", p, width = 9, height = 6.5, device = cairo_pdf)  # cairo handles non-ASCII text
