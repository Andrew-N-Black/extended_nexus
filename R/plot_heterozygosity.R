#Load libraries
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#Load libraries
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_nexus_new_plus_shotguns/grouse_extended_nexus_samples_v2.xlsx")
#Plot
ggplot(HET_FILT, aes(x=SPECIES, y=Heterozygosity, fill=SPECIES)) +geom_boxplot() +scale_fill_manual("", values=c("Tympanuchus cupido" = "goldenrod","Tympanuchus pallidicinctus"="brown","Tympanuchus phasianellus" = "black","Tympanuchus phasianellus/Tympanuchus cupido" = "grey")) + xlab("") + ylab("H") +theme_classic() +theme(axis.text.y = element_text(size=12)) +theme(legend.position="bottom") +theme(axis.text=element_text(size=14), axis.title=element_text(size=22, face="italic")) +theme(strip.text = element_text(size=18))+facet_wrap(~GROUP)+theme(axis.text.x = element_blank(), axis.title.x = element_blank())



#Test for normality
shapiro.test(HET_FILT$Heterozygosity)

	Shapiro-Wilk normality test

data:  HET_FILT$Heterozygosity
W = 0.83168, p-value < 2.2e-16

pairwise.wilcox.test(HET_FILT$Heterozygosity, HET_FILT$SPECIES, p.adjust.method = "BH")

	Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

data:  HET_FILT$Heterozygosity and HET_FILT$SPECIES 

                                            Tympanuchus cupido Tympanuchus pallidicinctus Tympanuchus phasianellus
Tympanuchus pallidicinctus                  1.4e-10            -                          -                       
Tympanuchus phasianellus                    1.9e-12            0.012                      -                       
Tympanuchus phasianellus/Tympanuchus cupido 1.000              0.193                      0.041    

library(dplyr)
library(rstatix)   # for effect sizes; install.packages("rstatix") if needed

# n, median and IQR per group
HET_FILT %>%
    group_by(SPECIES) %>%
    summarise(n = n(), median_H = median(Heterozygosity),
              IQR_H = IQR(Heterozygosity), mean_H = mean(Heterozygosity),
              sd_H = sd(Heterozygosity))

# Omnibus test, then pairwise tests on the three species only
kruskal.test(Heterozygosity ~ SPECIES, data = HET_FILT)
sp3 <- subset(HET_FILT, SPECIES != "Tympanuchus phasianellus/Tympanuchus cupido")
kruskal.test(Heterozygosity ~ SPECIES, data = sp3)
pairwise.wilcox.test(sp3$Heterozygosity, sp3$SPECIES, p.adjust.method = "BH")
wilcox_effsize(sp3, Heterozygosity ~ SPECIES)   # r: ~0.1 small, ~0.3 medium, ~0.5 large



wilcox_effsize(sp3, Heterozygosity ~ SPECIES)   # r: ~0.1 small, ~0.3 medium, ~0.5 large
# A tibble: 3 × 7
  .y.            group1                     group2                     effsize    n1    n2 magnitude
* <chr>          <chr>                      <chr>                        <dbl> <int> <int> <ord>    
1 Heterozygosity Tympanuchus cupido         Tympanuchus pallidicinctus   0.302    49   426 moderate 
2 Heterozygosity Tympanuchus cupido         Tympanuchus phasianellus     0.826    49    29 large    
3 Heterozygosity Tympanuchus pallidicinctus Tympanuchus phasianellus     0.129   426    29 small    
