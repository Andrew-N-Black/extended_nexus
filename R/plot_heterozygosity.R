#Load libraries
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files/nexus_metadata.xlsx")
#Plot
ggplot(HET_FILT, aes(x=SPECIES, y=Heterozygosity, fill=SPECIES)) +geom_boxplot() +scale_fill_manual("", values=c("Tympanuchus cupido" = "goldenrod","Tympanuchus pallidicinctus"="brown","Tympanuchus phasianellus" = "black","Tympanuchus phasianellus/Tympanuchus cupido" = "grey")) + xlab("") + ylab("H") +theme_classic() +theme(axis.text.y = element_text(size=12)) +theme(legend.position="bottom") +theme(axis.text=element_text(size=14), axis.title=element_text(size=22, face="italic")) +theme(strip.text = element_text(size=18))+facet_wrap(~GROUP)+theme(axis.text.x = element_blank(), axis.title.x = element_blank())



#Test for normality
shapiro.test(HET_FILT$Heterozygosity)

#data:  HET_FILT$Heterozygosity
#W = 0.77728, p-value < 2.2e-16

pairwise.wilcox.test(HET_FILT$Heterozygosity, HET_FILT$common, p.adjust.method = "BH")

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  HET_FILT$Heterozygosity and HET_FILT$common 

#GRPC    LEPC   STGR  
#LEPC          1.6e-08  -       -     
 # STGR        6.8e-14  0.0095  -     
 # STGR / GRPC 0.5663   0.1538  0.0129



library(dplyr)
library(rstatix)  

# n, median and IQR per group
HET_FILT %>%
  group_by(common) %>%
  summarise(n = n(), median_H = median(Heterozygosity),
            IQR_H = IQR(Heterozygosity), mean_H = mean(Heterozygosity),
            sd_H = sd(Heterozygosity))

# A tibble: 4 × 6
#common          n median_H     IQR_H  mean_H     sd_H
#<chr>       <int>    <dbl>     <dbl>   <dbl>    <dbl>
#  1 GRPC           49  0.00382 0.000363  0.00376 0.000275
#2 LEPC          426  0.00341 0.000358  0.00355 0.000535
#3 STGR           29  0.00333 0.0000860 0.00331 0.000125
#4 STGR / GRPC     2  0.00408 0.000544  0.00408 0.000769


# Omnibus test, then pairwise tests on the three species only
kruskal.test(Heterozygosity ~ common, data = HET_FILT)
#Kruskal-Wallis chi-squared = 48.349, df = 3, p-value = 1.795e-10
sp3 <- subset(HET_FILT, common != "STGR / GRPC")
kruskal.test(Heterozygosity ~ common, data = sp3)
#Kruskal-Wallis chi-squared = 46.217, df = 2, p-value = 9.207e-11
pairwise.wilcox.test(sp3$Heterozygosity, sp3$common, p.adjust.method = "BH")
#     GRPC    LEPC  
#LEPC 8.2e-09 -     
#  STGR 3.4e-14 0.0047

#P value adjustment method: BH 

wilcox_effsize(sp3, Heterozygosity ~ common)   # r: ~0.1 small, ~0.3 medium, ~0.5 large
# A tibble: 3 × 7
#.y.            group1 group2 effsize    n1    n2 magnitude
#* <chr>          <chr>  <chr>    <dbl> <int> <int> <ord>    
#  1 Heterozygosity GRPC   LEPC     0.268    49   426 small    
#2 Heterozygosity GRPC   STGR     0.757    49    29 large    
#3 Heterozygosity LEPC   STGR     0.132   426    29 small  

