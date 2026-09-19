#Load library
library(reshape2)
library(ggplot2)

#read in metadata
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_nexus_new_plus_shotguns/heterozygosity_extended_nexus.xlsx")

#Extract relevant information
sub<-metadata[,c("ID","GROUP","SPECIES","fROH_100kb-1Mb","fROH_1Mb","fROH_total")]

#Convert to long format
melt_data<-melt(sub,id.vars = c("ID","GROUP","SPECIES"))
melt_data$GROUP <- factor(melt_data$GROUP, levels = c("Allopatric","Sympatric"))

#Plot ROH by groups
ggplot(melt_data, aes(fill=SPECIES, y=value, x=reorder(ID,value))) +
    geom_bar(stat="identity") +
    facet_grid(variable~GROUP, scales="free") +
    theme_classic() +
    theme(
        panel.border = element_rect(color="black", fill=NA, linewidth=1),
        strip.background = element_rect(color="black", fill="grey90", linewidth=1),
        axis.line = element_line(color="black", linewidth=1)
    ) +scale_fill_manual("", values=c("Tympanuchus cupido" = "goldenrod","Tympanuchus pallidicinctus"="brown","Tympanuchus phasianellus" = "black","Tympanuchus phasianellus/Tympanuchus cupido" = "grey")) +
    ylab("fROH") + xlab("Sample (N=506)") +
    theme(legend.position="bottom") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=14, face="italic")) +
    theme(strip.text = element_text(size=16))+theme(panel.spacing = unit(0, "lines"))

#Test for normality
shapiro.test(sub$`fROH_100kb-1Mb`)


shapiro.test(sub$`fROH_1Mb`)


shapiro.test(sub$`fROH_total`)


#Pairwise test
pairwise.wilcox.test(sub$`fROH_100kb-1Mb`, sub$GRP, p.adjust.method = "BH")

