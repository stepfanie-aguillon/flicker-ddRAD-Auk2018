# This script contains the R code for analyses of FST and heterozygosity used in this manuscript:
#
# Authors: Aguillon SM, Campagna L, Harrison RG, Lovette IJ
# Year: 2018
# Title: A flicker of hope: Genomic data distinguish Northern Flicker taxa despite low levels of divergence
# Journal Info: The Auk, 135(3), 748-766
# DOI: 10.1642/AUK-18-7.1
#
# Please cite the paper if you use these scripts
#

# load packages
library(adegenet)
library(hierfstat)
library(ggplot2)

# set directory where STRUCTURE file is stored
setwd("~/Desktop/ddRAD2015/")

# load data in STRUCTURE format using adegenet
genofile <- read.structure("noMAF_denovo_structure_formatted.stru",onerowperind=FALSE,n.ind=38,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

# transform file into format for hierfstat
genofile_hier <- genind2hierfstat(genofile)



####### DIFFERENTIATION BETWEEN TAXA #######

####### PAIRWISE FST #######

# pairwise fst between the three flickers
pop_Fst <- pairwise.neifst(genofile_hier,diploid=TRUE)
#1      2      3
#1     NA 0.0181 0.1160
#2 0.0181     NA 0.1386
#3 0.1160 0.1386     NA




####### PER SNP FST #######

## red-shafted flicker versus yellow-shafted flicker ##

# load structure file with only RSFL and YSFL
genofile_RSFL_YSFL <- read.structure("noMAF_denovo_structure_formatted_RSFL_YSFL.stru",onerowperind=FALSE,n.ind=33,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

# run basic.stats in hierfstat
basic_stats_RSFL_YSFL <- basic.stats(genofile_RSFL_YSFL)
# write per SNP fst output to file
write.table(basic_stats_RSFL_YSFL$perloc,"PerLocus_summary_RSFL_YSFL_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)


## red-shafted flicker versus gilded flicker ##

# load structure file with only RSFL and GIFL
genofile_RSFL_GIFL <- read.structure("noMAF_denovo_structure_formatted_RSFL_GIFL.stru",onerowperind=FALSE,n.ind=18,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

# run basic.stats in hierfstat
basic_stats_RSFL_GIFL <- basic.stats(genofile_RSFL_GIFL)
# write per SNP fst output to file
write.table(basic_stats_RSFL_GIFL$perloc,"PerLocus_summary_RSFL_GIFL_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)


## yellow-shafted flicker versus gilded flicker ##

# load structure file with only YSFL and GIFL
genofile_YSFL_GIFL <- read.structure("noMAF_denovo_structure_formatted_YSFL_GIFL.stru",onerowperind=FALSE,n.ind=25,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

# run basic.stats in hierfstat
basic_stats_YSFL_GIFL <- basic.stats(genofile_YSFL_GIFL)
# write per SNP fst output to file
write.table(basic_stats_YSFL_GIFL$perloc,"PerLocus_summary_RSFL_GIFL_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)



##### PLOTTING FST FIGURE #####

# combine three output files into format that is readable by ggplot2 as:
# FSTP Comparison
# where FSTP is the fst of a given SNP and Comparison indates which two taxa are being compared

# read in file, omit NAs
FST <- read.table("./PerLocus_FSTP_combined_revisions_ggplot.txt",sep="",header=TRUE)
FST <- na.omit(FST)

# scale FST to 0 if a SNP has a negative FST
FST$FSTP_scaled <- ifelse(FST$FSTP<0,0,FST$FSTP)

# factor the comparison column for plotting
FST$Comparison2 <- factor(FST$Comparison,c("RY","RG","YG"))

# plot of fst
# Saved as PDFs with the dimensions: 5" (H) x 5" (L)
fst_plot <- ggplot() +
  geom_jitter(data=FST,aes(x=Comparison2,y=FSTP_scaled),alpha=0.2,color="black") +
  geom_hline(yintercept=0.0181,color="red") + #edit these in illustrator
  geom_hline(yintercept=0.1160,color="green") + #edit these in illustrator
  geom_hline(yintercept=0.1386,color="blue") + #edit these in illustrator
  labs(x="Comparison",y="Per SNP Fst") +
  theme_classic() +
  theme(axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10,color="black"))
fst_plot






####### HETEROZYGOSITY #######
# using adegenet

# separate genind object by population
genofile_sep <- seppop(genofile,res.type="genind",drop=FALSE)


# heterozygosity in red-shafted flicker
RSFL_het <- summary(genofile_sep$`1`)
head(RSFL_het$Hobs)
# write output
write.table(RSFL_het$Hobs,"RSFL_het_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)
# can also calculate expected heterozygosity if needed == Hexp

# heterozygosity in yellow-shafted flicker
YSFL_het <- summary(genofile_sep$`2`)
head(YSFL_het$Hobs)
# write output
write.table(YSFL_het$Hobs,"YSFL_het_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)

# heterozygosity in gilded flicker
GIFL_het <- summary(genofile_sep$`3`)
head(GIFL_het$Hobs)
# write output
write.table(GIFL_het$Hobs,"GIFL_het_noMAF.txt",sep="\t",quote=FALSE,row.names=TRUE)


##### PLOTTING HETEROZYGOSITY FIGURE #####

# combine three output files into format that is readable by ggplot2 as:
# het taxa
# where het is the heterozygosity of a given SNP and taxa is either RSFL, YSFL, or GIFL

# read in combined file
het <- read.table("./PerLocus_het_combined_revised_ggplot.txt",sep="",header=TRUE)
het$taxa2 <- factor(het$taxa,c("GIFL","RSFL","YSFL"))

# summaize results
library(dplyr)
het %>% group_by(taxa) %>% summarize(avg=mean(het))
#1   GIFL 0.09411218
#2   RSFL 0.12132004
#3   YSFL 0.11654063


# boxplot of heterozygosity
# Saved as PDFs with the dimensions: 5" (H) x 5" (L)
fig_colors <- c("#A020F0","#F21924","#FFFF00")
het_plot <- ggplot() +
  geom_boxplot(data=het,aes(x=taxa2,y=het,fill=taxa2)) +
  scale_fill_manual(values=fig_colors) +
  labs(x="Taxa",y="Heterozygosity") +
  theme_classic() +
  theme(legend.position="none",axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_blank(),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10,color="black"))
het_plot
