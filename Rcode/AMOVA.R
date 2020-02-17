# This script contains the R code to conduct the AMOVA used in this manuscript:
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
library(poppr)
library(ggplot2)

# set directory where STRUCTURE file is stored
setwd("~/Desktop/ddRAD2015/")

# load data in STRUCTURE format using adegenet
genofile <- read.structure("noMAF_denovo_structure_formatted.stru",onerowperind=FALSE,n.ind=38,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

##### AMOVA #####

poppr(genofile)
#Pop  N MLG eMLG       SE    H  G lambda E.5  Hexp     Ia   rbarD     File
#1     1 13  13   10 7.30e-08 2.56 13  0.923   1 0.137 491.86 0.05282 genofile
#2     2 20  20   10 0.00e+00 3.00 20  0.950   1 0.135 580.22 0.05268 genofile
#3     3  5   5    5 0.00e+00 1.61  5  0.800   1 0.111   8.84 0.00187 genofile
#4 Total 38  38   10 1.72e-06 3.64 38  0.974   1 0.138 477.18 0.03483 genofile

#locus_table(genofile)
#locus_table(genofile,pop = "1")

# transforms genind into genclone
strata(genofile) <- data.frame(genofile$pop)
genofile <- as.genclone(genofile)
genofile
#This is a genclone object
#-------------------------
#  Genotype information:
#
#  38 original multilocus genotypes
#38 diploid individuals
#16670 codominant loci
#
#Population information:
#
#  1 stratum - genofile.pop
#3 populations defined - 1, 2, 3

# outputs table of sample numbers
table(strata(genofile,~genofile.pop))
#3  1  2
#5 13 20


# runs AMOVA
amova_result <- poppr.amova(genofile,~genofile.pop,within=TRUE,cutoff=1)
amova_result

# permutation test of AMOVA results
amova_test <- randtest(amova_result,nrepet=999)
amova_test
plot(amova_test)

# AMOVA output looks like
#$componentsofcovariance
#                                                    Sigma          %
#Variations  Between genofile.pop                 62.20375   6.447319
#Variations  Between samples Within genofile.pop  60.45268   6.265823
#Variations  Within samples                      842.14394  87.286859
#Total variations                                964.80038 100.000000



##### PLOTTING FIGURES #####

# save AMOVA results as a separate file with the format (as below) to plot with ggplot2

###file AMOVA_ggplot.txt
#perc	statistic	group
#6.447319	FST	1
#6.265823	FIS	1
#87.286859	FIT	1

# read in ggplot formatted results file
amova <- read.table("./AMOVA_ggplot.txt",sep="",header=TRUE)
amova$statistic2 <- factor(amova$statistic,c("FST","FIS","FIT"))

# gray-scale barplot of AMOVA results
amova_plot <- ggplot() +
  geom_col(data=amova,aes(x=group,y=perc,fill=statistic2)) +
  labs(x="AMOVA",y="Percent Variation") +
  theme_classic() +
  scale_fill_brewer(palette = 6) +
  theme(axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10,color="black"))
amova_plot
