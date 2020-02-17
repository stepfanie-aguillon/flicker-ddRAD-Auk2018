
# This script contains the R code for BayeScan used in this manuscript:
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

# load data in STRUCTURE format
genofile <- read.structure("noMAF_denovo_structure_formatted.stru",onerowperind=FALSE,n.ind=38,n.loc=16670,col.lab=1,col.pop=2,row.marknames=1,NA.char="0",ask=FALSE)

# transform data into hierfstat format
genofile_hier <- genind2hierfstat(genofile)

# output file in the format for BayeScan
write.bayescan(genofile_hier,diploid=TRUE,fn="denovo_noMAF.bsc")



####################################
##### RUN BAYESCAN IN TERMINAL #####
####################################



##### PROCESS RESULTS IN R #####

# read file with per SNP FST
# obtained from basic.stats(genofile) in adegenet
fst_summary <- read.table("PerLocus_summary_noMAF.txt",sep="",header=TRUE)
# file form:
#Num	Locus	Ho	Hs	Ht	Dst	Htp	Dstp	Fst	Fstp	Fis	Dest
#1	1_49	0.0185	0.0191	0.0187	-4.00E-04	0.0185	-6.00E-04	-0.0213	-0.0324	0.0313	-6.00E-04
#2	1_97	0.0185	0.0191	0.0187	-4.00E-04	0.0185	-6.00E-04	-0.0213	-0.0324	0.0313	-6.00E-04
#3	7_15	0.0185	0.0192	0.0187	-4.00E-04	0.0185	-7.00E-04	-0.0235	-0.0357	0.0345	-7.00E-04
#...


# process results with false discovery rate of 0.01
bayescan_results_FDR0.01 <- plot_bayescan("denovo_noMAF_fst.txt",1,FDR=0.01,add_text=FALSE)

#list outliers
bayescan_results_FDR0.01$outliers
#[1]   157   544   870  1654  3195  3847  3849  4066  4068  4407  4629  5079  5087
#[14]  5299  5301  5304  7447  9428  9600  9601 10160 10672 11440 12118 12119 12120
#[27] 13149 13761 13768 14263 14416 14418 14421 14776 14908 15099 15549 15913 15914
#[40] 15941 16102 16103 16104 16105 16388 16390

# list number of outliers
bayescan_results_FDR0.01$nb_outliers
# 46

# combine locus name and BayeScan results
fst_summary_small <- fst_summary[,1:2]
loci_FDR0.01 <- cbind(bayescan_results_FDR0.01$outliers)
loci_nums_FDR0.01 <- merge(loci_FDR0.01,fst_summary_small,by.x="V1",by.y="Num")
colnames(loci_nums_FDR0.01) <- c("Num","Locus")

# list outliers by locus name
loci_nums_FDR0.01$Locus
#[1] 590_21    996_76    1361_103  2207_89   3542_73   4044_127  4213_118  4213_131
#[9] 4426_22   4583_103  4934_8    4934_60   5078_51   5078_57   5078_138  6521_123
#[17] 7982_98   7982_133  8325_87   8612_89   9134_52   9832_20   9832_48   9832_56
#[25] 11150_30  11958_28  12946_41  12946_77  12946_115 13181_108 13298_92  13409_128
#[33] 13768_88  14120_85  14120_109 14142_14  14306_65  14306_85  14306_96  14306_132
#[41] 15393_126 15393_130



##### RESULTS TABLE #####

# edit BayeScan results to name first column "Num"
# read in file as a table
BayeScan_results <- read.table("denovo_noMAF_fst_rename.txt",header=TRUE,sep="")

# merge results file with the outlier list
outlier_merge <- merge(loci_nums_FDR0.01,BayeScan_results,by="Num")
# write results table
write.table(outlier_merge,"BayeScan_results_outliers.txt",sep="\t",quote=FALSE,row.names=FALSE)





##### PLOTTING FIGURE #####
# BayeScan figure
bayescan_plot <- ggplot() +
  geom_point(data=bayescan_results_FDR0.01,aes(x=log10(qval),y=fst),color="black") +
  labs(x = "log10(q value)", y = "FST") +
  theme_classic() +
  scale_x_reverse() +
  theme(legend.position="none",axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10))
bayescan_plot
