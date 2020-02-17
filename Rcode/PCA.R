# This script contains the R code to make PCA plots used in this manuscript:
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
library(SNPRelate)
library(ggplot2)


# set directory where VCF files are stored
setwd("~/Desktop/ddRAD2015/")

# load SNP data in VCF format
snpgdsVCF2GDS(vcf.fn="./noMAF_denovo.vcf", out.fn="all_snps_denovo.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)

# summarize input file
snpgdsSummary("./all_snps_denovo.gds")

# open file
genofile <- snpgdsOpen("./all_snps_denovo.gds")
read.gdsn(index.gdsn(genofile,"sample.id"))

# get missing data rate across samples
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
# write to file
write.table(miss,"denovo_noMAF_missing.txt",sep="\t",quote=FALSE,row.names=TRUE)


##### PCA #####

# run PCA using SNPRelate
pca <- snpgdsPCA(gdsobj = genofile,autosome.only=FALSE)

# get percent variation explained for each PC axis
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# pull sample ID + first four PC axes
pca_coords <- data.frame(ID = pca$sample.id,
                         pc1 = pca$eigenvect[,1],    # the first eigenvector
                         pc2 = pca$eigenvect[,2],    # the second eigenvector
                         pc3 = pca$eigenvect[,3],
                         pc4 = pca$eigenvect[,4],
                         stringsAsFactors = FALSE)
head(pca_coords)

# add info on species and locality labels
sample_info <- read.table("./geographic_codes.txt",sep="",header=TRUE)

# merge pca results with sample info by ID number
pca_coords_merged <- merge(pca_coords,sample_info,by.x="ID")
# keep only desired columns
pca_coords_merged <- select(pca_coords_merged,ID,Species,Pop,pc1,pc2,pc3,pc4)



##### PLOTTING FIGURES #####
# Saved as PDFs with the dimensions: 3.5" (H) x 4" (L)

fig_colors <- c("#A020F0","#F21924","#FFFF00")

# scatterplot of PC1 versus PC2
pca_scatter1_2 <- ggplot() +
  geom_hline(aes(yintercept=0),color="gray") +
  geom_vline(aes(xintercept=0),color="gray") +
  geom_point(data=pca_coords_merged,aes(x=pc1,y=pc2,fill=Species),size=4,alpha=0.75,shape=21,stroke=0.2) +
  labs(x = "PC1 (5.62%)", y = "PC2 (3.57%)") +
  scale_fill_manual(values=fig_colors) +
  theme_classic() +
  theme(legend.position="none",axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10))
pca_scatter1_2


# scatterplot of PC2 versus PC3
pca_scatter2_3 <- ggplot() +
  geom_hline(aes(yintercept=0),color="gray") +
  geom_vline(aes(xintercept=0),color="gray") +
  geom_point(data=pca_coords_merged,aes(x=pc2,y=pc3,fill=Species),size=4,alpha=0.75,shape=21,stroke=0.2) +
  scale_fill_manual(values=fig_colors) +
  labs(x = "PC2 (3.57%)", y = "PC3 (3.35%)") +
  theme_classic() +
  theme(legend.position="none",axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10))
pca_scatter2_3


# plot of percent variation explained by each PC axis
PC_axis <- cbind(1:35)
plot_df <- as.data.frame(cbind(pc.percent,PC_axis))

var_plot <- ggplot() +
  geom_point(data=plot_df,aes(x=PC_axis,y=pc.percent),color="black") +
  labs(x = "Eigenvector", y = "Percent Variation Explained") +
  theme_classic() +
  theme(legend.position="none",axis.line.x=element_line(color="black"),axis.line.y=element_line(color="black"),  axis.title.x=element_text(face="bold",size=12),axis.title.y=element_text(face="bold",size=12),axis.text=element_text(size=10))
var_plot
