setwd("/Volumes/jaeyoonc/_THESIS/JCthesis")

# load packages
library(edgeR)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(tidyverse) # change from chr1 to 1; write_lines
# genome wide annotation for mouse
library(org.Mm.eg.db) 

# note: reading in of files is the same as for DMC

# read in txt file of samples key (same file as in 001_)
targets19 <- read.delim("targets19_noRR20.txt", row.names = "Sample", stringsAsFactors=FALSE)
# nrow(targets19) # 19 rows
# ncol(targets19) # 2 cols
Sample19 <- row.names(targets19) # no rr05, no rr20
# files is a list of all 31 sample file locations
files19 <- c(
  # 21 original batch
  "~/Desktop/THESIS/covs/G695_rr01_S15_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr02_S16_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr03_S17_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr04_S18_val_1_bismark_bt2_pe.bismark.cov",
  # remove rr05 outlier
  "~/Desktop/THESIS/covs/G695_rr06_S20_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr07_S21_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr08_S22_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr09_S23_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr10_S24_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr11_S25_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr12_S26_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr13_S27_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr14_S28_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr15_S29_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr17_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr18_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr19_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  # no rr20 outlier
  "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")

# yall is all 19 cov files, read and collate counts
yall19 <- readBismark2DGE(files19, sample.names=Sample19)
# set group name = population # need rep each 2x!!!
yall19$samples$group <- factor(rep(targets19$Population, each=2)) 
# also add "chr" to chr numbers
row.names(yall19$counts) <- paste0("chr",row.names(yall19$counts))
yall19$genes$Chr <- paste0("chr",yall19$genes$Chr)
nrow(yall19) # 10027860

# filter MT and Y chr
keep19 <- rep(TRUE, nrow(yall19))
Chr19 <- as.character(yall19$genes$Chr)
sort(unique(Chr19)) # there are no unassembled chr
# remove Y chr and mito DNA
keep19[Chr19=="chrY"] <- FALSE
keep19[Chr19=="chrMT"] <- FALSE
table(keep19) # removing 20971 False
yall19 <- yall19[keep19,, keep.lib.sizes=FALSE]
nrow(yall19) # 10,006,889 cytosines

# sort the DGEList into order, chr1 to X
ChrNames19 <- paste0("chr",c(1:19, "X"))
yall19$genes$Chr <- factor(yall19$genes$Chr, levels=ChrNames19)
o19 <- order(yall19$genes$Chr, yall19$genes$Locus)
yall19 <- yall19[o19,]

################################################################

# gene annotation -- comes from the edgeR package that includes mm10 genome
TSS19 <- nearestTSS(yall19$genes$Chr, yall19$genes$Locus, species="Mm") # mouse
# this function is using org.Mm.eg.db mm10
# add gene annotation to yall$genes
yall19$genes$EntrezID <- TSS19$gene_id
yall19$genes$Symbol <- TSS19$symbol
yall19$genes$Strand <- TSS19$strand
yall19$genes$Distance <- TSS19$distance
yall19$genes$Width <- TSS19$width
head(yall19$genes)

### up until here was the same as for DMCs

################# now define PCE ################# 

# PCE as 2kb upstream to 1kb downstream of TSS
InPCE19 <- yall19$genes$Distance >= -1000 & yall19$genes$Distance <= 2000
# filter
yIP19 <- yall19[InPCE19,,keep.lib.sizes=FALSE]
# how many genes? --> more filtering to come
length(unique(yIP19$genes$Symbol)) # 24632 

# now compute total counts for *each PCE*
# each row corresponds to not a CpG site but one row per gene (PCE)!
ypr19 <- rowsum(yIP19, yIP19$genes$Symbol, reorder=FALSE)

### filter
# save Me and Un for ease of access (odd columns are Me, even Un)
Methyl.pr19 <- gl(2,1,ncol(ypr19), labels=c("Me","Un"))
# filter out never or always methylated bc no info about differential methyl
Me.pr19 <- ypr19$counts[,Methyl.pr19=="Me"] # extract # Me
Un.pr19 <- ypr19$counts[,Methyl.pr19=="Un"] # extract # Un
HasBoth.pr19 <- rowSums(Me.pr19) > 0 & rowSums(Un.pr19) > 0
table(HasBoth.pr19) # FALSE (always methyl or always unmethyl) = 2029; TRUE = 22603 

# COVERAGE: set threshold of 5 occurences in all 19 samples (you can change this)
# coverage: returns each position and coverage for each sample
Coverage.pr19 <- ypr19$counts[, Methyl.pr19=="Me"] + ypr19$counts[, Methyl.pr19 =="Un"]
keep5.pr19 <- rowSums(Coverage.pr19 >= 5) == 19
table(keep5.pr19) # TRUE 17163 genes

# apply both filters
ypr5.19 <- ypr19[keep5.pr19 & HasBoth.pr19,, keep.lib.sizes=FALSE]

# total number of genes = 17155
length(unique(ypr5.19$genes$EntrezID)) # TOTAL 17155
nrow(ypr5.19$counts) # confirm 17155

# set library size to be equal for each pair (Me and Un) of libraries
# this is bc each pair is treated as a unit in analysis
# add up for each sample lib size for methylated with unmethylated
TotalLibSize.pr19 <- ypr5.19$samples$lib.size[Methyl.pr19=="Me"] + 
  ypr5.19$samples$lib.size[Methyl.pr19=="Un"]

# exploring this as global methylation
View(as.data.frame(ypr5.19$samples))

# apply library size
ypr5.19$samples$lib.size <- rep(TotalLibSize.pr19, each=2)

########
# design matrix *** is same steps as for DMC ***
designSL19 <- model.matrix(~0+Population, data=targets19)
colnames(designSL19) <- c("CONTROL", "CORT", "CRS", "VEHICLE")
design19 <- modelMatrixMeth(designSL19)
#  each DNA sample generates two counts, a count of methylated reads and a count of unmethylated reads, for each genomic locus for each sample. The function converts sample-level information about the treatment conditions to make an appropriate design matrix with two rows for each sample. Counts are assumed to be ordered as methylated and then unmethylated by sample.

# dispersion
yprdisp19 <- estimateDisp(ypr5.19, design=design19, trend="none")
yprdisp19$common.dispersion 
summary(yprdisp19$prior.df)

# visualise dispersion estimates with BCV plot
plotBCV(yprdisp, main = "Dispersion for PCE CpGs")

################ differential methylation analysis ################ 

fitpr19 <- glmFit(yprdisp19, design19)

# # comparison CORT
contr.CORT19 <- makeContrasts(CORT19=CORT-VEHICLE, levels=design19)
# # comparison CRS
contr.CRS19 <- makeContrasts(CRS19=CRS-CONTROL, levels=design19)

# test CORT
lrt.CORTpr19 <- glmLRT(fitpr19,contrast=contr.CORT19)
# test CRS
lrt.CRSpr19 <- glmLRT(fitpr19,contrast=contr.CRS19)

# view top differentially methylated CpG sites
topTags(lrt.CORTpr19)
topTags(lrt.CRSpr19)

# show number of differentially methylated CpG sites at bonferroni 5%
summary(decideTests(lrt.CORTpr19, method = "bonferroni")) # 0 NONE!

# for CRS
summary(decideTests(lrt.CRSpr19, method = "bonferroni")) # 1970 down, 119 up

# take just the 2089 diff methylated
DMC_CRSpr19 <- topTags(lrt.CRSpr19, n=2300)[1:2089,]
# separate for hypo & hypermethylated
hypo_CRSpr19 <- DMC_CRSpr19[DMC_CRSpr19$table$logFC < 0, ] # hypomethyl nrow(hypo_CRSpr19) # 1970
hyper_CRSpr19 <- DMC_CRSpr19[DMC_CRSpr19$table$logFC > 0, ] # hypermethyl nrow(hyper_CRSpr19) # 119

# export as dataframe including the symbol, location
all_DMCs <- as.data.frame(DMC_CRSpr19)
write.table(all_DMCs, "/Volumes/jaeyoonc/_THESIS/JCthesis/allCRS_2089DMCs_1970Hypo_119Hyper.txt")

# list of the 2089 gene symbols
hypo_CRSpr19_gene <- unique(rownames(hypo_CRSpr19)) # length(hypo_CRSpr19_gene) # 1970
hypo_CRSpr19_ENTREZID <- unique(hypo_CRSpr19$table$EntrezID) # length(hypo_CRSpr19_geneENTREZID) # 1970

hyper_CRSpr19_gene <- unique(rownames(hyper_CRSpr19)) # 119
hyper_CRSpr19_ENTREZID <- unique(hyper_CRSpr19$table$EntrezID) # length(hyper_CRSpr19_ENTREZID) # 119

# save gene names as files both as gene symbol and entrezID --> use for clusterProfiler
write_lines(hypo_CRSpr19_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/1970hypoCRS19genesPROMOTER19.txt")
write_lines(hypo_CRSpr19_ENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/1970hypoCRS19genesENTREZID.txt")
write_lines(hyper_CRSpr19_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/119hyperCRS19genesPROMOTER19.txt")
write_lines(hyper_CRSpr19_ENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/119hyperCRS19genesENTREZID.txt")

# save the original 17155 genes
allgenes <- rownames(ypr5.19$genes)
allgenesENTREZID <- unique(ypr5.19$genes$EntrezID) # length(allgenesENTREZID) 17155
length(allgenes) # 17155
write_lines(allgenes, "/Volumes/jaeyoonc/_THESIS/JCthesis/all17155genes_noRR20.txt")
write_lines(allgenesENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/all17155genes19_ENTREZID.txt")


# visualise MD plots -- THESE ARE FIGURE 13 IN THESIS
# CORT
pdf("CORTpce-logFC.pdf")
plotMD(lrt.CORTpr19, 
       main="CORT vs VEHICLE methylation in PCEs", 
       ylab="Change in methylation level (logFC)",
       ylim=c(-10,10),
       # col.main="#3E8C3D",
       cex.lab=1.3, # make fonts bigger
       cex.axis=1.1, 
       cex.main=1.7)
dev.off()

# CRS -- THESE ARE FIGURE 13 IN THESIS
pdf("CRSpce-logFC.pdf")
plotMD(lrt.CRSpr19, 
       main="CRS vs CONTROL methylation in PCEs",
       ylab="Change in methylation level (logFC)",
       ylim=c(-10,10),
       cex.lab=1.3, # make fonts bigger
       cex.axis=1.1, 
       cex.main=1.7)
dev.off()

#####################################

# by chromosome
ChrIndices <- list()
for (i in ChrNames19) ChrIndices[[i]] <- which(yprdisp19$genes$Chr==i)
# use fry rotational gene test
bychrCORTpr19 <- fry(yprdisp19, index=ChrIndices, design=design19, contrast=contr.CORT19)
bychrCRSpr19 <- fry(yprdisp19, index=ChrIndices, design=design19, contrast=contr.CRS19)

##################################
# hierarchical clustering

# load packages
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization

# export M values for all 17155 PCEs for each of the 19 samples
Mpr.vals <- data.frame(Mpr19)
write.txt(Mpr.vals, "/Volumes/jaeyoonc/_THESIS/Mvalues.txt")

# check names are same
nrow(Mpr.vals) - sum(rownames(Mpr.vals) == rownames(ypr5.19$genes)) # 0

# clustering -- THESE ARE FIGURE 12 IN THESIS
mCRSpr19 <- read.table("/Volumes/jaeyoonc/_THESIS/Mvalues.txt", header=TRUE)
nrow(mCRSpr19)
ncol(mCRSpr19)

# transpose to get into form for clustering
Mvals <- t(mCRSpr19)
Mvals <- scale(Mvals)

# calculate euclidean distance betweeen samples
d <- dist(Mvals, method = "euclidean")

# use agnes to find which clustering is best separating
hc2 <- agnes(Mvals, method = "complete")
# agglomerative coefficient
hc2$ac
# compare
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(Mvals, method = x)$ac
}
map_dbl(m, ac) # ward is best

# hierarchical clustering using ward
hc1 <- hclust(d, method = "ward" )
plot(hc1, cex = 0.6, hang= -1,
     main ="Hierarchical clustering by M-value",
     cex.main=0.9)
rect.hclust(hc1, k = 2, border = c("#8239B5", "#3E8C3D"))

# K means clustering with 2 groups
sub_grp <- cutree(hc1, k = 2)
fviz_cluster(list(data = Mvals, cluster = sub_grp, colour = c("#8239B5", "#3E8C3D")),
             labelsize = 8) + 
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  ggtitle("K-means clustering by M-value") +
  scale_colour_manual(values = c("#3E8C3D", "#8239B5")) +
  scale_fill_manual(values = c("#3E8C3D", "#8239B5")) 

