setwd("/Volumes/jaeyoonc/_THESIS/JCthesis")

# load packages
library(edgeR)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(tidyverse) # change from chr1 to 1; write_lines
library(qqman) # manhattan
# genome wide annotation for mouse
library(org.Mm.eg.db) 

# read in txt file of samples key
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
  # combined CRS
  "~/Desktop/THESIS/covs/combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr17_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr18_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr19_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  # no rr20
  "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")

# yall is all 19 cov files, read and collate counts
yall19 <- readBismark2DGE(files19, sample.names=Sample19)
# set group name = population # need rep each 2x!!!
yall19$samples$group <- factor(rep(targets19$Population, each=2)) 
# also add "chr" to chr numbers
row.names(yall19$counts) <- paste0("chr",row.names(yall19$counts))
yall19$genes$Chr <- paste0("chr",yall19$genes$Chr)
nrow(yall19) # 10027860 # vs nrow(yall) (with rr20) was # 10152084

# filter MT and Y chr
keep19 <- rep(TRUE, nrow(yall19))
Chr19 <- as.character(yall19$genes$Chr)
sort(unique(Chr19)) # there are no unassembled chr
# remove Y chr and mito DNA
keep19[Chr19=="chrY"] <- FALSE
keep19[Chr19=="chrMT"] <- FALSE
table(keep19) # removing 20971 False
yall19 <- yall19[keep19,, keep.lib.sizes=FALSE]
nrow(yall19) # 10,006,889 (from originally 10027860, keeping 99.7%)
# vs nrow(yall) with 20 samples (including rr20) was 10130650 

# sort the DGEList into order, chr1 to X
ChrNames19 <- paste0("chr",c(1:19, "X"))
yall19$genes$Chr <- factor(yall19$genes$Chr, levels=ChrNames19)
o19 <- order(yall19$genes$Chr, yall19$genes$Locus)
yall19 <- yall19[o19,]

################################################################

# gene annotation
TSS19 <- nearestTSS(yall19$genes$Chr, yall19$genes$Locus, species="Mm") # mouse
# this function is using org.Mm.eg.db
# add gene annotation to yall$genes
yall19$genes$EntrezID <- TSS19$gene_id
yall19$genes$Symbol <- TSS19$symbol
yall19$genes$Strand <- TSS19$strand
yall19$genes$Distance <- TSS19$distance
yall19$genes$Width <- TSS19$width
head(yall19$genes)

### up until here was the same as for DMCs

################# now define promoter ################# 

# promoter as 2kb upstream to 1kb downstream of TSS, Lisa ok'ed
InPromoter19 <- yall19$genes$Distance >= -1000 & yall19$genes$Distance <= 2000
# compare what happens: starting with 10027860 (vs 10130650 with rr20) CpG sites
sum(InPromoter19) # left with 2,283,385 (vs 2,294,885 with rr20) CpG sites meeting threshold of promoter
# 2283385/10,027,860*100 # thats ~22.8% left ok.
# filter
yIP19 <- yall19[InPromoter19,,keep.lib.sizes=FALSE]
# how many genes?
length(unique(yall19$genes$Symbol)) # 25889 (vs 25894 with rr20) genes
length(unique(yIP19$genes$Symbol)) # 24632 (vs 24661 with rr20) unique genes contained in the 2283614 promoters

# now compute total counts for *each promoter*
# each row corresponds to not a CpG site but one row per gene!
ypr19 <- rowsum(yIP19, yIP19$genes$Symbol, reorder=FALSE)
# ypr19$genes$Symbol <- NULL # redundant now since ypr$genes rowname is gene symbols

### filter
# save Me and Un for ease of access (odd columns are Me, even Un)
Methyl.pr19 <- gl(2,1,ncol(ypr19), labels=c("Me","Un"))
# filter out CpGs that are never or always methylated bc no info about differential methyl
Me.pr19 <- ypr19$counts[,Methyl.pr19=="Me"] # extract # Me
Un.pr19 <- ypr19$counts[,Methyl.pr19=="Un"] # extract # Un
HasBoth.pr19 <- rowSums(Me.pr19) > 0 & rowSums(Un.pr19) > 0
table(HasBoth.pr19) # FALSE (always methyl or always unmethyl) = 2029; TRUE = 22603 

# COVERAGE: set threshold of 5 occurences in all 19 samples
# coverage: returns each position and coverage for each sample
Coverage.pr19 <- ypr19$counts[, Methyl.pr19=="Me"] + ypr19$counts[, Methyl.pr19 =="Un"]

# ncol(Coverage.pr) is 20
keep5.pr19 <- rowSums(Coverage.pr19 >= 5) == 19
table(keep5.pr19) # FALSE (too low coverage) 7469; TRUE 17163 genes

# apply both filters
ypr5.19 <- ypr19[keep5.pr19 & HasBoth.pr19,, keep.lib.sizes=FALSE]

length(unique(ypr5.19$genes$EntrezID)) # 17155 (vs 17152) unique genes!
nrow(ypr5.19$counts) # confirm 17155
sum(unique(rownames(ypr5.19$genes)) == "Tert") # includes Tert
sum(unique(rownames(ypr5.19$genes)) == "Bdnf") # includes Bdnf

# set library size to be equal for each pair (Me and Un) of libraries
# this is bc each pair is treated as a unit in analysis
# add up for each sample lib size for methylated with unmethylated
TotalLibSize.pr19 <- ypr5.19$samples$lib.size[Methyl.pr19=="Me"] + 
  ypr5.19$samples$lib.size[Methyl.pr19=="Un"]
# example: before we do this, ** fyi these values outdated but method is same
# head(y$samples):
# rr01-Me lib.size is 1113815
# rr01-Un lib.size is 7079083
# after we do this (command y5$samples$lib.size <- rep(TotalLibSize, each=2))
# both rr01-Me and rr01-Un are 8192898, the sum

# exploring this more as global methylation
View(as.data.frame(ypr5.19$samples)) ### export to excel? ###

# apply library size
ypr5.19$samples$lib.size <- rep(TotalLibSize.pr19, each=2)

############ MDS PLOT #############
# M-value
Me.pr19 <- ypr5.19$counts[,Methyl.pr19=="Me"] # extract # Me
Un.pr19 <- ypr5.19$counts[,Methyl.pr19=="Un"] # extract # Un
Mpr19 <- log2(Me.pr19+2) - log2(Un.pr19+2) # calculate M-value, add 2 for zero frequency case
colnames(Mpr19) <- Sample19

# generate MDS plot
# set colors for visualisation
colors <- c(rep("blue",4),
            rep("black",5),
            rep("grey",5),
            rep("red",5)) # 5 CRSc
# plot with color legend (note: takes a long time)
pdf("2021.02.08_MDS_DMGs19promoters.pdf")
plotMDS(Mpr19, col=colors, gene.selection="pairwise", main="Promoters") + legend(
  "topleft",
  bty = "n",
  c("CORT", "VEHICLE", "CONTROL", "CRS"),
  fill = c("blue", "black", "grey", "red")
)
dev.off()

########
# design matrix *** is same as for edgeR all ***
designSL19 <- model.matrix(~0+Population, data=targets19)
# # colnames(designSL) # "PopulationCONTROL" "PopulationCORT" "PopulationCRS" "PopulationVEHICLE"
colnames(designSL19) <- c("CONTROL", "CORT", "CRS", "VEHICLE")
design19 <- modelMatrixMeth(designSL19)
#  each DNA sample generates two counts, a count of methylated reads and a count of unmethylated reads, for each genomic locus for each sample. The function converts sample-level information about the treatment conditions to make an appropriate design matrix with two rows for each sample. Counts are assumed to be ordered as methylated and then unmethylated by sample.
# nrow(design) # 38
# ncol(design) # 23: for 19 samples + 4 treatments

# dispersion
yprdisp19 <- estimateDisp(ypr5.19, design=design19, trend="none") # takes long long
yprdisp19$common.dispersion # 0.4373162 (vs 0.4515481) very high
summary(yprdisp19$prior.df) # 4.478 (vs 4.707) for all 
# ********* dispersion? ********** ? #
# Inf means all CpG-wise dispersions are equal to the common dispersion -- but not Inf

# visualise dispersion estimates with BCV plot
plotBCV(yprdisp, main = "Dispersion for promoter CpGs")
# compare to: plotBCV(ydisp, main = "Dispersion for all CpGs")

################ differential methylation analysis ################ 
fitpr19 <- glmFit(yprdisp19, design19)
# *** contrasts are same as before ***
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

# show number of differentially methylated CpG sites at FDR of 5%
summary(decideTests(lrt.CORTpr19, method = "fdr")) # 0 NONE! # also 0 with bonferroni # 0 sig, 17155 not sig

# from before doing promoters: Sec13, Nat14, Anks6
### note: Nat14 is the only one with distance -810
### Nat14 is one of the top hypomethylated CpGs found in one study, Littlejohn et al.

# for CRS
summary(decideTests(lrt.CRSpr19, method = "fdr")) # 1970 down, 119 up
# with rr20, had been 1568 down and 76 up

# take just the 2089 diff methylated
DMC_CRSpr19 <- topTags(lrt.CRSpr19, n=2300)[1:2089,]
# separate for hypo & hypermethylated
hypo_CRSpr19 <- DMC_CRSpr19[DMC_CRSpr19$table$logFC < 0, ] # hypomethyl nrow(hypo_CRSpr19) # 1970
hyper_CRSpr19 <- DMC_CRSpr19[DMC_CRSpr19$table$logFC > 0, ] # hypermethyl nrow(hyper_CRSpr19) # 119

# export as dataframe including the symbol, location, and 
all_DMCs <- as.data.frame(DMC_CRSpr19)
write.table(all_DMCs, "/Volumes/jaeyoonc/_THESIS/JCthesis/allCRS_2089DMCs_1970Hypo_119Hyper.txt")



# list of the 2089 gene symbols
hypo_CRSpr19_gene <- unique(rownames(hypo_CRSpr19)) # length(hypo_CRSpr19_gene) # 1970
hypo_CRSpr19_ENTREZID <- unique(hypo_CRSpr19$table$EntrezID) # length(hypo_CRSpr19_geneENTREZID) # 1970

hyper_CRSpr19_gene <- unique(rownames(hyper_CRSpr19)) # 119
hyper_CRSpr19_ENTREZID <- unique(hyper_CRSpr19$table$EntrezID) # length(hyper_CRSpr19_ENTREZID) # 119


# save gene names as files
write_lines(hypo_CRSpr19_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/1970hypoCRS19genesPROMOTER19.txt")
write_lines(hypo_CRSpr19_ENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/1970hypoCRS19genesENTREZID.txt")

write_lines(hyper_CRSpr19_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/119hyperCRS19genesPROMOTER19.txt")
write_lines(hyper_CRSpr19_ENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/119hyperCRS19genesENTREZID.txt")


# the original 17155 genes
allgenes <- rownames(ypr5.19$genes)
allgenesENTREZID <- unique(ypr5.19$genes$EntrezID) # length(allgenesENTREZID) 17155


length(allgenes) # 17155
write_lines(allgenes, "/Volumes/jaeyoonc/_THESIS/JCthesis/all17155genes_noRR20.txt")
write_lines(allgenesENTREZID, "/Volumes/jaeyoonc/_THESIS/JCthesis/all17155genes19_ENTREZID.txt")


# # lets check if there is any overlap between here and analysis with all CpG sites
# # for hypomethylated genes:
# sum(hypo_CRS_gene %in% hypo_CRSpr_gene) # 31 overlapping
# # for hypermethylated genes:
# sum(hyper_CRS_gene %in% hyper_CRSpr_gene) # 1 overlapping
# which(hyper_CRS_gene %in% hyper_CRSpr_gene)
# hyper_CRS_gene[111] # Pin1rt1

# visualise

tiff("2021.02.08.CORTpr_logFC.tiff",
     width=550, height=500)
plotMD(lrt.CORTpr19, main="CORT vs VEHICLE methylation in promoters", 
       ylim=c(-10,10))
dev.off()

tiff("2021.02.08.CRSpr19_logFC.tiff",
     width=550, height=500)
plotMD(lrt.CRSpr19, main="CRS vs CONTROL methylation in promoters",
       ylim=c(-10,10))
dev.off()


#####################################
#####################################
#####################################

# by chromosome
ChrIndices <- list()
for (i in ChrNames19) ChrIndices[[i]] <- which(yprdisp19$genes$Chr==i)
bychrCORTpr19 <- fry(yprdisp19, index=ChrIndices, design=design19, contrast=contr.CORT19)
bychrCRSpr19 <- fry(yprdisp19, index=ChrIndices, design=design19, contrast=contr.CRS19)

# --> these to excel ByChromosome

###################################################
###################################################

# make heatmap -- based on M values
# exporting the files here to make heatmap in separate R file

Mpr.vals <- data.frame(Mpr19)
# check names are same
nrow(Mpr.vals) - sum(rownames(Mpr.vals) == rownames(ypr5.19$genes)) # 0

# now we want to only select for CRS or CORT
# CRS
rrCRSpr <- c("rr11", "rr12", "rr13", "rr14", "rr15", "rr16c", "rr17c", "rr18c", "rr19c", "rr21c") # no rr20!
MCRSpr <- Mpr.vals[,rrCRSpr]
# transpose for heatmap
MCRSpr.map <- t(MCRSpr)
ncol(MCRSpr.map) # 17155 genes with CpGs in promoter
nrow(MCRSpr.map) # 10 conditions (rr's)
write.table(MCRSpr.map, "MvaluesCRSpr19_noRR20_forHeatmap.txt")

# CORT -- not expecting change. but bc of filtering, now have 3 more genes than before
# now 17155 genes without rr20 vs previously 17152 genes with rr20
rrCORTpr <- c("rr01","rr02","rr03","rr04","rr06","rr07","rr08","rr09","rr10")
MCORTpr <- Mpr.vals[,rrCORTpr]
# transpose for heatmap
MCORTpr.map <- t(MCORTpr)
ncol(MCORTpr.map) # 17155 genes with CpGs in promoter
nrow(MCORTpr.map) # 9 conditions (rr's)
write.table(MCORTpr.map, "MvaluesCORTpr19_noRR20_forHeatmap.txt")

# what about for all 19 samples?
MAll19pr <- Mpr.vals
# transpose for heatmap
MAll19pr <- t(MAll19pr)
ncol(MAll19pr) # 17155 genes with CpGs in promoter
nrow(MAll19pr) # 19 conditions (rr's)
write.table(MAll19pr, "Mvalues19pr_noRR20_forHeatmap.txt")

