setwd("/Volumes/jaeyoonc/_THESIS/JCthesis/")

library(edgeR)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
# genome wide annotation for mouse
library(org.Mm.eg.db) 
library(tidyverse) #change from chr1 to 1; vice versa etc
library(qqman) # manhattan

# read in txt file of samples key
targets19 <- read.delim("targets19_noRR20.txt", row.names = "Sample", stringsAsFactors=FALSE)
# nrow(targets19) # 19 rows
# ncol(targets19) # 2 cols
Sample19 <- row.names(targets19) # no rr05, no rr20
# files is a list of all 19 sample file locations
files19 <- c(
  # original batch except rr16-19 & rr21 CRSc
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
nrow(yall19) # 10,027,860
# vs nrow(yall) (with rr20) was # 10152084

# filter MT and Y chr
keep19 <- rep(TRUE, nrow(yall19))
Chr19 <- as.character(yall19$genes$Chr)
sort(unique(Chr19)) # there are no unassembled chr
# remove Y chr and mito DNA
keep19[Chr19=="chrY"] <- FALSE
keep19[Chr19=="chrMT"] <- FALSE
table(keep19) # removing 21434 False
yall19 <- yall19[keep19,, keep.lib.sizes=FALSE]
nrow(yall19) # 10006889 (from originally 10027860, keeping 99.7%)
# vs nrow(yall) here was 10130650 

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

# store methylation cols for ease later
Methylation19 <- gl(2,1,ncol(yall19), labels=c("Me","Un")) # length 38

# coverage: returns each position and coverage for each sample
Coverage19 <- yall19$counts[, Methylation19=="Me"] + yall19$counts[, Methylation19 =="Un"]
# class(Coverage19)

### FILTER COVERAGE: set threshold of 5 occurrences in all 19 samples
# ncol(Coverage19) is 19
# keep8.19 <- rowSums(Coverage19 >= 8) == 19 # sum(keep8.19) --> 225,717 loci
        # vs keep8 with 20 samples was 214,314 loci
keep5.19 <- rowSums(Coverage19 >= 5) == 19 # sum(keep5.19) --> 433,721 loci
        # vs keep5 with 20 samples was 422202 loci

####### interesting that with 19 samples the coverage is higher ######

# filter out CpGs that are never or always methylated bc no info about differential methyl
Me19 <- yall19$counts[,Methylation19=="Me"] # extract # Me
Un19 <- yall19$counts[,Methylation19=="Un"] # extract # Un
HasBoth19 <- rowSums(Me19) > 0 & rowSums(Un19) > 0

# apply both filters
y5.19 <- yall19[keep5.19 & HasBoth19,, keep.lib.sizes=FALSE]

# all the genes captured now with this coverage threshold
length(unique(y5.19$genes$Symbol)) # 17778 unique genes vs y5: 17577 unique genes
nrow(y5.19$genes) # 417,036 loci
sum(unique(y5.19$genes$Symbol) == "Tert") # includes Tert
sum(unique(y5.19$genes$Symbol) == "Bdnf") # includes Bdnf

# set library size to be equal for each pair (Me and Un) of libraries
# this is bc each pair is treated as a unit in analysis
# add up for each sample lib size for methylated with unmethylated
TotalLibSize19 <- y5.19$samples$lib.size[Methylation19=="Me"] +
  y5.19$samples$lib.size[Methylation19=="Un"]
# example: before we do this, ** fyi these values outdated but method is same
# head(y$samples):
# rr01-Me lib.size is 1113815
# rr01-Un lib.size is 7079083
# after we do this (command y5$samples$lib.size <- rep(TotalLibSize, each=2))
# both rr01-Me and rr01-Un are 8192898, the sum

# exploring this more as global methylation
as.data.frame(y5.19$samples) ### export to excel? ###

# apply library size
y5.19$samples$lib.size <- rep(TotalLibSize19, each=2)
y5.19$samples

### normalize reads

############ MDS PLOT #############
# M-value
Me19 <- y5.19$counts[,Methylation19=="Me"] # extract # Me
Un19 <- y5.19$counts[,Methylation19=="Un"] # extract # Un
M19 <- log2(Me19+2) - log2(Un19+2) # calculate M-value, add 2 for zero frequency case
colnames(M19) <- Sample19

# nrow(M19) # 417036

# generate MDS plot
# set colors for visualisation
colors19 <- c(rep("blue",4),
            rep("black",5),
            rep("grey",5),
            rep("red",5))
# plot with color legend (note: takes a long time)
pdf("2021.01.29_MDS_allCpGs_19samples.pdf")
plotMDS(M19, col=colors19, gene.selection="pairwise", main="All CpG sites") + legend(
  "topright",
  bty = "n",
  c("CORT", "VEHICLE", "CONTROL", "CRS"),
  fill = c("blue", "black", "grey", "red")
)
dev.off()

###################################################
###################################################

###################################################
###################################################


########
# design matrix
designSL19 <- model.matrix(~0+Population, data=targets19)
# colnames(designSL19) # "PopulationCONTROL" "PopulationCORT" "PopulationCRS" "PopulationVEHICLE"
colnames(designSL19) <- c("CONTROL", "CORT", "CRS", "VEHICLE")

design19 <- modelMatrixMeth(designSL19)
#  each DNA sample generates two counts, a count of methylated reads and a count of unmethylated reads, for each genomic locus for each sample. The function converts sample-level information about the treatment conditions to make an appropriate design matrix with two rows for each sample. Counts are assumed to be ordered as methylated and then unmethylated by sample.
# nrow(design19)
# ncol(design19)

# dispersion
ydisp19 <- estimateDisp(y5.19, design=design19, trend="none") # takes long long
ydisp19$common.dispersion # 0.4916338 very high
summary(ydisp19$prior.df) # infinity for all 
# Inf means all CpG-wise dispersions are equal to the common dispersion

################ differential methylation analysis ################ 
fit19 <- glmFit(ydisp19, design19)
# comparison CORT
contr.CORT19 <- makeContrasts(CORT19=CORT-VEHICLE, levels=design19)
# comparison CRS
contr.CRS19 <- makeContrasts(CRS19=CRS-CONTROL, levels=design19)

# test CORT
lrt.CORT19 <- glmLRT(fit19,contrast=contr.CORT19)
# test CRS
lrt.CRS19 <- glmLRT(fit19,contrast=contr.CRS19)

# view top differentially methylated CpG sites
topTags(lrt.CORT19)
topTags(lrt.CORT19, n=3)

# show number of differentially methylated CpG sites at FDR of 5%
summary(decideTests(lrt.CORT19, method = "bonferroni")) # only 3 -- all downregulate
# 3 genes for CORT are:
DMC_CORT19 <- topTags(lrt.CORT19)[1:3,] # same genes as before, as expected
# hypo_CORT <- DMC_CORT
# save gene names as file
# write_lines(hypo_CORT$table$Symbol, "/Volumes/jaeyoonc/_THESIS/JCthesis/3hypoCORTgenes.txt")

# Sec13, Nat14, Anks6
### note: Nat14 is the only one with distance -810
### Nat14 is one of the top hypomethylated CpGs found in one study, Littlejohn et al.

# for CRS
summary(decideTests(lrt.CRS19, method = "bonferroni")) # 790 down; 917 up, tot 1707
# take just the 1707 diff methylated
DMC_CRS19 <- topTags(lrt.CRS19, n=2000)[1:1707,]
# separate for hypo & hypermethylated
hypo_CRS19 <- DMC_CRS19[DMC_CRS19$table$logFC < 0, ] # hypomethyl nrow(hypo_CRS19) # 790
hyper_CRS19 <- DMC_CRS19[DMC_CRS19$table$logFC > 0, ] # hypermethyl nrow(hyper_CRS19) # 917

# list of the 735 gene symbols
hypo_CRS_gene <- unique(hypo_CRS19$table$Symbol)
length(hypo_CRS_gene) # 790 loci --> 617 unique genes

hyper_CRS_gene <- unique(hyper_CRS19$table$Symbol)
length(hyper_CRS_gene) # 917 loci --> 782 unique genes
# 
# # save gene names as files
# write_lines(hypo_CRS_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/313hypoCRSgenes.txt")
# write_lines(hyper_CRS_gene, "/Volumes/jaeyoonc/_THESIS/JCthesis/307hyperCRSgenes.txt")
# 
# sum(unique(y5$genes$Symbol) == "Pot1a") # includes Tert
# maybe not TERT but diff methylation of sth that regulates indirectly telomerase
# made list manually with some genes of interest from NCBI
telomeres <- c("Tert", "Terc", "Myc", "Pot1a", "Smg5", "Tep1", "Rtel1", "Terf2", "Terf1", 
               "Terb1", "Ctc1", "Dnmt1", "Telo2", "Sod2", "Terf2ip",
               "Terb2", "Cdk2", "Snai1", "Tet2", "Ten1",
               "Pot1b", "Tet1", "Tlq1", "Sde2", "Tinf2",
               "Dnmt3l","Tet3", "Tel-rs3", "Tel-rs4", "Tel-rs2",
               "Nek2", "Tel-rs1", "Tel-rs7", "Tep1", "Kmt5c")

telomeres %in% unique(y5$genes$Symbol)
telomeres %in% hypo_CRS_gene
telomeres %in% hyper_CRS_gene #!!! Terf2 is here!!!


"Pot1a" %in% hypo_CRS_gene
"Pot1a" %in% hyper_CRS_gene

"Bdnf" %in% hypo_CRS_gene
"Bdnf" %in% hyper_CRS_gene

"Mbtps1" %in% hypo_CRS_gene
"Mbtps1" %in% hyper_CRS_gene


# visualise
tiff("2021.02.08.CORT19_logFC.tiff", 
     width=550, height=500)
plotMD(lrt.CORT19, main="CORT vs VEHICLE
       methylation across all CpG sites", ylim=c(-10,10))
dev.off()

tiff("2021.02.08.CRS19_logFC.tiff",
     width=550, height=500)
plotMD(lrt.CRS19, main="CRS vs CONTROL
       methylation across all CpG sites")
dev.off()


#####################################
#####################################
#####################################

# making manhattan -- exported for use in the other R file
# need df x with bp, chr, p, gene (snp)

#CRS
chr <- lrt.CRS19$genes$Chr
chr <- str_remove(chr, "chr") # chr1 --> 1
bp <- lrt.CRS19$genes$Locus
p <- lrt.CRS19$table$PValue
gene <- lrt.CRS19$genes$Symbol
# length(chr) #417036
# length(bp)
# length(p) length of all is 417036
xcrs19 <- data.frame(chr,bp,p,gene)
writexl::write_xlsx(xcrs19,"/Volumes/jaeyoonc/_THESIS/JCthesis/ManhattanCRS19.xlsx")

# CORT -- no need to change bc CORT stays same.
chrcort <- lrt.CORT19$genes$Chr
chrcort <- str_remove(chrcort, "chr") # chr1 --> 1
bpcort <- lrt.CORT19$genes$Locus
pcort <- lrt.CORT19$table$PValue
genecort <- lrt.CORT19$genes$Symbol
# length(chrcort) # also 417036
# length(bp)
# length(p) length of all is 406421
xcort19 <- data.frame(chrcort,bpcort,pcort,genecort)
writexl::write_xlsx(xcort19,"/Volumes/jaeyoonc/_THESIS/JCthesis/ManhattanCORT19.xlsx")

#####################################
#####################################
#####################################

# by chromosome
ChrIndices <- list()
for (i in ChrNames19) ChrIndices[[i]] <- which(ydisp19$genes$Chr==i)
bychrCORT19 <- fry(ydisp19, index=ChrIndices, design=design19, contrast=contr.CORT19)
bychrCRS19 <- fry(ydisp19, index=ChrIndices, design=design19, contrast=contr.CRS19)

#####################################
# global methylation patterns around TSS
# set 20kb as distance around TSS
# CRS (the 5 CRSc sequences)
i19 <- abs(fit19$genes$Distance) < 20000
loCRS19 <- lowess(fit19$genes$Distance[i19], fit19$coefficients[i19, "CRS"], f=0.3) # f is smoothness, default 2/3
# visualise
tiff("2021.02.08.CRS19_DistTSSMethylation.tiff", 
     width=550, height=500)
  plot(loCRS19, type ="l", xlab="Distance to TSS", 
       ylab="Logit methylation level", 
       main="CRS",
       col="red",
       ylim=c(-6,3))
  abline(h=0, lty=2)
  abline(v=0, lty=2)
dev.off()

# CORT
loCORT19 <- lowess(fit19$genes$Distance[i19], fit19$coefficients[i19, "CORT"], f=0.3) # f is smoothness, default 2/3
# visualise
tiff("2021.02.08.CORT19_DistTSSMethylation.tiff", 
     width=550, height=500)
plot(loCORT19, type ="l", xlab="Distance to TSS", 
     ylab="Logit methylation level", 
     main="CORT",
     col="blue",
     ylim=c(-6,3))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()


# CONTROL
loCTRL19 <- lowess(fit19$genes$Distance[i19], fit19$coefficients[i19, "CONTROL"], f=0.3) # f is smoothness, default 2/3
# visualise
tiff("2021.02.08.CONTROL19_DistTSSMethylation.tiff", 
     width=550, height=500)
plot(loCTRL19, type ="l", xlab="Distance to TSS", 
     ylab="Logit methylation level", 
     main="CONTROL",
     col="grey",
     ylim=c(-6,3))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()



# VEHICLE
loVEH19 <- lowess(fit19$genes$Distance[i19], fit19$coefficients[i19, "VEHICLE"], f=0.3) # f is smoothness, default 2/3
# visualise
tiff("2021.02.08.VEHICLE19_DistTSSMethylation.tiff", 
     width=550, height=500)
plot(loVEH19, type ="l", xlab="Distance to TSS", 
     ylab="Logit methylation level", 
     main="VEHICLE",
     col="black",
     ylim=c(-6,3))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()


### methylation changes to TSS position!
# CORT
i.CORTvsVEH19 <- abs(lrt.CORT19$genes$Distance) < 80000
lo.CORTvsVEH19 <- lowess(lrt.CORT19$genes$Distance[i19], lrt.CORT19$table$logFC[i19], f=0.3) # f is smoothness, default 2/3

# don't need to actually repeat as CORT and VEH remain same as before
pdf("Final-CORTvsVEH_TSSglobalmethyl.pdf")
plot(lo.CORTvsVEH19, type ="l", xlab="Distance to TSS", 
     ylab="Change in methylation level (logFC)", 
     main="CORT vs VEHICLE",
     col="blue",
     ylim=c(-1.5,0.2),
     xlim=c(-22000,22000))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

i.CRSvsCTRL19 <- abs(lrt.CRS19$genes$Distance) < 80000
lo.CRSvsCTRL19 <- lowess(lrt.CRS19$genes$Distance[i.CRSvsCTRL19], lrt.CRS19$table$logFC[i.CRSvsCTRL19], f=0.3) # f is smoothness, default 2/3

pdf("Final-CRSvsCTRL_TSSglobalmethyl.pdf")
plot(lo.CRSvsCTRL19, type ="l", xlab="Distance to TSS", 
     ylab="Change in methylation level (logFC)", 
     main="CRS vs CONTROL",
     col="red",
     ylim=c(-1.5,0.2),
     xlim=c(-22000,22000))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

##########################################################
# differential methylation in PROMOTERS

hypo_CRS_InPromoter <- hypo_CRS19$table$Distance >= -1000 & hypo_CRS19$table$Distance <= 2000 
hypo_CRS_ip <- hypo_CRS19$table[hypo_CRS_InPromoter,,]
nrow(hypo_CRS_ip) # 163
hyper_CRS_InPromoter <- hyper_CRS19$table$Distance >= -1000 & hyper_CRS19$table$Distance <= 2000 
hyper_CRS_ip <- hyper_CRS19$table[hyper_CRS_InPromoter,,]
nrow(hyper_CRS_ip) # 55

telomeres %in% unique(y5.19$genes$Symbol)
telomeres %in% hyper_CRS19$table$Symbol #!!! Terf2 is here!!!
telomeres %in% hyper_CRS_ip # Terf2 no longer here

#############################################
# 2021.02.08. below code unedited.
# if wish to do heatmap for the DMCs for the 19 samples (417036 CpG sites),
# use "M19" M values

# make heatmap -- based on M values

Mvals <- data.frame(M)
nrow(Mvals) - sum(rownames(Mvals) == rownames(y5$counts)) # 0, loci are all same
# also make sure row names are the same for Mvals and gene names
nrow(Mvals) - sum(rownames(Mvals) == paste0("chr",rownames(y5$genes))) # 0 good
y5genes <- y5$genes$Symbol
# bind together so Mvals includes gene symbol names
Mvals <- (cbind(Mvals, y5genes))

# now we want to only select for CRS and the 735 differentially methylated genes
rrCRS <- c("rr11", "rr12", "rr13", "rr14", "rr15", "rr16c", "rr17c", "rr18c", "rr19c", "rr20", "rr21c", "y5genes")
MCRS <- Mvals[,rrCRS]
# select only for the loci that are DMCs
DMCposCRS <- rownames(DMC_CRS$table) # length(DMCposCRS) 735 yes
MCRS_DMC <- MCRS[DMCposCRS,] # only taking the 735 DMCs..
# length(unique(MCRS_DMC$y5genes)) -- there are 611 total diff methylated genes

MCRS_DMCmap <- MCRS_DMC %>% dplyr::select(-y5genes) # take out genes
MCRS_DMCmap <- t(MCRS_DMCmap)
ncol(MCRS_DMCmap) # 735 (DMCs)
nrow(MCRS_DMCmap) # 11 conditions (rr's)

pdf("2020.12.04_Heatmap_CRS_735DMCsonly.pdf")
Heatmap(MCRS_DMCmap, 
        show_row_names = TRUE, row_title = "Hippocampal sample",
        show_column_names = FALSE, column_title = "CRS vs CONTROL 735 DMCs")
dev.off()

# what if take all CpG sites, not just the DMCs? -- !!! this crashes
# !!! below code crashes too many CpGs
nrow(MCRS) # 406421 CpGs
MCRS_map <- MCRS %>% dplyr::select(-y5genes) # take out genes
MCRS_map <- t(MCRS_map)
ncol(MCRS_map) # 406421 CpGs
nrow(MCRS_map) # 11 conditions (rr's)

Heatmap(MCRS_map, 
        show_row_names = TRUE, row_title = "Mouse sample",
        show_column_names = FALSE, column_title = "DMCs")


