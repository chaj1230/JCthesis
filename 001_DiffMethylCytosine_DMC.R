setwd("/Volumes/jaeyoonc/_THESIS/JCthesis/")

# load packages
library(edgeR)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
# genome wide annotation for mouse
library(org.Mm.eg.db) 
library(tidyverse) #change from chr1 to 1; vice versa etc
library(qqman) # manhattan

# read in txt file of sample names
targets19 <- read.delim("targets19_noRR20.txt", row.names = "Sample", stringsAsFactors=FALSE)
# nrow(targets19) # 19 rows
# ncol(targets19) # 2 cols
Sample19 <- row.names(targets19)

# files is a list of all 19 sample file locations
files19 <- c(
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
  # remove rr20 outlier
  "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")

# yall is all 19 cov files, read and collate counts
yall19 <- readBismark2DGE(files19, sample.names=Sample19)

# set group name = population # need rep each 2x for normalization
yall19$samples$group <- factor(rep(targets19$Population, each=2)) 
# also add "chr" to chr numbers
row.names(yall19$counts) <- paste0("chr",row.names(yall19$counts))
yall19$genes$Chr <- paste0("chr",yall19$genes$Chr)
nrow(yall19) # 10,027,860 cytosines

# filter MT and Y chr
keep19 <- rep(TRUE, nrow(yall19))
Chr19 <- as.character(yall19$genes$Chr)
sort(unique(Chr19)) # no unassembled chr
# remove Y chr and mito DNA
keep19[Chr19=="chrY"] <- FALSE
keep19[Chr19=="chrMT"] <- FALSE
table(keep19) # removing 20971 cytosines
yall19 <- yall19[keep19,, keep.lib.sizes=FALSE]
nrow(yall19) # 10,006,889

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

### FILTER COVERAGE: set threshold of 5x in all samples
keep5.19 <- rowSums(Coverage19 >= 5) == 19 # 433,721 loci

# filter out CpGs that are never or always methylated bc no info about differential methyl
Me19 <- yall19$counts[,Methylation19=="Me"] # extract # Me
Un19 <- yall19$counts[,Methylation19=="Un"] # extract # Un
HasBoth19 <- rowSums(Me19) > 0 & rowSums(Un19) > 0

# apply both filters
y5.19 <- yall19[keep5.19 & HasBoth19,, keep.lib.sizes=FALSE]

# all the genes captured now with this coverage threshold
length(unique(y5.19$genes$Symbol)) # 17778 unique genes
nrow(y5.19$genes) # 417,036 loci

# set library size to be equal for each pair (Me and Un) of libraries
# this is bc each pair is treated as a unit in analysis
# add up for each sample lib size for methylated with unmethylated
TotalLibSize19 <- y5.19$samples$lib.size[Methylation19=="Me"] +
  y5.19$samples$lib.size[Methylation19=="Un"]

# apply library size for normalization
y5.19$samples$lib.size <- rep(TotalLibSize19, each=2)
y5.19$samples

############ MDS PLOT #############
# M-value
Me19 <- y5.19$counts[,Methylation19=="Me"] # extract # Me
Un19 <- y5.19$counts[,Methylation19=="Un"] # extract # Un
M19 <- log2(Me19+2) - log2(Un19+2) # calculate M-value, add 2 for zero frequency case
colnames(M19) <- Sample19

# generate MDS plot
# set colors for visualisation
colors19 <- c(rep("blue",4),
            rep("black",5),
            rep("grey",5),
            rep("red",5))

# plot with color legend
pdf("2021.01.29_MDS_allCpGs_19samples.pdf")
plotMDS(M19, col=colors19, gene.selection="pairwise", main="All CpG sites") + legend(
  "topright",
  bty = "n",
  c("CORT", "VEHICLE", "CONTROL", "CRS"),
  fill = c("blue", "black", "grey", "red")
)
dev.off()

################################################

# create design matrix
designSL19 <- model.matrix(~0+Population, data=targets19)
colnames(designSL19) <- c("CONTROL", "CORT", "CRS", "VEHICLE")
design19 <- modelMatrixMeth(designSL19)

# dispersion
ydisp19 <- estimateDisp(y5.19, design=design19, trend="none")
ydisp19$common.dispersion
summary(ydisp19$prior.df) # infinity for all 
# Inf means all CpG-wise dispersions are equal to the common dispersion

################ differential methylation analysis ################ 

# use generalized linear model to do differential methylation calling
fit19 <- glmFit(ydisp19, design19)
# comparison CORT
contr.CORT19 <- makeContrasts(CORT19=CORT-VEHICLE, levels=design19)
# comparison CRS
contr.CRS19 <- makeContrasts(CRS19=CRS-CONTROL, levels=design19)

# test CORT
lrt.CORT19 <- glmLRT(fit19,contrast=contr.CORT19)

# DMR::get all the hypo and hyper regardless of their p val
CORT_hypo_all <- lrt.CORT19[lrt.CORT19$table$logFC < 0,]$table # 241113
CORT_hyper_all <- lrt.CORT19[lrt.CORT19$table$logFC > 0,]$table #175923
# take the locus and p-value
CORT_hypo_all_df <- CORT_hypo_all %>% select(PValue)
CORT_hyper_all_df <- CORT_hyper_all %>% select(PValue)
# download for DMR calling with comb-p on command line
write.csv(CORT_hypo_all_df, "~/Desktop/jCORT_241113hypo.csv")
write.csv(CORT_hyper_all_df, "~/Desktop/jCORT_175923hyper.csv")

# test CRS
lrt.CRS19 <- glmLRT(fit19,contrast=contr.CRS19)
sum(lrt.CRS19$table$logFC < 0)

# DMR::get all the hypo and hyper regardless of their p val
CRS_hypo_all <- lrt.CRS19[lrt.CRS19$table$logFC < 0,]$table # 307321
CRS_hyper_all <- lrt.CRS19[lrt.CRS19$table$logFC > 0,]$table # 109715
# take the locus and p-value
CRS_hypo_all_df <- CRS_hypo_all %>% select(PValue)
CRS_hyper_all_df <- CRS_hyper_all %>% select(PValue)
# download for DMR calling with comb-p on command line
write.csv(CRS_hypo_all_df, "~/Desktop/jCRS_307321hypo.csv")
write.csv(CRS_hyper_all_df, "~/Desktop/jCRS_109715hyper.csv")

# view top differentially methylated DMCs
topTags(lrt.CORT19) # only 3
topTags(lrt.CORT19, n=3)

# show number of differentially methylated CpG sites at FDR of 5%
summary(decideTests(lrt.CORT19, method = "bonferroni")) # only 3 -- all hypo

# for CRS
summary(decideTests(lrt.CRS19, method = "bonferroni")) # 790 down; 917 up, tot 1707
# take just the 1707 diff methylated
DMC_CRS19 <- topTags(lrt.CRS19, n=2000)[1:1707,]
# separate for hypo & hypermethylated
hypo_CRS19 <- DMC_CRS19[DMC_CRS19$table$logFC < 0, ] # hypomethyl nrow(hypo_CRS19) # 790
hyper_CRS19 <- DMC_CRS19[DMC_CRS19$table$logFC > 0, ] # hypermethyl nrow(hyper_CRS19) # 917

hypo_CRS19df <- as.data.frame(hypo_CRS19)
hyper_CRS19df <- as.data.frame(hyper_CRS19)

# visualise MD plots
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
# manhattan plot
# need format: bp, chr, p, gene (identifier)
# load packages
library(dplyr)
library(qqman)

# CRS
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

### CRS - read in M values
xCRS <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/ManhattanCRS19.txt", header = TRUE)
colnames(xCRS) <- c("CHR", "BP", "P", "SNP")

# generate manhattan plot
# not annotating genes
pdf("/Volumes/jaeyoonc/_THESIS/ManhattanCRS_allCpGs.pdf",
    width=6.50,
    height=3.10)
manhattan(xCRS, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
          # set threshold as 9.05e-05 with multiple hypothesis test correction bonferroni
          # this corresponds to FDR of ~0.05
          suggestiveline = -log10(9.05e-05), 
          genomewideline= -log10(1e-07),
          ## (for readability, annotate gene names for the)
          ## (top 30 most significantly DMCs, corresponding p threshold = 1.55e-07)
          # annotatePval= 1e-07, # 9.039034e-05 is the p-value of the 735th
          # for readability, annotate the top for each chr only
          col = c("#d487e5", "#8239b5"),
          # annotateTop=TRUE,
          main="CRS vs CONTROL all CpG sites",
          cex=0.7) # font size
dev.off()

# CORT
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

### CORT - read in M values for all 400k cytosines
xCORT <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/ManhattanCORT19.txt", header = TRUE)
colnames(xCORT) <- c("CHR", "BP", "P", "SNP")

# generate manhattan plot
pdf("/Volumes/jaeyoonc/_THESIS/ManhattanCORT_allCpGs.pdf",
    width=6.50,
    height=3.10)
manhattan(xCORT, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
          # set threshold as 5e-07 with multiple hypothesis test correction bonferroni
          # this corresponds to FDR of ~0.05
          suggestiveline = -log10(5e-07), genomewideline=-log10(1e-07),
          # annotatePval=4e-07,
          col = c("#7AD56A", "#3e8c3d"),
          main="CORT vs VEHICLE all CpG sites",
          cex=0.7) # change size of font
dev.off()

# q-q plot
pdf("2021.02.13_CORT_Q-Q_allCpGs.pdf")
qq(xCORT$P, main="CORT vs VEHICLE Q-Q plot")
dev.off()


#####################################
# examine patterns of differential methylation

# by chromosome
ChrIndices <- list()
for (i in ChrNames19) ChrIndices[[i]] <- which(ydisp19$genes$Chr==i)
bychrCORT19 <- fry(ydisp19, index=ChrIndices, design=design19, contrast=contr.CORT19)
bychrCRS19 <- fry(ydisp19, index=ChrIndices, design=design19, contrast=contr.CRS19)

# global methylation patterns around TSS
# set 20kb as distance around TSS
i19 <- abs(fit19$genes$Distance) < 20000
# CRS
loCRS19 <- lowess(fit19$genes$Distance[i19], fit19$coefficients[i19, "CRS"], f=0.3) # f is smoothness, default 2/3
# visualise
tiff("CRS19_DistTSSMethylation.tiff", 
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
tiff("CORT19_DistTSSMethylation.tiff", 
     width=550, height=500)
plot(loCORT19, type ="l", 
     xlab="Distance to TSS", 
     ylab="Logit methylation level", 
     cex=100,
     main="CORT",
     cex.lab=1.5, # make fonts bigger
     cex.axis=1.5, 
     cex.main=2,
     col="#3E8C3D",
     lwd=3, # make line thicker
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

### methylation changes around TSS position!
# CORT
i.CORTvsVEH19 <- abs(lrt.CORT19$genes$Distance) < 80000
lo.CORTvsVEH19 <- lowess(lrt.CORT19$genes$Distance[i19], lrt.CORT19$table$logFC[i19], f=0.3) # f is smoothness, default 2/3

pdf("CORTvsVEH_TSSglobalmethyl.pdf")
plot(lo.CORTvsVEH19, type ="l", 
     xlab="Distance to TSS", 
     ylab="Change in methylation level (logFC)", 
     main="CORT vs VEHICLE",
     cex.lab=1.3, # make fonts bigger
     cex.axis=1.1, 
     cex.main=1.7,
     col="#3E8C3D",
     lwd=3, # make line thicker
     ylim=c(-1.5,0.2),
     xlim=c(-22000,22000))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

# CRS
i.CRSvsCTRL19 <- abs(lrt.CRS19$genes$Distance) < 80000
lo.CRSvsCTRL19 <- lowess(lrt.CRS19$genes$Distance[i.CRSvsCTRL19], lrt.CRS19$table$logFC[i.CRSvsCTRL19], f=0.3) # f is smoothness, default 2/3

pdf("CRSvsCTRL_TSSglobalmethyl.pdf")
plot(lo.CRSvsCTRL19, type ="l", 
     xlab="Distance to TSS", 
     ylab="Change in methylation level (logFC)", 
     main="CRS vs CONTROL",
     cex.lab=1.3, # make fonts bigger
     cex.axis=1.1, 
     cex.main=1.7,
     col="#8239B5",
     lwd=3, # make line thicker
     ylim=c(-1.5,0.2),
     xlim=c(-22000,22000))
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()