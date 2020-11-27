library(edgeR)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
# genome wide annotation for mouse
library(org.Mm.eg.db) 

# read in txt file of samples key
targets <- read.delim("targets.txt", row.names = "Sample", stringsAsFactors=FALSE)
# nrow(targets) #31 rows
# ncol(targets) # 2 cols
Sample <- row.names(targets)
# files is a list of all 31 sample file locations
files <- c(
  # 21 original batch
  "~/Desktop/THESIS/covs/G695_rr01_S15_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr02_S16_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr03_S17_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr04_S18_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr05_S19_val_1_bismark_bt2_pe.bismark.cov",
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
  "~/Desktop/THESIS/covs/G695_rr16_S30_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr17_S31_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr18_S32_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr19_S33_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr20_S34_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/G695_rr21_S35_val_1_bismark_bt2_pe.bismark.cov",
  # second batch CRS
  "~/Desktop/THESIS/covs/read2_covs/G695_rr16_S48_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/read2_covs/G695_rr17_S49_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/read2_covs/G695_rr18_S50_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/read2_covs/G695_rr19_S51_R1_001_val_1_bismark_bt2_pe.bismark.cov",                   
  "~/Desktop/THESIS/covs/read2_covs/G695_rr21_S52_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  # combined CRS
  "~/Desktop/THESIS/covs/combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr17_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr18_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr19_R1_001_val_1_bismark_bt2_pe.bismark.cov",
  "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")

# yall is all 31 cov files, read and collate counts
yall <- readBismark2DGE(files, sample.names=Sample)
yall$samples$group <- factor(targets$Population) # group name = population
# also add "chr" to chr numbers
row.names(yall$counts) <- paste0("chr",row.names(yall$counts))
yall$genes$Chr <- paste0("chr",yall$genes$Chr)

# nrow(yall) is 10766174
# filter MT and Y chr
keep <- rep(TRUE, nrow(yall))
Chr <- as.character(yall$genes$Chr)
sort(unique(Chr)) # there are no unassembled chr
# remove Y chr and mito DNA
keep[Chr=="Y"] <- FALSE
keep[Chr=="MT"] <- FALSE
table(keep) # FALSE is 22701; TRUE is 10743473 
yall <- yall[keep,, keep.lib.sizes=FALSE]
# nrow(yall) is 10743473

# sort the DGEList into order, chr1 to X
ChrNames <- paste0("chr",c(1:19, "X"))
yall$genes$Chr <- factor(yall$genes$Chr, levels=ChrNames)
o <- order(yall$genes$Chr, yall$genes$Locus)
yall <- yall[o,]

################################################################

# gene annotation
TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Mm") # mouse
# nrow(TSS) == nrow(yall) TRUE both 10743473
yall$genes$EntrezID <- TSS$gene_id
yall$genes$Symbol <- TSS$symbol
yall$genes$Strand <- TSS$strand
yall$genes$Distance <- TSS$distance
yall$genes$Width <- TSS$width
head(yall$genes)

          
          
