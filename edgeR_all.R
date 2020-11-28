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
# set group name = population # NOTE HERE paper is INCORRECT, need rep each 2x!!!
yall$samples$group <- factor(rep(targets$Population, each=2)) 
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

# store methylation cols for ease later
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))

# coverage: returns each position and coverage for each sample
Coverage <- yall$counts[, Methylation=="Me"] + yall$counts[, Methylation =="Un"]
class(Coverage)
head(Coverage[,1])

############################################################################3
###### VISUALISE COVERAGE ######
# convert to df for visualisation of coverage
cover <- as.data.frame(Coverage)

# get summary statistics for each
a <- summary(cover[,1][(cover[,1] > 0)])
for (i in 1:ncol(cover)){
  b <- summary(cover[,i][(cover[,i] > 0)])
  all <- rbind(a, b)
  a <- all
}
# remove the duplicate first row of rr01
coveragestats <- a[-1,]
rownames(coveragestats) <- Sample

# histogram density plot
# hist(cover[,i][   (cover[,i] > 0) & (cover[,i] < 50)   ])
# hist(cover[,13][   (cover[,13] > 0) & (cover[,13] < 50)   ])
# 
# d <- density(cover[,13][   (cover[,13] > 0) & (cover[,13] < 100)   ])
# d <- density(cover[,1][   (cover[,1] > 0) & (cover[,1] < 100)   ])
pdf("coveragedensity.pdf")
for (i in 1:ncol(cover)){
  # d <- density(cover[,i][   (cover[,i] > 0) & (cover[,i] < 100)   ])
  d <- density(cover[,i][   (cover[,i] > 0)   ])
  plot(d, main = Sample[i], xlim=c(1,100))
}
dev.off()

# cumulative
pdf("coveragecumulative.pdf")
for (i in 1:ncol(cover)){
  e <- ecdf(cover[,i][   (cover[,i] > 0)   ])
  plot(e, main = Sample[i], cex = 0, xlim=c(1, 10))
}
dev.off()

### FILTER COVERAGE
# change this with filter
y <- yall

##################################
# set library size to be equal for each pair (Me and Un) of libraries
# this is bc each pair is treated as a unit in analysis

# add up for each sample lib size for methylated with unmethylated
TotalLibSize <- y$samples$lib.size[Methylation=="Me"] +
                y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples

############ MDS PLOT #############
# M-value
Me <- y$counts[,Methylation=="Me"] # extract # Me
Un <- y$counts[,Methylation=="Un"] # extract # Un
M <- log2(Me+2) - log2(Un+2) # calculate M-value, add 2 for zero frequency case
colnames(M) <- Sample


# toy example
# toy <- M[1:6,1:6]
# toycolors <- c(rep("blue",3),
#                rep("black",3))
# plotMDS(toy, col=toycolors)
# plotMDS(toy, col=toycolors) + legend(
#   "topleft",
#   bty = "n",
#   c("toy", "toy2"),
#   fill = c("blue", "black")
# )

# generate MDS plot
# set colors for visualisation
colors <- c(rep("blue",5),
            rep("black",5),
            rep("grey",5),
            rep("red",6),
            rep("yellow",5),
            rep("orange", 5))
# plot with color legend
plotMDS(M, col=colors) + legend(
    "topleft",
    bty = "n",
    c("CORT", "VEHICLE", "CONTROL", "CRS", "CRS_2", "CRSc"),
    fill = c("green", "black", "grey", "red","blue","purple")
  )
  

  
