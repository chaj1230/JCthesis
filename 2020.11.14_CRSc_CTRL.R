# CRS combined vs CONTROL

library(methylKit)
library(genomation)

# crs.list has rr11-rr15 CTRL & rr16-rr19 and rr21 combined and rr20 original CRS
CRSc.list=list("~/Desktop/THESIS/covs/G695_rr11_S25_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr12_S26_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr13_S27_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr14_S28_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr15_S29_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/combined/G695_rr17_R1_001_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/combined/G695_rr18_R1_001_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/combined/G695_rr19_R1_001_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr20_S34_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")

CRSc=methRead(CRSc.list,
             sample.id = list("rr11", "rr12", "rr13", "rr14", "rr15", "rr16c", "rr17c", "rr18c", "rr19c", "rr20", "rr21c"),
             treatment = c(0,0,0,0,0,1,1,1,1,1,1),
             pipeline = "bismarkCoverage",
             header = FALSE,
             assembly="GRCm38")

getMethylationStats(CRSc[[2]],plot=FALSE,both.strands=FALSE)
# getMethylationStats(crs.comb[[6]],plot=TRUE,both.strands=FALSE)
# 
# getCoverageStats(crs.comb[[1]],plot=TRUE,both.strands=FALSE)
# 
# getCoverageStats(crs.comb[[10]],plot=TRUE,both.strands=FALSE)


# filter away low coverage
# discard bases that have coverage below 10X 
# discard bases with > 99.9th percentile of coverage in each sample
CRSc.filt=filterByCoverage(CRSc,lo.count=10,lo.perc=NULL,
                              hi.count=NULL,hi.perc=99.9)

# getMethylationStats(crs.filtered[[6]],plot=TRUE,both.strands=FALSE)
# 
# getCoverageStats(crs.filtered[[1]],plot=TRUE,both.strands=FALSE)
# 

############ COMPARATIVE ANALYSIS ############ 
# merge samples to ensure all bases covered
# merge all samples to one object for base-pair locations that are covered in all samples

# create methylBase obj with bases/regions covered in ALL samples
CRSc_unite=unite(CRSc, destrand=FALSE)
# relax: create methylBase obj with CpGs covered in AT LEAST ONE 1 sample per group
# CRSc_unite.min=unite(CRSc, min.per.group=1L)

# same thing for filtered:
CRSc.filt_unite=unite(CRSc.filt, destrand=FALSE)
# crs.filtered_unite.min=unite(crs.filtered,min.per.group=1L)


### sample correlation
# getCorrelation(crs_unite,plot=TRUE)

### cluster
clusterSamples(CRSc_unite, dist="correlation", method="ward", plot=TRUE)
# clusterSamples(crs_unite.min, dist="correlation", method="ward", plot=TRUE)

clusterSamples(CRSc.filt_unite, dist="correlation", method="ward", plot=TRUE)
# clusterSamples(CRSc.filt_unite, dist="correlation", method="ward", plot=TRUE)
# other clustering methods avail

PCASamples(CRSc_unite)
PCASamples(CRSc.filt_unite)


#################### Find DMRs ####################

# id DMR using logistic regression
CRSc_myDiff=calculateDiffMeth(CRSc_unite, mc.cores=2)

CRSc.filt_myDiff=calculateDiffMeth(CRSc.filt_unite, mc.cores=2) 

### view by chromosome
# for unfiltered
diffMethPerChr(CRSc_myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=10)
# generate dataframe, copy to excel
View(as.data.frame(diffMethPerChr(crs_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)))

# for filtered
diffMethPerChr(CRSc.filt_myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=10)
diffMethPerChr(CRSc.filt_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

View(as.data.frame(diffMethPerChr(crs_myDiff.filt,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)))


############# Annotating DMRs #############

# annotate differentially methylated CpGs with promoter/exon/intron
# with RefSeq bed from UCSC genome browser
genes=readTranscriptFeatures("~/Desktop/THESIS/GRCm38_geneAnnotation")

annotateWithGeneParts(as(CRSc.filt_myDiff,"GRanges"),genes)

# CpG island bed file from UCSC genome browser

cpgs=readFeatureFlank("~/Desktop/THESIS/GRCm38_CpGIslands",
                      feature.flank.name=c("CpGi","shores"))

annotateWithFeatureFlank(as(CRSc.filt_myDiff,"GRanges"),
                         cpgs$CpGi,cpgs$shores,
                         feature.name="CpGi",flank.name="shores")

