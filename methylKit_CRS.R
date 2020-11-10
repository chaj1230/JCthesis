library(methylKit)

# crs.list has rr11-rr05 CORT and rr06-rr10 VEHICLE
crs.list=list("~/Desktop/THESIS/covs/G695_rr11_S25_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr12_S26_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr13_S27_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr14_S28_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr15_S29_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr16_S30_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr17_S31_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr18_S32_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr19_S33_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr20_S34_val_1_bismark_bt2_pe.bismark.cov",
              "~/Desktop/THESIS/covs/G695_rr21_S35_val_1_bismark_bt2_pe.bismark.cov")

# crs is CONTROL and CRS
crs=methRead(crs.list,
              sample.id = list("rr11", "rr12", "rr13", "rr14", "rr15", "rr16", "rr17", "rr18", "rr19", "rr20", "rr21"),
              treatment = c(0,0,0,0,0,1,1,1,1,1,1),
              pipeline = "bismarkCoverage",
              header = FALSE,
              assembly="GRCm38")

getMethylationStats(crs[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(crs[[6]],plot=TRUE,both.strands=FALSE)
# note: rr20 (10th obj here) has significantly lower % methylation per base at 100% compared to the rest of CRS
# rr20 at 100% looks more like CTRL
getCoverageStats(crs[[1]],plot=TRUE,both.strands=FALSE)


# (2.5) filter away low coverage
# discards bases that have coverage below 10X 
# also discards bases with > 99.9th percentile of coverage in each sample
crs.filtered=filterByCoverage(crs,lo.count=10,lo.perc=NULL,
                               hi.count=NULL,hi.perc=99.9)

getMethylationStats(crs.filtered[[6]],plot=TRUE,both.strands=FALSE)

getCoverageStats(crs.filtered[[1]],plot=TRUE,both.strands=FALSE)


############ COMPARATIVE ANALYSIS ############ 
# merge samples to ensure all bases covered
# merge all samples to one object for base-pair locations that are covered in all samples

# create methylBase obj with bases/regions covered in ALL samples
crs_unite=unite(crs, destrand=FALSE)
# relax: create methylBase obj with CpGs covered in AT LEAST ONE 1 sample per group
crs_unite.min=unite(crs, min.per.group=1L)

# same thing for filtered:
crs.filtered_unite=unite(crs.filtered, destrand=FALSE)
crs.filtered_unite.min=unite(crs.filtered,min.per.group=1L)


### sample correlation
getCorrelation(crs_unite,plot=TRUE)

### cluster
clusterSamples(crs_unite, dist="correlation", method="ward", plot=TRUE)
clusterSamples(crs_unite.min, dist="correlation", method="ward", plot=TRUE)

clusterSamples(crs.filtered_unite, dist="correlation", method="ward", plot=TRUE)
clusterSamples(crs.filtered_unite.min, dist="correlation", method="ward", plot=TRUE)

PCASamples(crs_unite)
PCASamples(crs.filtered_unite)


#################### Find DMRs (3.6) ####################

# tests for differential methylation using logistic regression
crs_myDiff=calculateDiffMeth(crs_unite, mc.cores=2)

crs_myDiff.filt=calculateDiffMeth(crs.filtered_unite, mc.cores=2) 

### view by chromosome
# for unfiltered
diffMethPerChr(crs_myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
# generate dataframe, copy to excel
View(as.data.frame(diffMethPerChr(crs_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)))

# for filtered
diffMethPerChr(crs_myDiff.filt,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
View(as.data.frame(diffMethPerChr(crs_myDiff.filt,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)))





