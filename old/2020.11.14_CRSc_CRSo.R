# CRS original batch vs CRS combined

library(methylKit)
library(genomation)

# cort.list has rr01-rr05 CORT and rr06-rr10 VEHICLE
combined.list=list("~/Desktop/THESIS/covs/G695_rr16_S30_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr17_S31_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr18_S32_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr19_S33_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr21_S35_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/combined/G695_rr17_R1_001_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/combined/G695_rr18_R1_001_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/combined/G695_rr19_R1_001_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/combined/G695_rr21_R1_001_val_1_bismark_bt2_pe.bismark.cov")
               


# cort is CORT and VEHICLE
comb=methRead(combined.list,
              sample.id = list("rr16", "rr17", "rr18", "rr19", "rr21", "rr16c", "rr17c", "rr18c", "rr19c", "rr21c"),
              treatment = c(0,0,0,0,0,1,1,1,1,1),
              pipeline = "bismarkCoverage",
              header = FALSE,
              assembly="GRCm38")

getMethylationStats(comb[[5]],plot=TRUE,both.strands=FALSE)
getMethylationStats(comb[[10]],plot=TRUE,both.strands=FALSE)
getCoverageStats(comb[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(comb[[9]],plot=TRUE,both.strands=FALSE)




# (2.5) filter away low coverage
# discards bases that have coverage below 10X 
# also discards bases with > 99.9th percentile of coverage in each sample
comb.filt=filterByCoverage(comb,lo.count=10,lo.perc=NULL,
                               hi.count=NULL,hi.perc=99.9)
getMethylationStats(comb.filt[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(comb.filt[[9]],plot=TRUE,both.strands=FALSE)



############ COMPARATIVE ANALYSIS ############ 
# merge samples to ensure all bases covered
# merge all samples to one object for base-pair locations that are covered in all samples

# create methylBase obj with bases/regions covered in ALL samples
comb_unite=unite(comb, destrand=FALSE)
# relax: create methylBase obj with CpGs covered in AT LEAST ONE 1 sample per group
# cort_unite.min=unite(cort, min.per.group=1L)

# same thing for filtered:
comb.filt_unite=unite(comb.filt, destrand=FALSE)
#cort.filtered_unite.min=unite(cort.filtered,min.per.group=1L)

### sample correlation
#getCorrelation(comb_unite,plot=TRUE)

### cluster
clusterSamples(comb_unite, dist="correlation", method="ward", plot=TRUE)
#clusterSamples(co_unite.min, dist="correlation", method="ward", plot=TRUE)

clusterSamples(comb.filt_unite, dist="correlation", method="ward", plot=TRUE)
#clusterSamples(cort.filtered_unite.min, dist="correlation", method="ward", plot=TRUE)

PCASamples(comb_unite)
PCASamples(comb.filt_unite)

PCASamples(comb.filt_unite,transpose=FALSE)

PCASamples(cort_unite, screeplot=TRUE)


#################### Find DMRs (3.6) ####################

# tests for differential methylation using logistic regression
## (calculateDiffMeth takes long time)
comb_myDiff=calculateDiffMeth(comb_unite, mc.cores=2)

comb.filt_myDiff=calculateDiffMeth(comb.filt_unite, mc.cores=2)

comb_myDiff[order(comb_myDiff$meth.diff),]

### view by chromosome
# for unfiltered
diffMethPerChr(comb_myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=10)
diffMethPerChr(comb_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)
# generate dataframe, copy to excel
View(as.data.frame(diffMethPerChr(cort_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)))

# for filtered
diffMethPerChr(comb.filt_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(comb.filt_myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=10)
View(as.data.frame(diffMethPerChr(comb.filt_myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)))


############# Annotating differentially methylated bases or regions #############
# using genomation package here
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
annotateWithGeneParts(as(cort_myDiff.filt,"GRanges"),genes)

# CpG
annotateWithFeatureFlank(as(cort_myDiff.filt,"GRanges"),
                         cpgs$CpGi,cpgs$shores,
                         feature.name="CpGi",flank.name="shores")


