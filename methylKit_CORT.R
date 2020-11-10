# working, most updated code with methylKit_CRS
# later make a function with it and run it for CORT vs VEHICLE, CRS vs CORT

library(methylKit)

# cort.list has rr01-rr05 CORT and rr06-rr10 VEHICLE
cort.list=list("~/Desktop/THESIS/covs/G695_rr01_S15_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr02_S16_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr03_S17_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr04_S18_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr05_S19_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr06_S20_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr07_S21_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr08_S22_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr09_S23_val_1_bismark_bt2_pe.bismark.cov",
               "~/Desktop/THESIS/covs/G695_rr10_S24_val_1_bismark_bt2_pe.bismark.cov")

# cort is CORT and VEHICLE
cort=methRead(cort.list,
              sample.id = list("rr01", "rr02", "rr03", "rr04", "rr05", "rr06", "rr07", "rr08", "rr09", "rr10"),
              treatment = c(1,1,1,1,1,0,0,0,0,0),
              pipeline = "bismarkCoverage",
              header = FALSE,
              assembly="GRCm38")

getMethylationStats(cort[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(cort[[10]],plot=TRUE,both.strands=FALSE)
getCoverageStats(cort[[1]],plot=TRUE,both.strands=FALSE)





# (2.5) filter away low coverage
# discards bases that have coverage below 10X 
# also discards bases with > 99.9th percentile of coverage in each sample
cort.filtered=filterByCoverage(cort,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

getMethylationStats(cort.filtered[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(cort.filtered,plot=TRUE,both.strands=FALSE)



############ COMPARATIVE ANALYSIS ############ 
# merge samples to ensure all bases covered
# merge all samples to one object for base-pair locations that are covered in all samples

# create methylBase obj with bases/regions covered in ALL samples
cort_unite=unite(cort, destrand=FALSE)
# relax: create methylBase obj with CpGs covered in AT LEAST ONE 1 sample per group
cort_unite.min=unite(cort, min.per.group=1L)

# same thing for filtered:
cort.filtered_unite=unite(cort.filtered, destrand=FALSE)
cort.filtered_unite.min=unite(cort.filtered,min.per.group=1L)


### sample correlation
getCorrelation(cort_unite,plot=TRUE)

### cluster
clusterSamples(cort_unite, dist="correlation", method="ward", plot=TRUE)
clusterSamples(cort_unite.min, dist="correlation", method="ward", plot=TRUE)

clusterSamples(cort.filtered_unite, dist="correlation", method="ward", plot=TRUE)
clusterSamples(cort.filtered_unite.min, dist="correlation", method="ward", plot=TRUE)

PCASamples(cort_unite)
PCASamples(cort.filtered_unite)

PCASamples(cort_unite, screeplot=TRUE)


#################### Find DMRs (3.6) ####################

# tests for differential methylation using logistic regression
## (calculateDiffMeth takes long time)
myDiff=calculateDiffMeth(cort_unite) # 15 warnings: glm.fit: fitted probabilities numerically 0 or 1 occurred

myDiff.filt=calculateDiffMeth(cort.filtered_unite)




########### 2020.11.02 attempting to sort cov file by chr ###########
# # sort rr05 then read into file.list
# rr05_sort <- rr05[with(rr05, order(V1)), ]
# 
# # rename rows
# rownames(rr05_sort) <- c(1:nrow(rr05))
# 
# # save
# write.table(rr05_sort, file = "~/Desktop/THESIS/nomismatch.rr05.sorted.cov", 
#             quote = FALSE, row.names = FALSE, sep = "\t")
# 
# # here, went into Terminal and deleted the header column names in sorted.cov
# 
# # now make new file.list.sorted
# file.list_sort=list("~/Desktop/THESIS/nomismatch.rr05.sorted.cov")
# class(file.list_sort)
# 
# # this works!
# myobj_sort=methRead(file.list_sort,
#                     sample.id=list("rr05"),
#                     pipeline = "bismarkCoverage",
#                     header = FALSE,
#                     assembly="GRCm38",treatment=1)

# but -- think that coverage file does NOT need to be sorted for methylKit, so no need.

################# inputing BAM files #################

# my.meth.bam<-processBismarkAln(
#   "~/Desktop/THESIS/G695_rr01_S15_val_1_bismark_bt2_pe.sorted.bam",
#   sample.id = "rr01",
#   assembly = "mm10")
