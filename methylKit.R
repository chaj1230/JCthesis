# 2020.10.29 after meeting with Lisa

library(methylKit)
# read in rr05
rr05 <- read.table("~/Desktop/THESIS/nomismatch.rr05.cov")
nrow(rr05) # 4938416 rows
ncol(rr05) # 6 cols

file.list=list("~/Desktop/THESIS/nomismatch.rr05.cov")

myobj=methRead(file.list,
               sample.id=list("rr05"),
               pipeline = "bismarkCoverage",
               header = FALSE,
               assembly="GRCm38",treatment=1)
# (2.4) percent methylation statistics
getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)
# histogram for percent methylation distribution
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
# read coverage per base information
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)

# (2.5) filter away low coverage
# discards bases that have coverage below 10X 
# also discards bases with > 99.9th percentile of coverage in each sample
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
# again display graphs for methylation distribution and per base coverage
getMethylationStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)

