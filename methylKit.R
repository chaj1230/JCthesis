# 2020.10.29 after meeting with Lisa
# 2020.11.01 tried to sort by chr (V1) then feed into methRead but not working
# 2020.11.01 trying to sort by chr, saving as a new cov file, then re-reading into R, but again not working

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

getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)

getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)

# (2.5) filter away low coverage
# discards bases that have coverage below 10X 
# also discards bases with > 99.9th percentile of coverage in each sample
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

getMethylationStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)



# # sort rr05 then read into file.list
rr05_sort <- rr05[with(rr05, order(V1)), ]
write.table(rr05_sort, file = "~/Desktop/THESIS/nomismatch.rr05.sorted.cov")
rr05_sort <- read.table("~/Desktop/THESIS/nomismatch.rr05.sorted.cov")

nrow(rr05_sort) # 4938416 rows :) same as rr05
ncol(rr05_sort) # 6 cols :)
class(rr05)
class(rr05_sort)

# now make new file.list.sorted
file.list_sort=list("~/Desktop/THESIS/nomismatch.rr05.sorted.cov")
class(file.list_sort)

# again error...
myobj_sort=methRead(file.list_sort,
                    sample.id=list("rr05"),
                    pipeline = "bismarkCoverage",
                    header = FALSE,
                    assembly="GRCm38",treatment=1)




