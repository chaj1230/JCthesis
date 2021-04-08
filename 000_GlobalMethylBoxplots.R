library(dplyr)
library(ggplot2)

# read in global CpG methyl data from Bismark
globalmethyl <- read.table("~/Desktop/globalmethyl.txt", header=TRUE)
# order
globalmethyl$names <- factor(globalmethyl$Population , levels=c("VEHICLE", "CORT", "CONTROL", "CRS"))

# separate by CORT+VEHICLE and CRS+CONTROL
CORT <- globalmethyl[globalmethyl$Population %in% c("VEHICLE", "CORT"), ]
CORT$names <- factor(CORT$Population , levels=c("VEHICLE", "CORT"))

CRS <- globalmethyl[globalmethyl$Population %in% c("CONTROL", "CRS"), ]
CRS$names <- factor(CRS$Population , levels=c("CONTROL", "CRS"))

# boxplot CORT vs VEHICLE
pdf("cort.pdf",
    width=5,
    height=5)
boxplot(CORT$methylCpG ~ CORT$names, 
        data=CORT,
        scale_x_discrete(limits=CORT$Population),
        xlab="Treatment",
        ylab="Methylated CpG (%)",
        col=c("#7AD56A", "#3e8c3d"))
dev.off()

# boxplot CRS vs CONTROL
pdf("crs.pdf",
    width=5,
    height=5)
boxplot(methylCpG ~ Population, 
        data=CRS,
        scale_x_discrete(limits=CRS$Population),
        xlab="Treatment",
        ylab="Methylated CpG (%)",
        col=c("#d487e5", "#8239b5"),
        ylim=c(25,50))
dev.off()

# boxplot all 4 groups
pdf("all-globalmeth.pdf",
    width=7,
    height=5)
boxplot(globalmethyl$methylCpG ~ globalmethyl$names, 
        data=globalmethyl,
        scale_x_discrete(limits=globalmethyl$Population),
        xlab="Treatment",
        ylab="Methylated CpG (%)",
        col=c("#7AD56A", "#3e8c3d","#d487e5", "#8239b5"),
        ylim=c(25,50))
dev.off()

# kruskal test for significance
kruskal.test(methylCpG ~ Population, data = CORT) 
# 0.6242

kruskal.test(methylCpG ~ Population, data = CRS) 
# 0.0472