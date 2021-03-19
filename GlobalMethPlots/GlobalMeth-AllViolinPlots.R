setwd("/Volumes/jaeyoonc/_THESIS/JCthesis/GlobalMethPlots")

library(dplyr)
library(ggplot2)

### ORIGINAL BATCH ###
# read in data
original21 <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/GlobalMethPlots/2020.10.19.GlobalMeth_OriginalBatch.txt", header = TRUE)
# without the outlier sample 5, CORT
original20 <- original21[-10,]

# separate by CORT+VEHICLE and CRS+CONTROL
CORT <- original20[original20$treatment %in% c("CORT", "VEHICLE"), ]
CRS <- original20[original20$treatment %in% c("CRS", "CONTROL"), ]

# violin plot visualization of all 21 samples
p <- ggplot(original21, aes(x=as.factor(treatment), y=methylCpG, color=treatment)) +
  geom_violin() +
  ggtitle("Global methylation of original 21 samples") +
  ylab("Methylated CpG (%)") +
  # keep order so that x-axis not sorted alphabetically
  scale_x_discrete(limits=original21$treatment)
p + stat_summary(fun=median, geom="point", shape = 23, size=1, color="black")
# export graph
jpeg("ViolinPlots_original21.jpeg")
p + geom_boxplot(width=0.05) + 
  stat_summary(fun=median, geom="point", shape = 23, size=1, color="black") # add diamond for mean
dev.off()  

####### violin plot visualization with 20 samples, sample 5 outlier removed #######
po <- ggplot(original20, aes(x=treatment, y=methylCpG, color=treatment)) +
  geom_violin() +
  scale_x_discrete(limits=original20$treatment) +
  ggtitle("Global methylation of original 20 samples") +
  ylab("Methylated CpG (%)")
# add diamond for mean
po + stat_summary(fun=median, geom="point", shape = 23, size=1, color="black")
# export graph
jpeg("ViolinPlots_original20.jpeg")
po + geom_boxplot(width=0.05) + 
  stat_summary(fun=median, geom="point", shape = 23, size=1, color="black") # add diamond for mean
dev.off()  

## one-way ANOVA across 20 samples without sample 5 outlier
fit_original20 <- aov(original20$methylCpG ~ as.factor(original20$treatment))
summary(fit_original20) # p-value is 0.000104

# ANOVA between CORT and VEHICLE
fit_originalCORT <- aov(CORT$methylCpG ~ CORT$treatment)
summary(fit_originalCORT) # p-value is 0.584

# ANOVA between CRS and CONTROL
fit_originalCRS <- aov(CRS$methylCpG ~ CRS$treatment)
summary(fit_originalCRS) # p-value is 0.00343

#####################################################

## SECOND BATCH WITH CRS second batch only

# read in data
batch2.21 <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/GlobalMethPlots/2021.02.07.GlobalMeth_SecondBatch.txt", header = TRUE)
# without the outlier sample 5, CORT
batch2.20 <- batch2.21[-10,]

# separate by CORT+VEHICLE and CRS+CONTROL
batch2.CORT <- batch2.20[batch2.20$treatment %in% c("CORT", "VEHICLE"), ]
# batch2.CORT == CORT
batch2.CRSreseq <- batch2.20[batch2.20$treatment %in% c("CRS", "CONTROL"), ]

# violin plot for batch2 (no rr05) using reseq CRS
p.batch2 <- ggplot(batch2.20, aes(x=treatment, y=methylCpG, color=treatment)) +
  geom_violin() +
  scale_x_discrete(limits=batch2.20$treatment) +
  ggtitle("Global methylation of 19 samples with CRS second batch") +
  ylab("Methylated CpG (%)") +
  ylim(25, 53)
# view graph
p.batch2 + stat_summary(fun=median, geom="point", shape = 23, size=1, color="black")
# export graph
jpeg("ViolinPlots_Batch2_CRSreseqONLY.jpeg")
p.batch2 + geom_boxplot(width=0.05) + 
  stat_summary(fun=median, geom="point", shape = 23, size=1, color="black") # add diamond for mean
dev.off()  

## one-way ANOVA across 20 samples without sample 5 outlier
fit_batch2c <- aov(batch2.20$methylCpG ~ as.factor(batch2.20$treatment))
summary(fit_batch2c) # p-value is 1.52e-05

# ANOVA between CORT and VEHICLE
fit_batch2.CORT <- aov(batch2.CORT$methylCpG ~ batch2.CORT$treatment)
summary(fit_batch2.CORT) # p-value is 0.584

# ANOVA between CRS and CONTROL
fit_batch2.CRSreseq <- aov(batch2.CRSreseq$methylCpG ~ batch2.CRSreseq$treatment)
summary(fit_batch2.CRSreseq) # p-value is 0.000326


#####################################################

## SECOND BATCH WITH CRSc

# read in data
batch2.CRSc21 <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/GlobalMethPlots/2020.12.06.GlobalMeth_CRScombined.txt", header = TRUE)
# without the outlier sample 5, CORT
batch2.CRSc20 <- batch2.CRSc21[-10,]

# separate by CORT+VEHICLE and CRS+CONTROL
batch2.CORT <- batch2.CRSc20[batch2.CRSc20$treatment %in% c("CORT", "VEHICLE"), ]
  # batch2.CORT == CORT
batch2.CRSc <- batch2.CRSc20[batch2.CRSc20$treatment %in% c("CRS", "CONTROL"), ]

# violin plot for batch2.CRSc20 (no rr05)
p.batch2c <- ggplot(batch2.CRSc20, aes(x=treatment, y=methylCpG, color=treatment)) +
  geom_violin() +
  scale_x_discrete(limits=batch2.CRSc20$treatment) +
  ggtitle("Global methylation of 20 samples with CRS combined") +
  ylab("Methylated CpG (%)") +
  ylim(25, 53)
# view graph
p.batch2c + stat_summary(fun=median, geom="point", shape = 23, size=1, color="black")
# export graph
jpeg("ViolinPlots_Batch2_CRSc20.jpeg")
p.batch2c + geom_boxplot(width=0.05) + 
  stat_summary(fun=median, geom="point", shape = 23, size=1, color="black") # add diamond for mean
dev.off()  

## one-way ANOVA across 20 samples without sample 5 outlier
fit_batch2c <- aov(batch2.CRSc20$methylCpG ~ as.factor(batch2.CRSc20$treatment))
summary(fit_batch2c) # p-value is 0.39

# ANOVA between CORT and VEHICLE
fit_batch2.CORT <- aov(batch2.CORT$methylCpG ~ batch2.CORT$treatment)
summary(fit_batch2.CORT) # p-value is 0.584

# ANOVA between CRS and CONTROL
fit_batch2.CRSc <- aov(batch2.CRSc$methylCpG ~ batch2.CRSc$treatment)
summary(fit_batch2.CRSc) # p-value is 0.211

###
###
###

# without rr20, the only CRS that was not resequenced
batch2.CRSc19 <- batch2.CRSc20[-19,] # row 19 of CRSc20 corresponds to rr20

# separate by CORT+VEHICLE and CRS+CONTROL
batch2_19.CORT <- batch2.CRSc19[batch2.CRSc19$treatment %in% c("CORT", "VEHICLE"), ]
# batch2.CORT == CORT
batch2_19.CRSc <- batch2.CRSc19[batch2.CRSc19$treatment %in% c("CRS", "CONTROL"), ]

# violin plot for batch2.CRSc19 (no rr05 and no rr20)
p.batch2c19 <- ggplot(batch2.CRSc19, aes(x=treatment, y=methylCpG, color=treatment)) +
  geom_violin()+
  scale_x_discrete(limits=batch2.CRSc20$treatment) +
  ggtitle("Global methylation of 19 samples with CRS combined") +
  ylab("Methylated CpG (%)") +
  ylim(25, 53) +

# export graph
jpeg("ViolinPlot_Final19.jpeg")
p.batch2c19 + geom_boxplot(width=0.05) + 
  scale_color_manual(values=c("grey62", "blue", "red", "black")) +
  stat_summary(fun=median, geom="point", shape = 23, size=1, color="black")
  # add diamond for mean
dev.off()  






## one-way ANOVA across 19 samples without rr05 and rr20
fit_batch2c19 <- aov(batch2.CRSc19$methylCpG ~ as.factor(batch2.CRSc19$treatment))
summary(fit_batch2c19) # p-value is 0.0201

# ANOVA between CORT and VEHICLE
fit_batch2_19.CORT <- aov(batch2_19.CORT$methylCpG ~ batch2_19.CORT$treatment)
summary(fit_batch2_19.CORT) # p-value is 0.584, same as before

# ANOVA between CRS and CONTROL
fit_batch2_19.CRSc <- aov(batch2_19.CRSc$methylCpG ~ batch2_19.CRSc$treatment)
summary(fit_batch2_19.CRSc) # p-value is 0.0248


