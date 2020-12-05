library(dplyr)
library(qqman) # manhattan

# note: x has chrX replaced with 20
x <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/manhattan.txt", header = TRUE)
colnames(x) <- c("CHR", "BP", "P")
p01 <- x %>% filter(P < 0.05)

manhattan(p01, chr = "CHR", bp = "BP", p = "P")


xCORT <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/manhattanCORT.txt", header = TRUE)
colnames(xCORT) <- c("CHR", "BP", "P", "SNP")
# pCORT <- xCORT %>% filter(P < 0.05)
# manhattan(pCORT, chr = "CHR", bp = "BP", p = "P")

DMC_CORT <- c("Sec13", "Nat14", "Anks6")
manhattan(xCORT, chr = "CHR", bp = "BP", p = "P", snp = "SNP", highlight = "Sec13")
manhattan(xCORT, chr = "CHR", bp = "BP", p = "P", snp = "SNP", highlight = "Nat14")
manhattan(xCORT, chr = "CHR", bp = "BP", p = "P", snp = "SNP", highlight = "Anks6")

#########

xCRS <- read.table("/Volumes/jaeyoonc/_THESIS/JCthesis/manhattancrs.txt", header = TRUE)
colnames(xCRS) <- c("CHR", "BP", "P", "SNP")

pCRS <- xCRS %>% filter(P < 0.3)

manhattan(pCRS, chr = "CHR", bp = "BP", p = "P", snp = "SNP")
qq(xCRS$P)

hist(xCRS$P,
     main="P-values of 406421 differentially 
     methylated CpGs in CRS",
     xlab = "p-values",
     ylab = "frequency",
     col = "red")


hist(xCORT$P,
     main="P-values of 406421 differentially 
     methylated CpGs in CORT",
     xlab = "p-values",
     ylab = "frequency",
     col = "blue")

