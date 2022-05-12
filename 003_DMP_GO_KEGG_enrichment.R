setwd("/Volumes/jaeyoonc/_THESIS/JCthesis")
# note: need source f(x) of clusterProfiler from their GitHub for dotplot font formatting changes
# this is called 003_SOURCE_FUNCTION_clusterProfiler-dependencyfunctions.R in GITHUB

# FIGURES 14 AND 15 IN THESIS

# load required packages
library(clusterProfiler)
library(org.Mm.eg.db) 
library(stringr)
# for dependencies
library(ggplot2)
library(DOSE)

# read in the 1970 hypomethylated DMP genes
genes1970 <- readLines("1970hypoCRS19genesENTREZID.txt") # this is in github

# read in the 119 hypermethylated DMP genes
hyper119 <- readLines("119hyperCRS19genesENTREZID.txt")

# read in all 17155 analyzed genes
allgenes <- readLines("all17155genes19_ENTREZID.txt")
# take RANDOM SAMPLE of 1970 of these 
set.seed(0)
random1970 <- sample(allgenes, 1970, replace=FALSE)
random119 <- sample(allgenes, 119, replace=FALSE)
length(random1970)

#####################

# write a function to get GO term enrichment per category
GOenrich <- function(ENTREZIDs, category){
  # run enrichGO with parameters
  ego <- enrichGO(gene          = ENTREZIDs,
                  # universe is all mm10 genes for enrichment bkground
                  OrgDb         = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont           = category,
                  # bonferroni p < 0.05 and q < 0.05
                  pAdjustMethod = "bonferroni",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}


######################### HYPO ############################

####### GO cellular component
egoCC <- GOenrich(genes1970, "CC") # 16 GO terms
CC.df <- head(egoCC, n=100)
# export
write.csv(CC.df, "/Volumes/jaeyoonc/_THESIS/JCthesis/cellularComponent.csv")

# visualize dotplot
pdf("CC.pdf", height=5.5,width=6)
dotplot(egoCC, 
        title = "GO Cellular Component",
        label_format = 30,
        font.size=12)
dev.off()

# random CC
randomCC <- GOenrich(random1970, "CC")
# 1 enriched term: "GO:0005667" "transcription regulator complex"
# p adjusted after bonferroni: 0.0355
# generation: 59/1766


####### GO biological process
egoBP <- GOenrich(genes1970, "BP") # 44 GO terms

# simplify
# reduce redundancy; GO terms with semantic similarity higher than 0.7 are redundant
egoBPsimp <- clusterProfiler::simplify(
  x=egoBP, 
  cutoff=0.7, 
  by="p.adjust", 
  select_fun=min)

nrow(egoBPsimp) # 23
BP.df.simp <- head(egoBPsimp, n=100)
# cutoff 0.7 gives 23 GO terms (from initially 44)

# export
write.csv(BP.df.simp, "/Volumes/jaeyoonc/_THESIS/JCthesis/biologicalProcessSIMP.csv")

# visualize all 44 BP terms
pdf("BP-all.pdf")
dotplot(egoBP, 
        title = "GO Biological Process",
        label_format = 30,
        font.size=12)
dev.off()

# visualize 23 nonredundant GO terms
pdf("BP-simp.pdf",
    height=9,
    width=6)
dotplot(egoBPsimp, 
        title = "GO Biological Process",
        label_format = 30,
        font.size=12)
dev.off()

# BP random
randomBP <- GOenrich(random1970, "BP")
randomBP # 2 terms but they're the same essentially:
# "regulation of chromosome organization" 
# "positive regulation of chromosome organization"
# simplify:
randomBPsimp <- clusterProfiler::simplify(
  x=randomBP, 
  cutoff=0.7, 
  by="p.adjust", 
  select_fun=min)
randomBPsimp # "regulation of chromosome organization"
# adjusted p-val=0.0222


####### GO molecular function
egoMF <- GOenrich(genes1970, "MF")
head(egoMF, 50) # 5 terms
egoMF 

# visualize
pdf("MF.pdf", height=3.5, width=6)
dotplot(egoMF, 
        title = "GO Molecular Function",
        label_format = 25,
        font.size=12)
dev.off()

# export
MF.df <- head(egoMF, n=100)
write.csv(MF.df, "/Volumes/jaeyoonc/_THESIS/JCthesis/molecularFunction.csv")

# random MF
randomMF <- GOenrich(random1970, "MF")
randomMF # nothing

##### KEGG terms
kegg <- enrichKEGG(gene       = genes1970,
                   organism     = 'mmu',
                   pAdjustMethod = "bonferroni",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
kegg.df <- head(kegg, n=100)

# export
write.csv(kegg.df, "/Volumes/jaeyoonc/_THESIS/JCthesis/KEGG.csv")
nrow(kegg.df) # 19 KEGG terms
# visualize dotplot
pdf("kegg.pdf",height=9,width=6)
dotplot(kegg, showCategory = 100, title = "KEGG Pathways",
        label_format = 26, font.size = 12)
dev.off()

# for random KEGG
randokegg <- enrichKEGG(gene       = random1970,
                        organism     = 'mmu',
                        pAdjustMethod = "bonferroni",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
# "Human papillomavirus infection" 
# "Thyroid hormone signaling pathway" 
# "N-Glycan biosynthesis"

######################### HYPER ############################

# GO terms
hyperCC <- GOenrich(hyper119, "CC") # 0
hyperBP <- GOenrich(hyper119, "BP") # 0
hyperMF <- GOenrich(hyper119, "MF") # 0

# for random 119 genes
randomCC <- GOenrich(random119, "CC") # 0
randomBP <- GOenrich(random119, "BP") # 0
randomMF <- GOenrich(random119, "MF") # 0

# KEGG
hyperkegg <- enrichKEGG(gene         = hyper119,
                   organism     = 'mmu',
                   pAdjustMethod = "bonferroni",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) # 0
