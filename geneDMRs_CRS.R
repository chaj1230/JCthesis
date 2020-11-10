# install
if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "ffbase", "genomation", "KEGG.db", 
                       "pheatmap", "plotrix", "qqman", "RCircos", "VennDiagram", "org.Mm.eg.db"))
# using the annotation dataset for enrichment, e.g., "org.Mm.eg.db" of mouse
source("https://install-github.me/xiaowangCN/GeneDMRs")



