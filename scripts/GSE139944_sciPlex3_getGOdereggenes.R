#### Packages ####
library(glue)
library(dplyr)

#### Config ####
setwd("~/mrc/project/scrna-seq")

# Parameters
cell.type <- ""
treatment <- ""

alpha <- 0.05

#### Load data ####
gene.set <- read.table(glue("processed/GSE139944/deseq/DEGtable/sciPlex3_DESeq_{cell.type}_{treatment}_vs_Vehicle_0.txt"), header=TRUE, sep='\t',
                       row.names=1, check.names=FALSE)

#### Filtering ####
updownreg.set <- gene.set %>% filter(padj < alpha)

#### Save to list ####
fileConn <- file("processed/GSE139944/GO/GSE139944_sciPlex3_{treatment}_upregset.txt")
writeLines(rownames(updownreg.set %>% filter(log2FoldChange > 0)), fileConn)
close(fileConn)

fileConn <- file("processed/GGSE139944/GO/GSE139944_sciPlex3_{treatment}_downregset.txt")
writeLines(rownames(updownreg.set %>% filter(log2FoldChange < 0)), fileConn)
close(fileConn)

fileConn <- file("processed/GGSE139944/GO/GSE139944_sciPlex3_{treatment}_updownregset.txt")
writeLines(rownames(updownreg.set), fileConn)
close(fileConn)


