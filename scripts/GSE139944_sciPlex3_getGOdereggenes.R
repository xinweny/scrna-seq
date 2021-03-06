#### Packages ####
library(glue)
library(dplyr)

#### Config ####
setwd("~/mrc/project/scrna-seq")

# Parameters
cell.type <- "K562"
treatment <- "Trichostatin A_10000"

alpha <- 0.05

#### Load data ####
gene.set <- read.table(glue("./data/GSE139944/deseq/DEGtable/sciPlex3_DESeq_{cell.type}_{treatment}_vs_Vehicle_0.txt"), header=TRUE, sep='\t',
                       row.names=1, check.names=FALSE)

rownames(gene.set) <- gsub("\\.[0-9_A-Z]+$", "", rownames(gene.set))

#### Filtering ####
updownreg.set <- gene.set %>% filter(padj < alpha)

#### Save to list ####
fileConn <- file(glue("./data/GSE139944/GO/gene_sets/sciPlex3_upregset_{cell.type}_{treatment}.txt"))
writeLines(rownames(updownreg.set %>% filter(log2FoldChange > 0)), fileConn)
close(fileConn)

fileConn <- file(glue("./data/GSE139944/GO/gene_sets/sciPlex3_downregset_{cell.type}_{treatment}.txt"))
writeLines(rownames(updownreg.set %>% filter(log2FoldChange < 0)), fileConn)
close(fileConn)

fileConn <- file(glue("./data/GSE139944/GO/gene_sets/sciPlex3_updownregset_{cell.type}_{treatment}.txt"))
writeLines(rownames(updownreg.set), fileConn)
close(fileConn)


