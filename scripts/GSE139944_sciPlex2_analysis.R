#### Packages ####
library(monocle3)
library(Matrix)

#### Functions ####
aggregate_cols <- function(counts, col.data) {
  t.counts <- as.data.frame(t(counts))
  t.counts$treatment <- col.data$top_oligo
  
  t.counts.agg <- aggregate(t.counts[1:length(t.counts) - 1], 
                            list(t.counts$treatment), 
                            mean)
  
  counts.agg <- t(t.counts.agg)
  colnames(counts.agg) <- counts.agg[1, ]
  counts.agg <- as.data.frame(counts.agg[-1, ])
  
  counts.agg.mat <- as.matrix(sapply(counts.agg, as.numeric))
  rownames(counts.agg.mat) <- rownames(counts.agg)
  
  return(counts.agg.mat)
  # 
}

#### Config ####
setwd("~/mrc/project/scrna-seq")

#### Load data ####
cds <- readRDS("./processed/GSE139944/GSM4150377_sciPlex2_cds.RDS")
proteo.list <- read.csv("~/mrc/project/rna-seq/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

# Extract relevant data from S4 Object
counts <- cds@assays$data$counts
col.data <- data.frame(cds@colData)

# Format gene names in count matrix
rownames(counts) <- gsub("\\.[0-9_A-Z]+$", "", rownames(counts))

#### Filtering ####
# Get list of CORE proteostasis genes
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for only HDAC inhibitors
filt.col.data <- col.data[grep("SAHA", col.data$top_oligo), ]

# Filter for genes and cells of interest
filt.counts <- as.data.frame(as.matrix(counts[which(rownames(counts) %in% proteo.genes),
                             which(colnames(counts) %in% filt.col.data$Cell)]))

# Collapse matrix by average expression per treatment per well
agg.filt.counts <- aggregate_cols(filt.counts, filt.col.data)


