# https://github.com/cole-trapnell-lab/sci-plex

#### Packages ####
library(scran)
library(monocle3)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(gdata)
library(vioplot)
library(glue)

#### Functions ####
dm_analysis <- function(counts, top.n) {
  means <- rowMeans(counts)
  cv2 <- apply(counts, 1, var) / (means ^ 2)
  
  dm <- DM(means, cv2)
  dm <- dm[!(is.na(dm))]
  dm <- dm[order(abs(dm), decreasing=TRUE)]
  
  return(head(dm, top.n))
}

correlation_analysis <- function(counts, top.n.dm) {
  return(rcorr(t(as.matrix(counts[names(top.n.dm), ])),
               type="spearman"))
}

dcorr <- function(corr) {
  sqrt((1 - upperTriangle(corr$r)) / 2)
}

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin")

#### Load data ####
col.data <- read.csv("./processed/GSE139944/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
protein.coding <- read.csv("./data/gencode.v27.transcripts.bed", 
                           sep="\t", quote='', header=FALSE) %>% 
  filter(V9 == "protein_coding") %>%
  .$V7 %>% unique()

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./processed/GSE139944/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name), ])

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose, rownames(cell.freq))]

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(readRDS("./processed/GSE139944/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, rownames(filt.col.data)]

filt.gene.data <- fData(readRDS("./processed/GSE139944/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, ]

rownames(filt.counts) <- filt.gene.data$gene_short_name
rownames(filt.gene.data) <- filt.gene.data$gene_short_name

#### DM and correlation analysis ####
dm.corr.list <- list()

# Parameters
cell.type <- "K562"
product <- "Luminespib"
d <- c(10, 100, 1000, 10000)

top.n <- 500

vehicle.cells <- filt.col.data[filt.col.data$vehicle, ] %>%
  filter(cell_type == cell.type) %>%
  rownames()
vehicle.counts <- filt.counts[, vehicle.cells]

# RPM normalisation
vehicle.rpm <- vehicle.counts / (colSums(vehicle.counts) / 10^6)

# Exclude genes with low mean expression
vehicle.rpm <- vehicle.rpm[rowMeans(vehicle.rpm) > 10, ]

# Get distance-to-median metric
dm.corr.list[["Vehicle_0_topdm"]] <- dm_analysis(vehicle.rpm, top.n)

# Get correlation metric
dm.corr.list[["Vehicle_0_corr"]] <- correlation_analysis(vehicle.rpm,
                                                         dm.corr.list[["Vehicle_0_topdm"]])

# Repeat for product
for (i in 1:length(d)) {
  sample.cells <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
    filter(cell_type == cell.type, dose == d[i]) %>% 
    rownames()
  sample.counts <- filt.counts[, sample.cells]
  
  # RPM normalisation
  sample.rpm <- sample.counts / (colSums(sample.counts) / 10^6)
  
  # Exclude genes with low mean expression
  sample.rpm <- sample.rpm[rowMeans(sample.rpm) > 10, ]
  
  # Get distance-to-median metric
  dm.corr.list[[paste0(product, "_", d[i], "_topdm")]] <- dm_analysis(sample.rpm, top.n)
  
  # Get correlation metric
  dm.corr.list[[paste0(product, "_", d[i], "_corr")]] <- correlation_analysis(sample.rpm,
                                                                              dm.corr.list[[paste0(product, "_", d[i], "_topdm")]])
}

# Plot violin plot
png(glue("./processed/GSE139944/transcriptional_noise/sciPlex3_DM_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
vioplot(dcorr(dm.corr.list[["Vehicle_0_corr"]]), 
        dcorr(dm.corr.list[[paste0(product, "_10_corr")]]),
        dcorr(dm.corr.list[[paste0(product, "_100_corr")]]),
        dcorr(dm.corr.list[[paste0(product, "_1000_corr")]]),
        dcorr(dm.corr.list[[paste0(product, "_10000_corr")]]),
        names=c("Vehicle_0", 
                paste0(product, "_10"),
                paste0(product, "_100"),
                paste0(product, "_1000"),
                paste0(product, "_10000")), 
        col="orange")
title(main=glue("{cell.type} transcriptional heterogeneity in {product} treatments
                for the top N={top.n} highly variable genes"),
      ylab="Pairwise distance measure")
dev.off()

