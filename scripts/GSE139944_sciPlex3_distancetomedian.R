# https://github.com/cole-trapnell-lab/sci-plex

#### Packages ####
library(scran)
library(monocle3)
library(dplyr)
library(ggplot2)

#### Functions ####

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

#### DM analysis ####
# Parameters
cell.type <- "K562"
product <- "Luminespib"
d <- 10000

sample.cells <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
  filter(cell_type == cell.type, dose == d) %>% 
  rownames()

sample.counts <- filt.counts[, sample.cells]

# Normalisation by RPM
sample.rpm <- sample.counts / (colSums(sample.counts) / 10^6)

# Distance to median metric
means.sample <- rowMeans(sample.rpm)

cv2.sample <- apply(sample.rpm, 1, var) / (means.sample ^ 2)

dm.stat.sample <- DM(means.sample, cv2.sample)
dm.stat.sample <- dm.stat.sample[!(is.na(dm.stat.sample))]
dm.stat.sample <- dm.stat.sample[order(abs(dm.stat.sample), decreasing=TRUE)]
