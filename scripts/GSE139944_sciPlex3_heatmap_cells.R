#### Packages ####
library(monocle3)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(glue)
library(gridExtra)
library(grid)

#### Functions ####

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

set.seed(1)

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin",
                "Trichostatin A")

#### Load data ####
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
protein.coding <- read.csv("./data/GSE139944/data/gencode.v27.transcripts.bed", 
                           sep="\t", quote='', header=FALSE) %>% 
  filter(V9 == "protein_coding") %>%
  .$V7 %>% unique()

proteo.list <- read.csv("./data/GSE139944/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./data/GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name), ])

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose, rownames(cell.freq))]

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(readRDS("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, rownames(filt.col.data)]

filt.gene.data <- fData(readRDS("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, ]

# Set rownames for matrix and gene metadata
rownames(filt.counts) <- filt.gene.data$gene_short_name
rownames(filt.gene.data) <- filt.gene.data$id

# RPM normalisation
filt.rpms <- filt.counts / (colSums(filt.counts) / 10^6)

#### Parameters ####
cell.type <- "K562"
product <- "Alvespimycin"
d <- c(10, 100, 1000, 10000)

ncells.vehicle <- 250
nrep.vehicle <- 3

vehicles <- paste0("Vehicle_0_rep", 1:nrep.vehicle)
samples <- c(vehicles, paste0(product, "_", d))

#### Coefficient of variation ####
# Get list of cells for each sample
sample.cells.list <- list()

for (i in 1:length(samples)) {
  if (grepl("Vehicle_0", samples[i])) {
    # Randomly sample vehicle cells
    sample.cells.list[[samples[i]]] <- col.data[col.data$vehicle, ] %>%
      filter(cell_type == cell.type, dose == 0) %>%
      sample_n(ncells.vehicle) %>%
      rownames()
  } else {
    sample.cells.list[[samples[i]]] <- col.data[grep(product, col.data$product_name), ] %>%
      filter(cell_type == cell.type, dose == d[i - nrep.vehicle]) %>% 
      rownames() 
  }
}

# Get cell number information
ncells <- c(rep(ncells.vehicle, nrep.vehicle),
            stack(cell.freq)[grep(glue("{cell.type}_{product}.+_10+$"), 
                                  stack(cell.freq)$ind), 
                             "values"])
names(ncells) <- paste0(cell.type, "_", samples)

#### Heatmap ####
# Select proteostasis genes
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]
proteo.genes <- filt.gene.data[proteo.genes, "gene_short_name"]
proteo.genes <- proteo.genes[!is.na(proteo.genes)]

for (i in 1:length(samples)) {
  mat <- as.matrix(filt.rpms[proteo.genes, 
                             sample.cells.list[[samples[i]]]])
  
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]
  
  png(file=glue("./data/GSE139944/heatmap/sciPlex3_heatmap_cells_custom_{cell.type}_{samples[i]}.png"), 
      width=20000, height=20000, res=300)
  heatmap.2(log10(mat + 1),
            Rowv=TRUE,
            Colv=TRUE,
            main=glue("RPMs of {cell.type} {samples[i]} cells (N={ncol(mat)}) 
                      for n={nrow(mat)} proteostasis genes"),
            dendrogram="both",
            scale="none",
            col=magma(299),
            breaks=seq(0, 5, length.out=300),
            labCol=colnames(mat), 
            labRow=rownames(mat),
            cexRow=0.5,
            cexCol=1,
            srtCol=45,
            trace="none",
            density.info="none",
            key=TRUE, keysize=0.3,
            lhei=c(1, 45), lwid=c(1, 10),
            key.xlab=glue("log10(RPM + 1)"),
            margins=c(15, 5))
  dev.off()
}

# Colour labeling of samples
d.cols <- c('blue', 'green', 'orange', 'red') # List of available colours for labeling inhibitors

select.cells <- c()
for (i in 1:length(d)) {
  select.cells <- c(select.cells, sample.cells.list[[samples[i + nrep.vehicle]]])
}

col.cols <- rep('', length(select.cells))
d.col.strs <- c()

d.col.list <- vector(mode="list", length=length(d))
names(d.col.list) <- d

for (i in 1:length(d)) {
  col.cols[select.cells %in% rownames(filt.col.data[filt.col.data$dose == d[i], ])] <- d.cols[i]
  d.col.list[[d[i]]] <- d.cols[i]
  d.col.strs[i] <- glue("{d[i]} ({d.cols[i]})")
}

mat <- as.matrix(filt.rpms[proteo.genes,select.cells])

mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

png(file=glue("./data/GSE139944/heatmap/sciPlex3_heatmap_cells_custom_{cell.type}_{product}.png"), 
    width=40000, height=20000, res=300)
heatmap.2(log10(mat + 1),
          Rowv=TRUE,
          Colv=TRUE,
          main=glue("RPMs of {cell.type} cells (N={ncol(mat)}) treated with {product} 
                    for doses {paste(d.col.strs, sep=', ', collapse=', ')}  
                    for n={nrow(mat)} proteostasis genes"),
          dendrogram="both",
          scale="none",
          col=magma(299),
          labCol=colnames(mat), 
          labRow=rownames(mat),
          cexRow=0.5,
          cexCol=1,
          srtCol=45,
          colCol=col.cols,
          trace="none",
          density.info="none",
          key=TRUE, keysize=0.3,
          lhei=c(1, 45), lwid=c(1, 10),
          key.xlab=glue("log10(RPM + 1)"),
          margins=c(15, 5))
dev.off()
