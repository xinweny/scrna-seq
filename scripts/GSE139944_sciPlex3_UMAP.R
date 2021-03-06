#### Packages ####
library(monocle3)
library(glue)
library(ggplot2)

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

#### Load data ####
cds <- readRDS(glue("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
proteo.list <- read.csv("./data/GSE139944/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

#### Parameters ####
cells <- c("MCF7") # K562, A549, MCF7
treatments <- c("Mocetinostat")
# HSP90i - Luminespib, Alvespimycin, Tanespimycin
# HDACi - Abexinostat, Tubastatin A, Divalproex Sodium, Sodium Phenylbutyrate, Dacinostat, Droxinostat, Tucidinostat, PCI-34051, Resminostat, Panobinostat, Trichostatin A, Quisinostat, MC1568, Belinostat, Mocetinostat, ITSA-1, TMP195, Givinostat, Entinostat, M344, Tacedinaline, AR-42, Pracinostat
doses <- c(10, 100, 1000, 10000) # 10, 100, 1000, 10000
use.genes <- c()

include.vehicle <- FALSE # TRUE or FALSE
preprocess.method <- "PCA" # PCA or LSI

#### Formatting ####
if (include.vehicle) {
  # Add vehicle cells
  treatments <- append(treatments, "Vehicle")
  doses <- append(doses, 0)
}

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

# Output frequency table of no. of cells per sample
cell.freq <- table(col.data$cell_product_dose)
col.data$cell_freq <- cell.freq[match(col.data$cell_product_dose, rownames(cell.freq))]

# Reformat sample name with cell frequency info
col.data$cell_product_dose <- paste0(col.data$cell_product_dose, " (", col.data$cell_freq, ")")

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./data/GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Get list of CORE proteostasis genes for human/mouse
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for selected cell types, inhibitors and doses
filt.col.data <- with(col.data,
                      col.data[grepl(paste(cells, sep="|", collapse="|"), col.data$cell_type) &
                               grepl(paste(treatments, sep="|", collapse="|"), col.data$product_name) &
                               grepl(paste(paste("^", doses, "$", sep=""), sep="|", collapse="|"), col.data$dose), ])

filt.counts <- counts(cds)[grep("ENSG", rownames(counts(cds))), 
                           rownames(filt.col.data)]

# Filter out mouse genes and non-targets
filt.cds <- new_cell_data_set(expression_data=filt.counts,
                              cell_metadata=filt.col.data[match(colnames(filt.counts), rownames(filt.col.data)), ],
                              gene_metadata=fData(cds)[grep("ENSG", rownames(fData(cds))), ]) 

#### Normalisation and pre-processing ####
if (length(use.genes) > 0) {
  filt.cds <- preprocess_cds(filt.cds,
                             method=preprocess.method,
                             use_genes=use.genes,
                             norm_method="log",
                             pseudo_count=1)
} else {
  filt.cds <- preprocess_cds(filt.cds,
                             method=preprocess.method,
                             norm_method="log",
                             pseudo_count=1)
}

#### Clustering ####
# Reduce dimensions using UMAP
filt.cds <- reduce_dimension(filt.cds,
                             reduction_method="UMAP",
                             preprocess_method=preprocess.method)

# Cluster cells
filt.cds <- cluster_cells(filt.cds,
                          reduction_method="UMAP",
                          num_iter=1)

#### UMAP ####
# Sample metadata
agg.samples <- unique(pData(filt.cds)$cell_product_dose)
sample.metadata <- data.frame(row.names=agg.samples,
                              cell_type=unlist(lapply(strsplit(agg.samples, "_"), `[[`, 1)),
                              product=unlist(lapply(strsplit(agg.samples, "_"), `[[`, 2)),
                              dose=gsub("\\s\\([0-9]+\\)", "", unlist(lapply(strsplit(agg.samples, "_"), `[[`, 3))),
                              n_cells=as.vector(cell.freq[match(agg.samples, names(cell.freq))]))

# Plot and save UMAP
png(file=glue("./data/GSE139944/UMAP/sciPlex3_UMAP_{paste(cells, sep='-', collapse='-')}_{paste(treatments[!treatments %in% c('Vehicle')], sep='-', collapse='-')}{if (include.vehicle) '-Vehicle' else ''}.png"), 
    width=3000, height=2000, res=300)
plot_cells(filt.cds,
           color_cells_by="cell_product_dose",
           group_cells_by="cluster",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           cell_size=if (include.vehicle) 0.5 else 1) +
  ggtitle(glue("UMAP of (N={ncol(filt.cds)}) {paste(cells, sep=', ', collapse=', ')} cells 
               for {paste(treatments, sep=', ', collapse=', ')}
               at doses {paste(doses[!doses %in% c(0)], sep=', ', collapse=', ')}")) +
  scale_color_manual(breaks=rownames(sample.metadata[order(sample.metadata$dose), ]),
                     values=if (include.vehicle) c("grey", "yellow", "orange", "red", "darkred") else c("yellow", "orange", "red", "darkred"))
dev.off()

