#### Packages ####
library(monocle3)
library(gplots)
library(viridis)
library(glue)

#### Functions ####
aggregate_cols <- function(counts, col.data) {
  t.counts <- as.data.frame(t(counts))
  
  t.counts.agg <- aggregate(t.counts, 
                            list(col.data$product_dose_rep), 
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
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
pscount <- 1
cell.type <- "K562" # K562, A549, MCF7

#### Load data ####
cds <- readRDS(glue("./processed/GSE139944/GSM4150378_sciPlex3_{cell.type}_24hrs.RDS"))
proteo.list <- read.csv("~/mrc/project/rna-seq/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

# Extract relevant data from S4 Object
counts <- exprs(cds)
col.data <- data.frame(pData(cds))

# Format sample names
col.data$product_dose_rep <- paste0(col.data$product_dose, "_", col.data$replicate)

#### Normalisation ####
# Log normalisation
counts <- normalized_counts(cds, norm_method="log", pseudocount=pscount)

# Format gene names in count matrix
rownames(counts) <- gsub("\\.[0-9_A-Z]+$", "", rownames(counts))

#### Filtering ####
# Get list of CORE proteostasis genes
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for selected inhibitors
inhibitors <- c("HSP90")
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$target),])

# Filter for genes and cells of interest
filt.counts <- as.data.frame(as.matrix(counts[which(rownames(counts) %in% proteo.genes),
                                              which(colnames(counts) %in% filt.col.data$cell)]))

# Collapse matrix by average expression per treatment per replicate
agg.filt.counts <- aggregate_cols(filt.counts, filt.col.data)

#### Visualisation ####
# Initialisation of colour labeling
inh.cols <- c('red', 'blue', 'green')
col.cols <- rep('black', ncol(agg.filt.counts))
inh.col.strs <- c()

# Colour column labels depending on the inhibitor target
inh.col.list <- vector(mode="list", length=length(inhibitors))
names(inh.col.list) <- inhibitors

for (i in 1:length(inhibitors)) {
  col.cols[colnames(agg.filt.counts) %in% filt.col.data[grep(inhibitors[i], filt.col.data$target), 'product_dose_rep']] <- inh.cols[i]
  inh.col.list[[inhibitors[i]]] <- inh.cols[i]
  inh.col.strs[i] <- glue("{inhibitors[i]} ({inh.cols[i]})")
}

# Plot heatmap
png(file=glue("processed/GSE139944/sciPlex3_{cell.type}_proteostasis_heatmap_withinhibitors.png"), 
    width=6000, height=4000, res=300)
heatmap.2(agg.filt.counts,
          Rowv=TRUE,
          Colv=TRUE,
          main=glue("Mean log-normalised counts per cell of CORE proteostasis genes (n={nrow(agg.filt.counts)})
                    in {cell.type} cells treated with Vehicle (black) or inhibitors of
                    {paste(inh.col.strs, sep=', ', collapse=', ')} (N={ncol(agg.filt.counts)})"),
          dendrogram="both",
          scale="none",
          col=magma(299),
          labCol=colnames(agg.filt.counts), labRow=FALSE,
          colCol=col.cols,
          srtCol=45,
          trace="none",
          key=TRUE, 
          density.info="none",
          keysize=0.5,
          margins=c(15, 2),
          key.xlab="Mean log-normalised counts per cell")
dev.off()

