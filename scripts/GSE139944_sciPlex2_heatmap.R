#### Packages ####
library(monocle3)
library(gplots)
library(viridis)
library(glue)

#### Functions ####
aggregate_cols <- function(counts, col.data) {
  t.counts <- as.data.frame(t(counts))
  
  t.counts.agg <- aggregate(t.counts, 
                            list(col.data$top_oligo), 
                            mean)
  
  counts.agg <- t(t.counts.agg)
  colnames(counts.agg) <- counts.agg[1, ]
  counts.agg <- as.data.frame(counts.agg[-1, ])
  
  counts.agg.mat <- as.matrix(sapply(counts.agg, as.numeric))
  rownames(counts.agg.mat) <- rownames(counts.agg)
  
  return(counts.agg.mat)
}

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
pscount <- 1

#### Load data ####
cds <- readRDS("./GSE139944/data/GSM4150377_sciPlex2_cds.RDS")
proteo.list <- read.csv("./GSE139944/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

# Extract relevant data from S4 Object
counts <- exprs(cds)
col.data <- data.frame(pData(cds))

#### Normalisation ####
# Normalise and aggregate
counts <- aggregate_gene_expression(cds, 
                                    cell_group_df=col.data[, c('Cell', 'top_oligo')], 
                                    norm="size_only")

# Format gene names in count matrix
rownames(counts) <- gsub("\\.[0-9_A-Z]+$", "", rownames(counts))

#### Filtering ####
# Get list of CORE proteostasis genes
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for selected inhibitors
inhibitors <- c("SAHA", "BMS")
filt.col.data <- col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$top_oligo), ]

# Filter for genes and cells of interest
filt.counts <- as.data.frame(as.matrix(counts[which(rownames(counts) %in% proteo.genes),
                             which(colnames(counts) %in% filt.col.data$top_oligo)]))

#### Visualisation ####
# Initialisation of colour labeling
inh.cols <- c('red', 'blue', 'green')
col.cols <- rep('black', ncol(filt.counts))
inh.col.strs <- c()

# Colour column labels depending on the inhibitor target
inh.col.list <- vector(mode="list", length=length(inhibitors))
names(inh.col.list) <- inhibitors

for (i in 1:length(inhibitors)) {
  col.cols[colnames(filt.counts) %in% filt.col.data[grep(inhibitors[i], filt.col.data$top_oligo), 'top_oligo']] <- inh.cols[i]
  inh.col.list[[inhibitors[i]]] <- inh.cols[i]
  inh.col.strs[i] <- glue("{inhibitors[i]} ({inh.cols[i]})")
}


# Plot heatmap
png(file="./GSE139944/heatmap/sciPlex2_proteostasis_heatmap_inhibitors.png", 
    width=6000, height=4000, res=300)
heatmap.2(log2(as.matrix(filt.counts) + pscount),
          Rowv=TRUE,
          Colv=TRUE,
          main=glue("Log mean UMI counts per cell of CORE proteostasis genes (n={nrow(filt.counts)}) in A549 cells
                    treated with different doses of {paste(inh.col.strs, sep=', ', collapse=', ')}
                    across scRNA-seq wells (N={ncol(filt.counts)})"),
          dendrogram="both",
          scale="none",
          col=magma(299),
          labCol=colnames(filt.counts), labRow=FALSE,
          colCol=col.cols,
          srtCol=45,
          trace="none",
          key=TRUE, 
          density.info="none",
          keysize=0.5,
          margins=c(5, 2),
          key.xlab=glue("log2(Mean UMI counts per cell + {pscount})"))
dev.off()

