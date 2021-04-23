#### Packages ####
library(monocle3)
library(gplots)
library(viridis)
library(glue)

#### Functions ####
aggregate_cols <- function(counts, col.data) {
  t.counts <- as.data.frame(t(counts))
  
  t.counts.agg <- aggregate(t.counts, 
                            list(col.data$cell_product_dose_rep), 
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
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin")
pscount <- 1 # for log2 scaling

#### Load data ####
cds <- readRDS(glue("./processed/GSE139944/GSM4150378_sciPlex3_cds_all_cells.RDS"))
col.data <- read.csv("./processed/GSE139944/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
proteo.list <- read.csv("~/mrc/project/rna-seq/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

#### Pre-filtering ####
# Get list of CORE proteostasis genes for human/mouse
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name),])

# Filter out mouse genes and non-HSP90 targets
cds <- new_cell_data_set(counts(cds)[grep("ENSG", rownames(counts(cds))), 
                                     rownames(col.data) %in% rownames(filt.col.data)],
                         cell_metadata=pData(cds)[rownames(col.data) %in% rownames(filt.col.data), ],
                         gene_metadata=fData(cds)[grep("ENSG", rownames(fData(cds))), ]) 

# Extract relevant data from CDS object
row.data <- data.frame(fData(cds))

#### Normalisation ####
# Format sample names
filt.col.data$cell_product_dose_rep <- paste0(filt.col.data$cell_type, "_",
                                              filt.col.data$product_name, "_", 
                                              filt.col.data$dose, "_",
                                              filt.col.data$replicate)

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose_rep)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose_rep, rownames(cell.freq))]

# Reformat sample name with cell frequency info
filt.col.data$cell_product_dose_rep <- paste0(filt.col.data$cell_product_dose_rep, " (", filt.col.data$cell_freq, ")")

# Obtain, normalise and aggregate counts
counts <- aggregate_gene_expression(cds, 
                                    cell_group_df=filt.col.data[, c('cell', 'cell_product_dose_rep')], 
                                    norm="size_only")

# Format gene names in count matrix
rownames(counts) <- gsub("\\.[0-9_A-Z]+$", "", rownames(counts))
row.data$id <- gsub("\\.[0-9_A-Z]+$", "", rownames(row.data))

#### Filtering ####
# Filter for genes of interest
filt.counts <- as.data.frame(as.matrix(counts[which(rownames(counts) %in% proteo.genes), ]))

# Format row labels 
rownames(filt.counts) <- row.data[match(rownames(filt.counts), row.data$id), 'gene_short_name']

#### Visualisation ####
# Initialisation of colour labeling
inh.cols <- c('red', 'blue', 'green') # List of available colours for labeling inhibitors
col.cols <- rep('black', ncol(filt.counts)) # Initialise colours of all columns to black
inh.col.strs <- c() # Initialise vector for labeling in title

# Colour column labels depending on the inhibitor and set title labeling
inh.col.list <- vector(mode="list", length=length(inhibitors))
names(inh.col.list) <- inhibitors

for (i in 1:length(inhibitors)) {
  col.cols[colnames(filt.counts) %in% filt.col.data[grep(inhibitors[i], filt.col.data$product_name), 'cell_product_dose_rep']] <- inh.cols[i]
  inh.col.list[[inhibitors[i]]] <- inh.cols[i]
  inh.col.strs[i] <- glue("{inhibitors[i]} ({inh.cols[i]})")
}

# Plot heatmap
png(file=glue("processed/GSE139944/sciPlex3_allcells_HSP90i_proteostasis_heatmap.png"), 
    width=12000, height=24000, res=300)
heatmap.2(log2(as.matrix(filt.counts) + pscount),
          Rowv=TRUE,
          Colv=FALSE,
          main=glue("Log mean UMI counts per cell of CORE proteostasis genes (n={nrow(filt.counts)})
                    in K562, A549 and MCF7 cells treated with Vehicle (black) or HSP90i inhibitors
                    {paste(inh.col.strs, sep=', ', collapse=', ')} (N={ncol(filt.counts)})"),
          dendrogram="row",
          scale="none",
          col=magma(299),
          labCol=colnames(filt.counts), labRow=rownames(filt.counts),
          colCol=col.cols,
          srtCol=45,
          trace="none",
          density.info="none",
          key=TRUE, keysize=0.3,
          lhei=c(1, 45), lwid=c(1, 10),
          key.xlab=glue("log2(Mean UMI counts per cell + {pscount})"),
          margins=c(15, 5))
dev.off()

