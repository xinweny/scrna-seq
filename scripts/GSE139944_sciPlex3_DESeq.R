#### Packages ####
library(monocle3)
library(DESeq2)
library(glue)
library(SingleCellExperiment)
library(Matrix.utils)
library(dplyr)

#### Functions ####
add_ensembl_symbol <- function (table) {
  genes <- row.names(table)
  
  if (grepl("ENSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "hsapiens_gene_ensembl"
    symbol <- "hgnc_symbol"
  } else if (grepl("ENSMUSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "mmusculus_gene_ensembl"
    symbol <- "mgi_symbol"
  }
  
  mart <- useDataset(ensemblDataset, useMart("ENSEMBL_MART_ENSEMBL", host="http://www.ensembl.org"))
  geneList <- getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id", symbol),
                    values=genes,
                    mart=mart)
  
  geneList <- distinct(geneList, ensembl_gene_id, .keep_all=TRUE)
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin", 
                "Trichostatin A", "Belinostat", "Mocetinostat")

#### Load data ####
cds <- readRDS("./processed/GSE139944/GSM4150378_sciPlex3_cds_all_cells.RDS")
col.data <- read.csv("./processed/GSE139944/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')

# Format sample names
col.data$cell_product_dose_rep <- paste0(col.data$cell_type, "_",
                                         col.data$product_name, "_", 
                                         col.data$dose, "_",
                                         col.data$replicate)

#### Pre-filtering ####
# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name),])

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose_rep)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose_rep, rownames(cell.freq))]

# Reformat sample name with cell frequency info
# filt.col.data$cell_product_dose_rep <- paste0(filt.col.data$cell_product_dose_rep, " (", filt.col.data$cell_freq, ")")

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(cds)[grep("ENSG", rownames(counts(cds))), 
                           rownames(col.data) %in% rownames(filt.col.data)]

filt.gene.data <- fData(cds)[grep("ENSG", rownames(fData(cds))), ]

rownames(filt.counts) <- filt.gene.data$id
rownames(filt.gene.data) <- filt.gene.data$id

filt.cds <- new_cell_data_set(filt.counts,
                              cell_metadata=filt.col.data[match(colnames(filt.counts), rownames(filt.col.data)), ],
                              gene_metadata=filt.gene.data) 

# Create Single Cell Experiment object
sce <- SingleCellExperiment(assays=list(counts=counts(filt.cds)),
                            colData=filt.col.data)

#### QC and filtering ####
# Remove lowly expressed genes
# sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Aggregate across samples
pb <- t(aggregate.Matrix(t(counts(sce)), 
                         groupings=colData(sce)[, c('cell_product_dose_rep')], 
                         fun="sum"))

rownames(pb) <- gsub("\\.[0-9_A-Z]+$", "", filt.gene.data$id)
pb <- as.matrix(aggregate.Matrix(pb, rownames(pb), fun="sum"))

#### Sample metadata ####
agg.samples <- unique(filt.col.data$cell_product_dose_rep)
sample.metadata <- data.frame(row.names=agg.samples,
                              cell_type=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "cell_type"]),
                              product=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "product_name"]),
                              dose=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "dose"]),
                              rep=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "replicate"]),
                              n_cells=as.vector(cell.freq[match(agg.samples, names(table(filt.col.data$cell_product_dose_rep)))]),
                              target=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "target"]))

sample.metadata[sample.metadata$target == "HSP (e.g. HSP90)", "target"] <- "HSP90"

#### DESeq2 ####
# Select conditions of interest
cell.type <- "MCF7"
treatment <- "Mocetinostat"
doses <- c(10, 100, 1000, 10000)

alpha <- 0.05

# Filter for conditions of interest
for (i in 1:length(doses)) {
  filt.sample.metadata <- with(sample.metadata,
                               sample.metadata[grepl(cell.type, sample.metadata$cell_type) &
                                               grepl(paste0("Vehicle|", treatment), sample.metadata$product), ])
  filt.sample.metadata <- filt.sample.metadata[filt.sample.metadata$dose == 0 |
                                               filt.sample.metadata$dose == doses[i], ]
  
  filt.pb <- as.matrix(pb[, rownames(filt.sample.metadata)])
  
  # Create DESeq2 object
  deseq.coldata <- data.frame(row.names=colnames(filt.pb),
                              condition=c(rep("Vehicle_0", 2), rep(paste0(treatment, "_", doses[i]), 2)))
  
  dds <- DESeqDataSetFromMatrix(countData=filt.pb,
                                colData=deseq.coldata,
                                design=~condition)
  
  dds$condition <- relevel(dds$condition, ref="Vehicle_0")
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  res <- results(dds, 
                 contrast=c("condition", paste0(treatment, "_", doses[i]), "Vehicle_0"),
                 alpha=alpha)
  
  nUp <- nrow(filter(as.data.frame(res), padj < alpha & log2FoldChange > 0))
  nDown <- nrow(filter(as.data.frame(res), padj < alpha & log2FoldChange < 0))
  
  # Plot and save MA plot
  png(glue("processed/GSE139944/deseq/MAplot/sciPlex3_MAplot_{cell.type}_{treatment}_{doses[i]}_vs_Vehicle_0.png"))
  DESeq2::plotMA(res, 
                 main=glue("{cell.type} {treatment}_{doses[i]} vs. Vehicle_0
          n={nUp + nDown}, UP={nUp}, DOWN={nDown}
          p < {alpha}"))
  dev.off()
  
  # Sort and save DESeq2 output table
  res <- add_ensembl_symbol(res)
  de.genes <- as.data.frame(res) %>% arrange(padj, desc(log2FoldChange))
  write.table(de.genes,
              file=glue("processed/GSE139944/deseq/DEGtable/sciPlex3_DESeq_{cell.type}_{treatment}_{doses[i]}_vs_Vehicle_0.txt"),
              row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE) 
}
