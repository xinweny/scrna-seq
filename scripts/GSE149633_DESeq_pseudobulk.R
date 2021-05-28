#### Packages ####
library(Matrix)
library(Matrix.utils)
library(SingleCellExperiment)
library(DESeq2)
library(dplyr)
library(biomaRt)
library(glue)

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

#### Load data ####
mat <- readMM("./data/GSE149633/data/GSE149633_matrix.mtx")

colnames(mat) <- scan("./data/GSE149633/data/GSE149633_barcodes.tsv", character(), quote="")
rownames(mat) <- scan("./data/GSE149633/data/GSE149633_genes.tsv", character(), quote="")

# Set up cell metadata
coldata <- data.frame(id=colnames(mat),
                      condition=gsub("Mono-|[12]$", "", unlist(lapply(strsplit(colnames(mat), "_"), `[[`, 1))),
                      replicate=paste0("rep", substring(unlist(lapply(strsplit(colnames(mat), "_"), `[[`, 1)),
                                                        nchar(unlist(lapply(strsplit(colnames(mat), "_"), `[[`, 1))))),
                      cell_type=unlist(lapply(strsplit(colnames(mat), "_"), `[[`, 3)))
rownames(coldata) <- coldata$id
coldata$id <- NULL

# Format sample names
coldata$sample <- paste0(coldata$cell_type, "_",
                         coldata$condition, "_",
                         coldata$replicate)

#### Pre-processing ####
sce <- SingleCellExperiment(assays=list(counts=mat),
                            colData=coldata)

# Remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Aggregate across samples
pb <- t(aggregate.Matrix(t(counts(sce)), 
                         groupings=colData(sce)[, c('sample')], 
                         fun="sum"))
pb <- as.matrix(aggregate.Matrix(pb, rownames(pb), fun="sum"))

#### DESeq2 ####
# Set parameters
treatment <- "Tuft_Tm-25h"
control <- "Tuft_DMSO"

alpha <- 0.05

# Create DESeq2 object
deseq.coldata <- data.frame(row.names=colnames(pb),
                            condition=format_condition(colnames(pb)))

dds <- DESeqDataSetFromMatrix(countData=pb,
                              colData=deseq.coldata,
                              design=~condition)

dds$condition <- relevel(dds$condition, ref=control)

# PCA plot
rld <- vst(dds, blind=TRUE)

png("./data/GSE149633/PCA/GSE149633_PCA_aggregated.png")
plotPCA(rld)
dev.off()

# DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", treatment, control),
               alpha=alpha)

deGenes <- as.data.frame(res) %>% arrange(padj, desc(log2FoldChange))

# Add gene symbol
deGenes <- add_ensembl_symbol(deGenes)

# Save results table to output
write.table(deGenes,
            file=glue("./data/GSE149633/deseq/DEGtable/GSE149633_DESeq_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#### Visualisation ####
# Plot MA plot
nUp <- nrow(filter(deGenes, padj < alpha & log2FoldChange > 0))
nDown <- nrow(filter(deGenes, padj < alpha & log2FoldChange < 0))

png(glue("./data/GSE149633/deseq/MAplot/GSE149633_MAplot_{treatment}_vs_{control}.png"))
DESeq2::plotMA(res, 
               main=glue("GSE149633: {treatment} vs. {control}
          N={nUp + nDown}, Nup={nUp}, Ndown={nDown}
          p < {alpha}"))
dev.off()
