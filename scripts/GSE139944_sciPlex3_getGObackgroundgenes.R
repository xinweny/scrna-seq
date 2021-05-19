#### Packages ####
library(DESeq2)
library(glue)
library(dplyr)
library(SingleCellExperiment)
library(Matrix.utils)
library(EDASeq)

setwd("~/mrc/project/scrna-seq")

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
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

#### Load data ####
cds <- readRDS("./GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS")
col.data <- read.csv("GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
gene.lengths <- read.csv("./GSE139944/data/GRCh38_genelengths.txt",
                         sep="\t")

# Format sample names
col.data$cell_product_dose_rep <- paste0(col.data$cell_type, "_",
                                         col.data$product_name, "_", 
                                         col.data$dose, "_",
                                         col.data$replicate)

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

filt.col.data <- col.data[col.data$vehicle, ]

control.counts <- counts(cds)[grep("ENSG", rownames(counts(cds))), 
                              rownames(filt.col.data)]

#### Aggregate matrix ####
# Create Single Cell Experiment object
sce <- SingleCellExperiment(assays=list(counts=control.counts),
                            colData=filt.col.data)

# Aggregate cells by sample
pb <- as.matrix(t(aggregate.Matrix(t(counts(sce)), 
                                   groupings=colData(sce)[, c('cell_product_dose_rep')], 
                                   fun="sum")))

# Format gene ID
rownames(pb) <- gsub("\\.[0-9_A-Z]+$", "", rownames(pb))

# Aggregate genes with the same ID
pb <- as.matrix(aggregate.Matrix(pb, rownames(pb), fun="sum"))

#### DESeq2 ####
# Create DDS object
colData <- data.frame(row.names=colnames(pb),
                      condition=colnames(pb))
dds <- DESeqDataSetFromMatrix(countData=pb,
                              colData=colData,
                              design= ~ 1)

# Add gene length info
mcols(dds)$basepairs <- gene.lengths[match(rownames(dds), gene.lengths$Geneid), "Length"]

# Calculate FPKM
dds.fpkm <- fpkm(dds)

# Save to output
write.table(dds.fpkm, file=glue("./GSE139944/GO/sciPlex3_controlFPKM.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#### Filtering ####
cell.types <- c("K562", "A549", "MCF7")
fpkm.thresh <- 1

for (i in 1:length(cell.types)) {
  dds.fpkm.filt <- dds.fpkm[which(rowMeans(as.data.frame(dds.fpkm) %>% dplyr::select(matches(cell.types[i]))) > fpkm.thresh), ]
  fileConn <- file(glue("./GSE139944/GO/gene_sets/sciPlex3_BackgroundGeneSet_{cell.types[i]}.txt"))
  writeLines(rownames(dds.fpkm.filt), fileConn)
  close(fileConn)
}

