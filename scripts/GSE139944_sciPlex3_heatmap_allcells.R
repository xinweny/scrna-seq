#### Packages ####
library(monocle3)
library(gplots)
library(viridis)
library(glue)
library(ggfortify)

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
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin", 
                "Trichostatin A", "Belinostat", "Mocetinostat")
pscount <- 1 # for log2 scaling

#### Load data ####
cds <- readRDS("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS")
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
proteo.list <- read.csv("./data/GSE139944/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")

#### Pre-filtering ####
# Keep valid cells
col.data <- col.data[scan("./data/GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Get list of CORE proteostasis genes for human/mouse
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name),])

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

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(cds)[grep("ENSG", rownames(counts(cds))), 
                          rownames(filt.col.data)]
filt.cds <- new_cell_data_set(filt.counts,
                              cell_metadata=filt.col.data[match(colnames(filt.counts), rownames(filt.col.data)), ],
                              gene_metadata=fData(cds)[grep("ENSG", rownames(fData(cds))), ]) 

# Extract relevant data from CDS object
row.data <- data.frame(fData(filt.cds))

#### Normalisation ####
# Obtain, normalise and aggregate counts
agg.counts <- aggregate_gene_expression(filt.cds, 
                                        cell_group_df=pData(filt.cds)[, c('cell', 'cell_product_dose_rep')], 
                                        norm="size_only",
                                        scale_agg_values=FALSE)

# Format gene names in count matrix
rownames(agg.counts) <- gsub("\\.[0-9_A-Z]+$", "", rownames(agg.counts))
row.data$id <- gsub("\\.[0-9_A-Z]+$", "", rownames(row.data))

#### Filtering ####
# Filter for genes of interest
agg.filt.counts <- as.data.frame(as.matrix(agg.counts[which(rownames(agg.counts) %in% proteo.genes), ]))

# Format row labels 
rownames(agg.filt.counts) <- row.data[match(rownames(agg.filt.counts), row.data$id), 'gene_short_name']

# Filter out genes with low expression
agg.filt.counts <- agg.filt.counts[rowMeans(agg.filt.counts) > 0.1, ]

#### Heatmap ####
# Sample metadata
agg.samples <- unique(filt.col.data$cell_product_dose_rep)
sample.metadata <- data.frame(row.names=agg.samples,
                              cell_type=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "cell_type"]),
                              product=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "product_name"]),
                              dose=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "dose"]),
                              rep=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "replicate"]),
                              n_cells=as.vector(cell.freq[match(agg.samples, names(table(filt.col.data$cell_product_dose_rep)))]),
                              target=as.vector(filt.col.data[match(agg.samples, filt.col.data$cell_product_dose_rep), "target"]))

sample.metadata[sample.metadata$target == "HSP (e.g. HSP90)", "target"] <- "HSP90"

sample.metadata <- sample.metadata %>%
  mutate(cell_type=factor(cell_type, level=c("K562", "A549", "MCF7")),
         product=factor(product, level=unique(sample.metadata$product)),
         rep=factor(rep, level=c("rep1", "rep2")),
         target=factor(target, level=c("HSP90", "HDAC", "Vehicle"))) %>%
  arrange(cell_type, target, product, dose, rep) %>%
  mutate_if(is.factor, as.character)

# Reorder matrix
agg.filt.counts <- agg.filt.counts[, rownames(sample.metadata)]

# Initialisation of colour labeling
targets <- c("HSP90", "HDAC")
target.cols <- c('red', 'blue') # List of available colours for labeling targets
col.cols <- rep('black', ncol(agg.filt.counts)) # Initialise colours of all columns to black

for (i in 1:length(targets)) {
  col.cols[colnames(agg.filt.counts) %in% filt.col.data[grep(targets[i], filt.col.data$target), 'cell_product_dose_rep']] <- target.cols[i]
}

# Plot heatmap
agg.matrix <- log2(as.matrix(agg.filt.counts) + pscount)

png(file=glue("./data/GSE139944/heatmap/sciPlex3_allcells_HSP90i_HDACi_proteostasis_heatmap.png"), 
    width=18000, height=24000, res=300)
heatmap.2(agg.matrix,
          Rowv=TRUE,
          Colv=FALSE,
          main=glue("Log mean UMI counts per cell of CORE proteostasis genes (n={nrow(agg.matrix)})
                    in K562, A549 and MCF7 cells treated with Vehicle (black), HSP90i (red) or HDACi (blue) (N={ncol(agg.matrix)})"),
          dendrogram="row",
          scale="none",
          col=magma(299),
          breaks=seq(0, 2, length.out=300),
          labCol=colnames(agg.matrix), 
          labRow=rownames(agg.matrix),
          cexRow=1,
          colCol=col.cols,
          srtCol=45,
          trace="none",
          density.info="none",
          key=TRUE, keysize=0.3,
          lhei=c(1, 45), lwid=c(1, 10),
          key.xlab=glue("log10(Mean UMI counts per cell + {pscount})"),
          margins=c(15, 5))
dev.off()

#### PCA plot #### ####
agg.samples <- unique(filt.col.data$cell_product_dose_rep)
agg.metadata <- data.frame(row.names=agg.samples,
                           cell_type=unlist(lapply(strsplit(agg.samples, "_"), `[[`, 1)),
                           product=unlist(lapply(strsplit(agg.samples, "_"), `[[`, 2)),
                           dose=unlist(lapply(strsplit(agg.samples, "_"), `[[`, 3)),
                           rep=gsub("\\s\\([0-9]+\\)", "", unlist(lapply(strsplit(agg.samples, "_"), `[[`, 4))),
                           n_cells=gsub("rep[1-2]\\s\\(|\\)", "", unlist(lapply(strsplit(agg.samples, "_"), `[[`, 4))))

agg.metadata <- agg.metadata[match(rownames(t(agg.matrix)), rownames(agg.metadata)), ]

# Plot PCA
png(file=glue("./data/GSE139944/PCA/sciPlex3_PCA_HSP90i.png"), 
    width=3000, height=2000, res=300)
pca <- prcomp(t(agg.matrix))
autoplot(pca, data=agg.metadata,
         colour="product",
         shape="cell_type",
         size="dose",
         alpha=0.5)
dev.off()


