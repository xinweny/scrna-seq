# https://github.com/cole-trapnell-lab/sci-plex

#### Packages ####
library(scran)
library(monocle3)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(gdata)
library(vioplot)
library(glue)

#### Functions ####
dm_analysis <- function(counts) {
  means <- rowMeans(counts)
  cv2 <- apply(counts, 1, var) / (means ^ 2)
  
  dm <- DM(means, cv2)
  
  return(dm)
}

correlation_analysis <- function(counts, dm, top.n) {
  return(rcorr(t(as.matrix(counts[names(head(dm, top.n)), ])),
               type="spearman"))
}

dcorr <- function(corr) {
  sqrt((1 - upperTriangle(corr$r)) / 2)
}

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin")

#### Load data ####
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
protein.coding <- read.csv("./data/GSE139944/data/gencode.v27.transcripts.bed", 
                           sep="\t", quote='', header=FALSE) %>% 
  filter(V9 == "protein_coding") %>%
  .$V7 %>% unique()

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

rownames(filt.counts) <- filt.gene.data$gene_short_name
rownames(filt.gene.data) <- filt.gene.data$gene_short_name

#### DM and correlation analysis ####
dm.corr.list <- list()

# Parameters
cell.type <- "K562"
product <- "Alvespimycin"
d <- c(10, 100, 1000, 10000)

top.n <- 500

vehicle.cells <- filt.col.data[filt.col.data$vehicle, ] %>%
  filter(cell_type == cell.type) %>%
  rownames()
vehicle.counts <- filt.counts[, vehicle.cells]

# RPM normalisation
vehicle.rpm <- vehicle.counts / (colSums(vehicle.counts) / 10^6)

# Get distance-to-median metric
dm.vehicle <- dm_analysis(vehicle.rpm)

dm.df <- data.frame(mean_rpm=rowMeans(vehicle.rpm),
                    dm=dm.vehicle,
                    gene_short_name=names(dm.vehicle),
                    sample=rep("Vehicle_0", length(dm.vehicle)))

# Get correlation metric
# dm.corr.list[["Vehicle_0_topncorr"]] <- correlation_analysis(vehicle.rpm,
#                                                              dm.corr.list[["Vehicle_0_dm"]],
#                                                              top.n)

# Repeat for product
for (i in 1:length(d)) {
  sample.cells <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
    filter(cell_type == cell.type, dose == d[i]) %>% 
    rownames()
  sample.counts <- filt.counts[, sample.cells]
  
  # RPM normalisation
  sample.rpm <- sample.counts / (colSums(sample.counts) / 10^6)
  
  # Get distance-to-median metric
  dm.sample <- dm_analysis(sample.rpm)
    
  dm.df <- rbind(dm.df,
                 data.frame(mean_rpm=rowMeans(sample.rpm),
                            dm=dm.sample,
                            gene_short_name=names(dm.sample),
                            sample=rep(paste0(product, "_", d[i]), length(dm.sample))))
  
  # Get correlation metric
  # dm.corr.list[[paste0(product, "_", d[i], "_topncorr")]] <- correlation_analysis(sample.rpm,
  #                                                                                 dm.df[dm.df$sample == paste0(product, "_", d[i]), "dm"],
  #                                                                                 top.n)
}

# Filter out genes not expressed in any sample
keep.genes <- dm.df %>%
  group_by(gene_short_name) %>%
  dplyr::summarize(mean=mean(mean_rpm)) %>%
  filter(mean > 0) %>%
  .$gene_short_name %>% unique()

dm.df <- dm.df[dm.df$gene_short_name %in% keep.genes, ]

# Set order of samples
dm.df$sample <- factor(dm.df$sample, levels=c("Vehicle_0",
                                              paste0(product, "_10"),
                                              paste0(product, "_100"),
                                              paste0(product, "_1000"),
                                              paste0(product, "_10000")))


#### Visualisation ####
# Plot box plot for DM metric
png(glue("./data/GSE139944/transcriptional_noise/DM/sciPlex3_DM_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
ggplot(dm.df, aes(x=sample, y=dm)) +
  geom_boxplot() +
  xlab("") +
  ylab("Distance-to-median metric (DM)") +
  ggtitle(glue("Distance-to-median metric of genes across {product} treatments
               in {cell.type} cells"))
dev.off()

# Plot DM vs. DM of 2 conditions
dose <- 100

cells.vehicle <- stack(cell.freq)[grep(paste0(cell.type, "_Vehicle_0"), stack(cell.freq)$ind), "values"]
cells.treatment <- stack(cell.freq)[grep(glue("{cell.type}_{product}.+_{dose}$"), stack(cell.freq)$ind), "values"]

dm.dm <- data.frame(dm.df[dm.df$sample == "Vehicle_0", "dm"],
                    dm.df[dm.df$sample == paste0(product, "_", dose), "dm"],
                    dm.df[dm.df$sample == "Vehicle_0", "gene_short_name"])
names(dm.dm) <- c("vehicle_dm", "treatment_dm", "gene_short_name")

dm.dm$direction <- ifelse(dm.dm$treatment_dm > dm.dm$vehicle_dm, "UP", "DOWN")

png(glue("./data/GSE139944/transcriptional_noise/DMvsDM/sciPlex3_DMvsDM_{cell.type}_{product}_{dose}_vs_Vehicle_0.png"),
    width=3000, height=2000, res=300)
dm.dm %>% drop_na(vehicle_dm, treatment_dm) %>%
  ggplot(aes(x=vehicle_dm, y=treatment_dm)) +
  geom_abline(intercept=0, 
              slope=1,
              linetype="dashed", size=0.5, color="gray") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(aes(fill=direction),
             colour="black", pch=21) +
  xlab("DM (Vehicle_0)") + 
  ylab(glue("DM ({product}_{dose})")) +
  labs(title=glue("Distance-to-median metrics (DM) of N={dm.dm %>% drop_na(vehicle_dm, treatment_dm) %>% nrow()} genes in {cell.type} cells"),
       subtitle=glue("{product}_{dose} (nCells={cells.treatment}) vs. Vehicle_0 (nCells={cells.vehicle})
                     nUP={dm.dm %>% drop_na(vehicle_dm, treatment_dm) %>% filter(direction == 'UP') %>% nrow()}, nDOWN={dm.dm %>% drop_na(vehicle_dm, treatment_dm) %>% filter(direction == 'DOWN') %>% nrow()}"),
       fill="Direction (Treatment vs. Vehicle)")
dev.off()

# Plot violin plot for pairwise distance measure for top n genes
# png(glue("./GSE139944/transcriptional_noise/sciPlex3_top{top.n}pairwiseDM_{cell.type}_{product}.png"),
#     width=3000, height=2000, res=300)
# vioplot(dcorr(dm.corr.list[["Vehicle_0_topncorr"]]), 
#         dcorr(dm.corr.list[[paste0(product, "_10_topncorr")]]),
#         dcorr(dm.corr.list[[paste0(product, "_100_topncorr")]]),
#         dcorr(dm.corr.list[[paste0(product, "_1000_topncorr")]]),
#         dcorr(dm.corr.list[[paste0(product, "_10000_topncorr")]]),
#         names=c("Vehicle_0", 
#                 paste0(product, "_10"),
#                 paste0(product, "_100"),
#                 paste0(product, "_1000"),
#                 paste0(product, "_10000")), 
#         col="orange")
# title(main=glue("{cell.type} transcriptional heterogeneity in {product} treatments
#                 for the top N={top.n} highly variable genes"),
#       ylab="Pairwise distance measure")
# dev.off()

