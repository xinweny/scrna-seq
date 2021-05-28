#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

#### Load data ####
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

# Keep valid cells
col.data <- col.data[scan("./data/GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Obtain cell frequencies
cell.freq <- table(col.data$cell_product_dose)

# Format into output table
samples <- names(cell.freq)

sample.metadata <- data.frame(sample=samples,
                              cell_type=as.vector(col.data[match(samples, col.data$cell_product_dose), "cell_type"]),
                              product=as.vector(col.data[match(samples, col.data$cell_product_dose), "product_name"]),
                              dose=as.vector(col.data[match(samples, col.data$cell_product_dose), "dose"]),
                              target=as.vector(col.data[match(samples, col.data$cell_product_dose), "target"]),
                              pathway=as.vector(col.data[match(samples, col.data$cell_product_dose), "pathway"]),
                              pathway_level_1=as.vector(col.data[match(samples, col.data$cell_product_dose), "pathway_level_1"]),
                              pathway_level_2=as.vector(col.data[match(samples, col.data$cell_product_dose), "pathway_level_2"]),
                              n_cells=as.vector(cell.freq[match(samples, names(table(col.data$cell_product_dose)))]))

# Save to output
write.table(sample.metadata,
            file="./data/GSE139944/data/sciPlex3_samplemetadata.txt",
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

