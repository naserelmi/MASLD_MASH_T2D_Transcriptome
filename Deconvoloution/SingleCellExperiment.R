# create singlecellexperiment reference with GSE115469
cell <- read.delim("GSE115469_CellClusterType.txt", header = T, sep = "\t")
table(cell$CellType)
count <- read.delim("GSE115469_Data.csv", header = T, sep = ",")
# convert log2CPM to CPM
path <- file.path("GSE115469_Data.csv")
data <- read.csv(path, row.names = 1, header = TRUE, sep = ",")
# Load the limma package
library(limma)
# Extract Log2CPM values
log2cpm_values <- data[, -1]  # Assuming the first column contains gene names
small_constant <- 0.001
# Convert Log2CPM to CPM using voom
cpm_values <- 2^(log2cpm_values)
# Convert all ones to zeros because log0 == 1
cpm_values[cpm_values == 1] <- 0
# Optionally, set gene names as row names
rownames(cpm_values) <- log2cpm_data[, 1]
# Ensure gene names are unique
unique_gene_names <- make.names(gene_names, unique = TRUE)
# Replace '0' and '1' in unique_gene_names with more descriptive names if necessary
unique_gene_names[unique_gene_names == "0"] <- "Gene0"
unique_gene_names[unique_gene_names == "1"] <- "Gene1"
# Assign unique gene names as row names
rownames(cpm_values) <- unique_gene_names
# rownames of cell must be equal to colnames of count
all(rownames(cell)==colnames(cpm_values))
# or all(nrow(cell)==ncol(count))
# if FALSE
exp <- file.path("GSE115469_Data.csv")
count <- as.matrix(read.table(exp, header = TRUE, sep = ",", row.names = 1, as.is = TRUE))
all(rownames(cell)==colnames(count))
# Create SingleCellExperiment
library(SingleCellExperiment)
GSE115469_single_reference <- SingleCellExperiment(
  assays = list(counts = count),
  rowData = DataFrame(gene.symbol = rownames(count)),
  colData = DataFrame(sample = cell$CellName, CellType = cell$CellType)
)

saveRDS(GSE115469_single_reference, file = "GSE115469_single_reference.rds")