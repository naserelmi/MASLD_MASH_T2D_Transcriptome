# Load packages
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
# Read the CSV file that contain just ensembl IDs as a data frame
data <- read.delim("ENSG_ID_from_new_raw_data.csv", header = F)
data
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol <- mapIds(org.Hs.eg.db,
                 keys = data$V1,
                 keytype = "ENSEMBL",
                 column = "SYMBOL")

gene_symbol <- cbind(data, symbol)
write.csv(gene_symbol, file = "gene_symbol_new_raw_data")
# check the possible duplication of gene symbols 
# check randomly conversion of gene symbol

# convert to HGNC symbol
symbol <- mapIds(org.Hs.eg.db, attributes = c("ensembl_gene_id", "hgnc_symbol"), keys = data$V1, keytype = "ENSEMBL", column = "SYMBOL")

gene_symbol <- cbind(data, symbol)
write.csv(gene_symbol, file = "HGNC_gene_symbol")



