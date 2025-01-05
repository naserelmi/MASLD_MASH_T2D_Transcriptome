# Set working directory
setwd('Path to count and metadata files')

# Load count and metadata
count_data <- read.csv('Counts_raw.csv', header=TRUE)
meta_data <- read.csv('Phenotype-Data.csv', header=TRUE)

# Remove first column and set rownames to EntrezGeneID
countdata <- count_data[,-1]
rownames(countdata) <- count_data[,1]
colnames(countdata) <- meta_data$Patient_ID

# Load edgeR library for RNA-seq analysis
library(edgeR)

# Calculate Counts Per Million (CPM) for filtering
myCPM <- cpm(countdata)

# Plot CPM vs raw counts for visual inspection
plot(myCPM[,1], countdata[,1], xlab="CPM", ylab="Raw Count", main=colnames(myCPM)[1], ylim=c(0,50), xlim=c(0,3))
abline(v=0.5)  # Add vertical line at 0.5 CPM
abline(v=1)    # Add vertical line at 1 CPM
abline(h=10)   # Add horizontal line at 10 raw counts

# Filter genes with CPM > 1 in at least 20 samples
isexpr <- rowSums(myCPM > 1) >= 20
summary(isexpr)

# Subset the count data to keep highly expressed genes
counts.keep <- count_data[isexpr,]
dim(countdata)
dim(counts.keep)

write.csv(counts.keep,'raw_filtered_count.csv',row.names = FALSE)
# Load DESeq2 library for differential expression analysis
library(DESeq2)

# Convert metadata columns to factors
meta_data$Diagnosis_compact <- factor(meta_data$Diagnosis_compact)
meta_data$Diabet <- factor(meta_data$Diabet)
meta_data$SEX <- factor(meta_data$SEX)
meta_data$Simplified_class <- factor(meta_data$Simplified_class)

# Differential expression analysis: Diabetic vs Non-diabetic
dds_Diabet <- DESeqDataSetFromMatrix(counts.keep, colData=meta_data, design=~0+Diabet)
dds_Diabet <- DESeq(dds_Diabet)
res_Diabet <- results(dds_Diabet, contrast=c("Diabet", "Diabetic", "Non Diabetic"))
deg_total_diabet <- subset(as.data.frame(res_Diabet), pvalue<0.05 & (log2FoldChange>0.6 | log2FoldChange< -0.6))

# Annotate gene symbols
library(EnsDb.Hsapiens.v86)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_total_diabet), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_total_diabet <- cbind(geneIDs1, deg_total_diabet)
write.csv(deg_total_diabet, "deg_whole_diabet_deseq2.csv")

# Differential expression analysis: Normal+Steatosis vs (Advanced/Non-advanced Fibrosis)
dds_L1 <- DESeqDataSetFromMatrix(counts.keep, colData=meta_data, design=~Diabet+L1)
dds_L1 <- DESeq(dds_L1)
res_L1 <- results(dds_L1, contrast=c("L1", "1", "0"))
deg_L1 <- subset(as.data.frame(res_L1), padj<0.05 & (log2FoldChange>1 | log2FoldChange< -1))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_L1), keytype="GENEID", columns=c("SYMBOL","GENEID"))
library(dplyr)
deg_L1$Ensemble_ID <- rownames(deg_L1)
s1 <- merge(geneIDs1, deg_L1, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s1, "deg_L1_deseq2.csv")

# Differential expression analysis: Advanced Fibrosis vs Non-advanced Fibrosis
L2_index <- (meta_data$Simplified_class=="Advanced_fibrosis" | meta_data$Simplified_class=="Non_advanced_Fibrosis")
samples_L2 <- counts.keep[,L2_index]
meta_data_L2 <- meta_data[L2_index,]
dds_L2 <- DESeqDataSetFromMatrix(samples_L2, colData=meta_data_L2, design=~Diabet+Simplified_class)
dds_L2 <- DESeq(dds_L2)
res_L2 <- results(dds_L2, contrast=c("Simplified_class", "Advanced_fibrosis", "Non_advanced_Fibrosis"))
deg_L2 <- subset(as.data.frame(res_L2), padj<0.05 & (log2FoldChange>1 | log2FoldChange< -1))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_L2), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_L2$Ensemble_ID <- rownames(deg_L2)
s2 <- merge(geneIDs1, deg_L2, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s2, "deg_L2_deseq2.csv")

# Differential expression analysis: Advanced fibrosis with Diabetics vs Healthy
advanced_Diabeic_index <- (meta_data$Simplified_class=="Advanced_fibrosis" & meta_data$Diabet=="Diabetic")
Healthy_NonDiabetic_index <- (meta_data$Diagnosis=="NORMAL" & meta_data$Diabet=="Non Diabetic")
ADV_D_H_index <- advanced_Diabeic_index | Healthy_NonDiabetic_index
samples_ADV_D_H <- counts.keep[,ADV_D_H_index]
meta_data_ADV_D_H <- meta_data[ADV_D_H_index,]
dds_ADV_D_H <- DESeqDataSetFromMatrix(samples_ADV_D_H, colData=meta_data_ADV_D_H, design=~Simplified_class)
dds_ADV_D_H <- DESeq(dds_ADV_D_H)
res_ADV_D_H <- results(dds_ADV_D_H, contrast=c("Simplified_class", "Advanced_fibrosis", "Normal"))
deg_ADV_D_H <- subset(as.data.frame(res_ADV_D_H), padj<0.01 & (log2FoldChange>1 | log2FoldChange< -1))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_ADV_D_H), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_ADV_D_H$Ensemble_ID <- rownames(deg_ADV_D_H)
s3 <- merge(geneIDs1, deg_ADV_D_H, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s3, "deg_ADV_D_H_deseq2.csv")

# Differential expression analysis: Advanced fibrosis vs Diabetic Healthy
Healthy_Diabetic_index <- (meta_data$Diagnosis=="NORMAL" & meta_data$Diabet=="Diabetic")
ADV_D_H_D_index <- advanced_Diabeic_index | Healthy_Diabetic_index
samples_ADV_D_H_D <- counts.keep[,ADV_D_H_D_index]
meta_data_ADV_D_H_D <- meta_data[ADV_D_H_D_index,]
dds_ADV_D_H_D <- DESeqDataSetFromMatrix(samples_ADV_D_H_D, colData=meta_data_ADV_D_H_D, design=~Simplified_class)
dds_ADV_D_H_D <- DESeq(dds_ADV_D_H_D)
res_ADV_D_H_D <- results(dds_ADV_D_H_D, contrast=c("Simplified_class", "Advanced_fibrosis", "Normal"))
deg_ADV_D_H_D <- subset(as.data.frame(res_ADV_D_H_D), padj<0.01 & (log2FoldChange>2 | log2FoldChange< -2))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_ADV_D_H_D), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_ADV_D_H_D$Ensemble_ID <- rownames(deg_ADV_D_H_D)
s4 <- merge(geneIDs1, deg_ADV_D_H_D, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s4, "deg_ADV_D_H_D_deseq2.csv")

# Differential expression analysis: Advanced non-diabetic vs Healthy
ADV_nonDiabeic_index <- (meta_data$Simplified_class=="Advanced_fibrosis" & meta_data$Diabet=="Non Diabetic")
ADV_H_index <- ADV_nonDiabeic_index | Healthy_NonDiabetic_index
samples_ADV_H <- counts.keep[,ADV_H_index]
meta_data_ADV_H <- meta_data[ADV_H_index,]
dds_ADV_H <- DESeqDataSetFromMatrix(samples_ADV_H, colData=meta_data_ADV_H, design=~Simplified_class)
dds_ADV_H <- DESeq(dds_ADV_H)
res_ADV_H <- results(dds_ADV_H, contrast=c("Simplified_class", "Advanced_fibrosis", "Normal"))
deg_ADV_H <- subset(as.data.frame(res_ADV_H), padj<0.01 & (log2FoldChange>2 | log2FoldChange< -2))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_ADV_H), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_ADV_H$Ensemble_ID <- rownames(deg_ADV_H)
s5 <- merge(geneIDs1, deg_ADV_H, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s5, "deg_ADV_H_deseq2.csv")

# Differential expression analysis: Non-advanced vs Healthy
Non_ADV_nonDiabeic_index <- (meta_data$Simplified_class=="Non_advanced_Fibrosis" & meta_data$Diabet=="Non Diabetic")
NonADV_H_index <- Non_ADV_nonDiabeic_index | Healthy_NonDiabetic_index
samples_NonADV_H <- counts.keep[,NonADV_H_index]
meta_data_NonADV_H <- meta_data[NonADV_H_index,]
dds_NonADV_H <- DESeqDataSetFromMatrix(samples_NonADV_H, colData=meta_data_NonADV_H, design=~Simplified_class)
dds_NonADV_H <- DESeq(dds_NonADV_H)
res_NonADV_H <- results(dds_NonADV_H, contrast=c("Simplified_class", "Non_advanced_Fibrosis", "Normal"))
deg_NonADV_H <- subset(as.data.frame(res_NonADV_H), padj<0.01 & (log2FoldChange>2 | log2FoldChange< -2))

# Annotate gene symbols
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(deg_NonADV_H), keytype="GENEID", columns=c("SYMBOL","GENEID"))
deg_NonADV_H$Ensemble_ID <- rownames(deg_NonADV_H)
s6 <- merge(geneIDs1, deg_NonADV_H, by.x="GENEID", by.y="Ensemble_ID", all=TRUE)
write.csv(s6, "deg_NonADV_H_deseq2.csv")


library(limma)

# check your samples grouping
my_meta<-meta_data[,]
my_meta$Diagnosis == group_Diag
my_meta$ICD_250==group_Diabetus
design <- model.matrix(~0+ group_Diabetus+group_Diag)
colnames(design)<- gsub("group","",colnames(design))
colnames(design)<- gsub("_","",colnames(design))
colnames(design)<- gsub("Diag","",colnames(design))
design
# calculate normalization factors between libraries
nf <- calcNormFactors(mycounts)
#
# normalise the read counts with 'voom' function
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
# extract the normalised read counts
counts.voom <- y$E

# Save normalized expression data
write.table(counts.voom, file="Normal.counts.voom.txt", row.names=TRUE, quote=FALSE, sep="\t")