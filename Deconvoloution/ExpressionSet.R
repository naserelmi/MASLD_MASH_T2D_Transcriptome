library(Biobase)
# Assay data
exp <- file.path("Counts.csv")
exprs <- as.matrix(read.table(exp, header = TRUE, sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE))
class(exprs)
dim(exprs)
colnames(exprs)
# Phenotypic data
pDataFile <- file.path("PhenoData-Data.csv")
pData <- read.csv(pDataFile, header=TRUE, sep=",")
pData
dim(pData)
rownames(pData)
summary(pData)
# the number of rows of phenotypic data should match the number of columns of expression data
all(rownames(pData)==colnames(exprs))
# if FALSE, run these command
sample_ID <- pData$sample_ID
rownames(pData) <- sample_ID
all(rownames(pData)==colnames(exprs))
# specify group names (in your file)
metadata <- data.frame(labelDescription=c("sample_ID", "Diabetic_Advanced/Diabetic_Non_advanced/Diabetic_Normal/Non_Diabetic_Advanced/Non_Diabetic_Non_advanced/Non_Diabetic_Normal"), row.names= c("sample_ID", "group"))
# phenoData
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
head(pData(phenoData))
# An ExpressionSet object is created by assembling its component parts and calling the ExpressionSet constructor:
ExpressionSet.ob <- ExpressionSet(assayData=exprs, phenoData=phenoData)
saveRDS(ExpressionSet.ob, file = "bulk_expression_set.rds")
