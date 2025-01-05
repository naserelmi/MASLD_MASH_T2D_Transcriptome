
setwd('F:/NAFLD/wgcna')
#setwd('D:/Users/Farnaz/Desktop/CBB/NAFLD/wgcna')
dir()
#selection of highly variated genes
#exp_data<-read.csv('orenge.csv',header=TRUE)
exp_data<-read.table('Normal.counts.voom.txt',header=TRUE)
traitData = read.csv("Phenotype-Data.csv");
sample_names<-traitData$Patient_ID

#function to keep rows with maximum n number of NA
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

#---------------
exp_data<-as.data.frame(exp_data)
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
# Take a quick look at what is in the data set:
dim(exp_data);
names(exp_data);
colnames(exp_data)<-c("gene",sample_names)
#exp_data0<-t(exp_data[,2:17397])
rownames(exp_data)<-exp_data[,1]
exp_data<-exp_data[,-1]
exp_data0_m<-as.matrix(exp_data)

#colnames(exp_data0) = exp_data$Patient_ID;
#rownames(exp_data0) = names(exp_data)[2:17397];
#install.packages('matrixStats')
library('matrixStats')

row_mean<-rowMeans(exp_data0_m)
row_var<-rowVars(exp_data0_m)
row_dispersion<-row_var/row_mean
#exp_data0<-exp_data
exp_data0_d<-as.data.frame(exp_data0_m)
exp_data0_d$dispersion<-row_dispersion
exp_data_s<-exp_data0_d[order(-exp_data0_d$dispersion),]
hist(exp_data_s$dispersion,100)
exp_data_selected<-exp_data_s[1:17396,-dim(exp_data_s)[2]]

setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
datExpr0<-t(exp_data_selected)
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

dev.off()
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 200, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
#rownames(datExpr0)<-sample_names
keepSamples=as.logical(rep("TRUE",nrow(datExpr0)))
#keepSamples2= (clust==0 | clust==1 | clust==2 | clust==3)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


traitData = read.csv("Phenotype-Data.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -1];
allTraits = allTraits;
dim(allTraits)
names(allTraits)

orenge_Edge<-cbind(t(exp_data),allTraits)
write.csv(orenge_Edge,"orenge_edgeR.csv",row.names = TRUE,col.names=TRUE)
orenge_Edge_dis<-orenge_Edge[orenge_Edge$Dis_binary==1,]
write.csv(orenge_Edge_dis,"orenge_edgeR_dis.csv",row.names = TRUE,col.names=TRUE)

# Form a data frame analogous to expression data that will hold the clinical traits.

nafldSamples = rownames(datExpr0);
traitRows = match(nafldSamples, allTraits$Patient_ID);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "nafld_voom-01-dataInput.RData")






#section 2 (network reconstruction and module detection)
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "nafld_voom-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


softPower = 3;
#softPower = 10;
adjacency = adjacency(datExpr, power = softPower,type="unsigned");
#adjacency_p_s = adjacency(datExpr[,1:1000], power = softPower,type="unsigned");
#adjacency_s_s = adjacency(datExpr, power = softPower,type="unsigned",corOptions=list(use = 'p', method = 'spearman'));
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 3, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "liver_TOM-blockwise",
                         verbose = 3)
table(bwnet$colors)

bwModuleColors = labels2colors(bwnet$colors)

# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================






#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes

MEList = moduleEigengenes(datExpr, colors = bwModuleColors)

MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


#=====================================================================================


MEDissThres = 0.01
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, bwModuleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, file = "nafld_voom-02-networkConstruction-blockwise.RData")


#section 3

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "nafld_voom-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "nafld_voom-02-networkConstruction-blockwise.RData");
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p",method="spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar =  c(6, 12, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable Diagnosis containing the weight column of datTrait
Diagnosis = as.data.frame(datTraits$Diagnosis);
names(Diagnosis) = "Diagnosis"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Diagnosis, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Diagnosis), sep="");
names(GSPvalue) = paste("p.GS.", names(Diagnosis), sep="");

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "purple"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Disease status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "greenyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Disease status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(Ensemble = colnames(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Diagnosis, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Diagnosis));
geneInfo = geneInfo0[geneOrder, ]

#Ensemble to Gene symbol

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemble_genes<-colnames(datExpr)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensemble_genes,mart= mart)


library(EnsDb.Hsapiens.v86)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensemble_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs1)<-c("gene_symbol","ensembl_gene_id")
final_gene_info<-merge(geneIDs1,geneInfo,by.x="ensembl_gene_id",by.y="Ensemble" ,all=TRUE)

#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(final_gene_info, file = "geneInfo_voom_Diagnosis.csv")

#======================================================================================


# Define variable ICD250 containing the weight column of datTrait
Diagnosis = as.data.frame(datTraits$ICD_250);
names(Diagnosis) = "ICD_250"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Diagnosis, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Diagnosis), sep="");
names(GSPvalue) = paste("p.GS.", names(Diagnosis), sep="");
# Create the starting data frame
geneInfo0 = data.frame(Ensemble = colnames(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Diagnosis, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.ICD_250));
geneInfo = geneInfo0[geneOrder, ]

#Ensemble to Gene symbol

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemble_genes<-colnames(datExpr)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensemble_genes,mart= mart)


library(EnsDb.Hsapiens.v86)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensemble_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs1)<-c("gene_symbol","ensembl_gene_id")
final_gene_info_icd<-merge(geneIDs1,geneInfo,by.x="ensembl_gene_id",by.y="Ensemble" ,all=TRUE)


write.csv(final_gene_info_icd, file = "geneInfo_voom_ICD_250.csv")
