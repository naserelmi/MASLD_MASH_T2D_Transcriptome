setwd('path_to_WGCNA_results')
gene_info_diagnosis<-read.csv('geneInfo_voom_Diagnosis.csv',header=TRUE)
gene_info_ICD<-read.csv('geneInfo_voom_ICD_250.csv',header = TRUE)
View(gene_info_diagnosis)
mgreenyellow_diag<-gene_info_diagnosis[gene_info_diagnosis$moduleColor=='greenyellow' & gene_info_diagnosis$GS.Diagnosis>=0.5 & gene_info_diagnosis$MM.greenyellow>0.7,]
View(mgreenyellow_diag)

mgreenyellow_ICD<-gene_info_ICD[gene_info_ICD$moduleColor=='greenyellow',]
mgreenyellow_ICD2<-mgreenyellow_ICD[mgreenyellow_diag$ensembl_gene_id %in% mgreenyellow_ICD$ensembl_gene_id,]
View(mgreenyellow_ICD)
mgreenyellow_ICD_filtered<-mgreenyellow_ICD[mgreenyellow_ICD$ensembl_gene_id %in% mgreenyellow_diag$ensembl_gene_id,]
merged_greenyellow<-cbind(mgreenyellow_diag[,1:6],mgreenyellow_ICD_filtered[5:6])
write.csv(merged_greenyellow,"merged_greenyellow_WGCNA.csv")