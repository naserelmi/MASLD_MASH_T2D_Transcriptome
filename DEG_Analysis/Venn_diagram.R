
library(purrr)
library(RVenn)
library(ggplot2)

#read data
deg_ADV_D_H<-read.csv("deg_ADV_D_H_deseq2.csv")
deg_ADV_D_H_D<-read.csv("deg_ADV_D_H_D_deseq2.csv")
deg_ADV_H<-read.csv("deg_ADV_H_deseq2.csv")
# Extract GENEID columns
genes_ADV_D_H <- deg_ADV_D_H$GENEID
genes_ADV_D_H_D <- deg_ADV_D_H_D$GENEID
genes_ADV_H <- deg_ADV_H$GENEID




list_ADV<-list(genes_ADV_D_H,genes_ADV_D_H_D,genes_ADV_H)
names(list_ADV)<-c('Advanced_Fibrosis with Diabet vs non-Diabetic Healthy','Advanced_Fibrosis with Diabetic vs Diabetic Healthy', 'Non Diabetic Advanced-Fibrosis vs Non Diabetic Healthy')
ADV_model = Venn(list_ADV)
overlap(ADV_model)
#Pairwise differences

ADV_setdif<-discern_pairs(ADV_model, slice = 1:3)
ADV_diabet<-as.array(unlist(ADV_setdif[[3]]))
a_ADV<-c(ADV_diabet)
b_ADV<-setdiff(a_ADV,rownames(deg_ADV_D_H))

genes<-b_ADV
#library(EnsDb.Hsapiens.v86)
#G_list<- ensembldb::select(EnsDb.Hsapiens.v86, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

library(dplyr)
deg_ADV_D_H_D_venn<-deg_ADV_D_H_D
#deg_ADV_D_H_D_venn$Ensemble_ID<-rownames(deg_ADV_D_H_D)
#s9<-merge(G_list,deg_ADV_D_H_D_venn,by.x="GENEID",by.y="Ensemble_ID",all=TRUE)
#s9<-s9[s9$GENEID %in% genes,]
#write.csv(s9,"deg_ADV_D_H_D_venn_deseq2.csv")
names(list_ADV) <- c('Advanced_Fibrosis with Diabet\nvs non-Diabetic Healthy',
                     'Advanced_Fibrosis with Diabetic\nvs Diabetic Healthy',
                     'Non Diabetic Advanced-Fibrosis\nvs Non Diabetic Healthy')
ggvenn::ggvenn(list_ADV,fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),stroke_size = 0.5, set_name_size = 8,fill_alpha = 0.5,text_size = 7)

dev.off()
###################
#read data
deg_non_ADV_D_H<-read.csv("deg_non_ADV_D_H_deseq2.csv")
deg_non_ADV_D_H_D<-read.csv("deg_non_ADV_D_H_D_deseq2.csv")
deg_non_ADV_H<-read.csv("deg_non_ADV_H_deseq2.csv")
# Extract GENEID columns
genes_Non_ADV_D_H <- deg_non_ADV_D_H$GENEID
genes_non_ADV_D_H_D <- deg_non_ADV_D_H_D$GENEID
genes_non_ADV_H <- deg_non_ADV_H$GENEID


list_non_ADV<-list(genes_Non_ADV_D_H,genes_non_ADV_D_H_D,genes_non_ADV_H)
names(list_non_ADV) <- c('Non_Advanced_Fibrosis with Diabet\nvs non-Diabetic Healthy',
                         'Non_Advanced_Fibrosis with Diabetic\nvs Diabetic Healthy',
                         'Non Diabetic non_Advanced-Fibrosis\nvs Non Diabetic Healthy')
non_ADV_model = Venn(list_non_ADV)
overlap(non_ADV_model)
#Pairwise differences

non_ADV_setdif<-discern_pairs(non_ADV_model, slice = 1:3)

non_ADV_diabet<-as.array(unlist(non_ADV_setdif[[3]]))
a_non_ADV<-c(non_ADV_diabet)
b_non_ADV<-setdiff(a_non_ADV,rownames(deg_non_ADV_D_H))

genes<-b_non_ADV
library(EnsDb.Hsapiens.v86)
G_list<- ensembldb::select(EnsDb.Hsapiens.v86, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

library(dplyr)
deg_non_ADV_D_H_D_venn<-deg_non_ADV_D_H_D
deg_non_ADV_D_H_D_venn$Ensemble_ID<-rownames(deg_non_ADV_D_H_D)
s10<-merge(G_list,deg_non_ADV_D_H_D_venn,by.x="GENEID",by.y="Ensemble_ID",all=TRUE)
s10<-s10[s10$GENEID %in% genes,]

write.csv(s10,"deg_non_ADV_D_H_D_venn_deseq2.csv")
ggvenn::ggvenn(list_non_ADV,fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),stroke_size = 0.5, set_name_size = 8,fill_alpha = 0.5,text_size = 7)

dev.off()

#intersection of advanced and non-advanced
genes_common<-intersect(s9$GENEID,s10$GENEID)
genes_specific_advanced<-setdiff(s9$GENEID,s10$GENEID)

s10_common<-s10[s10$GENEID %in% genes_common,]
s9_common<-s9[s9$GENEID %in% genes_common,]
s9_uniq<-s9[s9$GENEID %in% genes_specific_advanced,]

write.csv(s9_common,'advanced_common_with_nonADV_venn.csv')
write.csv(s10_common,'nonadvanced_common_with_ADV_venn.csv')

write.csv(s9_uniq,'advanced_uniq_venn.csv')

#trim data for AI based on DEG
L1_orange<-read.csv('orenge_Layer1.csv',header = T)
L2_orange<-read.csv('orenge_edgeR_layer2.csv',header = T)
deg_diabet<-read.csv("deg_whole_diabet_qval_deseq2.csv",header = T)
deg_L1<-read.csv("deg_L1_qval_deseq2.csv",header = T)
deg_L2<-read.csv("deg_L2_qval_deseq2.csv",header = T)
aa<-as.array(colnames(L1_orange))
L1_trim_diabet<-L1_orange[,c(deg_diabet$GENEID,aa[17398:17408])]
L1_trim_L1<-L1_orange[,c(deg_L1$GENEID,aa[17398:17408])]
L2_trim_L2<-L2_orange[,c(deg_L2$GENEID,aa[17398:17408])]

write.csv(L1_trim_diabet,"L1_orange_Diabet_Deg.csv")
write.csv(L1_trim_L1,"L1_orange_L1_Deg.csv")
write.csv(L2_trim_L2,"L2_orange_L2_Deg.csv")