library(xbioc)
library(AnnotationDbi)
# Bulk RNA-seq data
bulk.eset = readRDS("bulk_expression_set.rds")
bulk.mtx = exprs(bulk.eset)
bulk.eset
# clinical conditions
table(bulk.eset$group)
bulk.Diabetic_Advanced.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Diabetic_Advanced']
bulk.Diabetic_Non_Advanced.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Diabetic_Non_Advanced']
bulk.Diabetic_Normal.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Diabetic_Normal']
bulk.Non_Diabetic_Advanced.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Non_Diabetic_Advanced']
bulk.Non_Diabetic_Non_Advanced.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Non_Diabetic_Non_Advanced']
bulk.Non_Diabetic_Normal.mtx = exprs(bulk.eset)[, bulk.eset$group == 'Non_Diabetic_Normal']
# scRNA-seq data
library(SingleCellExperiment)
GSE115469.sce <- readRDS("GSE115469_single_reference.rds")
GSE115469.sce
table(GSE115469.sce$CellType)
# Estimate cell type proportions
library(MuSiC)
library(MuSiC2)
est1 <- music2_prop_t_statistics(bulk.control.mtx = bulk.Non_Diabetic_Normal.mtx, bulk.case.mtx = bulk.Non_Diabetic_Advanced.mtx, select.ct = c('alpha-beta_T_Cells', 'Central_venous_LSECs', 'Cholangiocytes', 'Erythroid_Cells', 'gamma-delta_T_Cells_1', 'gamma-delta_T_Cells_2', 'Hepatic_Stellate_Cells', 'Hepatocyte_1', 'Hepatocyte_2', 'Hepatocyte_3', 'Hepatocyte_4', 'Hepatocyte_5', 'Hepatocyte_6', 'Inflammatory_Macrophage', 'Mature_B_Cells', 'NK-like_Cells', 'Non-inflammatory_Macrophage', 'Periportal_LSECs', 'Plasma_Cells', 'Portal_endothelial_Cells'), sc.sce = GSE115469.sce, clusters = 'CellType', samples = 'sample', n_resample=5, sample_prop=0.25,cutoff_c=0.05,cutoff_r=0.01)
est2 <- music2_prop_t_statistics(bulk.control.mtx = bulk.Non_Diabetic_Normal.mtx, bulk.case.mtx = bulk.Non_Diabetic_Non_Advanced.mtx, select.ct = c('alpha-beta_T_Cells', 'Central_venous_LSECs', 'Cholangiocytes', 'Erythroid_Cells', 'gamma-delta_T_Cells_1', 'gamma-delta_T_Cells_2', 'Hepatic_Stellate_Cells', 'Hepatocyte_1', 'Hepatocyte_2', 'Hepatocyte_3', 'Hepatocyte_4', 'Hepatocyte_5', 'Hepatocyte_6', 'Inflammatory_Macrophage', 'Mature_B_Cells', 'NK-like_Cells', 'Non-inflammatory_Macrophage', 'Periportal_LSECs', 'Plasma_Cells', 'Portal_endothelial_Cells'), sc.sce = GSE115469.sce, clusters = 'CellType', samples = 'sample', n_resample=5, sample_prop=0.25,cutoff_c=0.05,cutoff_r=0.01)
est3 <- music2_prop_t_statistics(bulk.control.mtx = bulk.Diabetic_Normal.mtx, bulk.case.mtx = bulk.Diabetic_Advanced.mtx, select.ct = c('alpha-beta_T_Cells', 'Central_venous_LSECs', 'Cholangiocytes', 'Erythroid_Cells', 'gamma-delta_T_Cells_1', 'gamma-delta_T_Cells_2', 'Hepatic_Stellate_Cells', 'Hepatocyte_1', 'Hepatocyte_2', 'Hepatocyte_3', 'Hepatocyte_4', 'Hepatocyte_5', 'Hepatocyte_6', 'Inflammatory_Macrophage', 'Mature_B_Cells', 'NK-like_Cells', 'Non-inflammatory_Macrophage', 'Periportal_LSECs', 'Plasma_Cells', 'Portal_endothelial_Cells'), sc.sce = GSE115469.sce, clusters = 'CellType', samples = 'sample', n_resample=5, sample_prop=0.25,cutoff_c=0.05,cutoff_r=0.01)
est4 <- music2_prop_t_statistics(bulk.control.mtx = bulk.Diabetic_Normal.mtx, bulk.case.mtx = bulk.Diabetic_Non_advanced.mtx, select.ct = c('alpha-beta_T_Cells', 'Central_venous_LSECs', 'Cholangiocytes', 'Erythroid_Cells', 'gamma-delta_T_Cells_1', 'gamma-delta_T_Cells_2', 'Hepatic_Stellate_Cells', 'Hepatocyte_1', 'Hepatocyte_2', 'Hepatocyte_3', 'Hepatocyte_4', 'Hepatocyte_5', 'Hepatocyte_6', 'Inflammatory_Macrophage', 'Mature_B_Cells', 'NK-like_Cells', 'Non-inflammatory_Macrophage', 'Periportal_LSECs', 'Plasma_Cells', 'Portal_endothelial_Cells'), sc.sce = GSE115469.sce, clusters = 'CellType', samples = 'sample', n_resample=5, sample_prop=0.25,cutoff_c=0.05,cutoff_r=0.01)
est.prop = est1$Est.prop
est$Est.prop
prop_all = cbind('proportion'=c(est.prop), 'sampleID'=rep(rownames(est.prop),times=ncol(est.prop)), 'celltype'=rep(colnames(est.prop), each=nrow(est.prop)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
## you should specify the numbers of each groups by look at prop_all and pData files. for example Non_Diabetic_Normal groups are from 1 to 51
# specify group names
prop_all$group = ifelse(prop_all$sampleID %in% seq(from=1, to=51, by=1), 'Non_Diabetic_Normal', 'Non_Diabetic_Advanced')
prop_all$group = ifelse(prop_all$sampleID %in% seq(from=1, to=91, by=1), 'Non_Diabetic_Normal', 'Non_Diabetic_Non_Advanced')
prop_all$group = ifelse(prop_all$sampleID %in% seq(from=1, to=24, by=1), 'Diabetic_Advanced', 'Diabetic_Normal')
prop_all$group = ifelse(prop_all$sampleID %in% seq(from=1, to=37, by=1), 'Diabetic_Normal', 'Diabetic_Non_Advanced')

cols <-c("alpha-beta_T_Cells" = "aquamarine", "Central_venous_LSECs" = "blue3", "Cholangiocytes" = "chartreuse4", "Erythroid_Cells" = "darkcyan", "gamma-delta_T_Cells_1" = "darkorange1", "gamma-delta_T_Cells_2" = "deepskyblue3", "Hepatic_Stellate_Cells" = "hotpink1", "Hepatocyte_1" = "gold2", "Hepatocyte_2" = "lavender", "Hepatocyte_3" = "lightcyan3", "Hepatocyte_4" = "deepskyblue1", "Hepatocyte_5" = "darkolivegreen3", "Hepatocyte_6" = "magenta2", "Inflammatory_Macrophage" = "olivedrab",  "Non-inflammatory_Macrophage" = "lightslateblue", "Mature_B_Cells" = "palegreen2", "NK-like_Cells" = "salmon2", "Periportal_LSECs" = "thistle", "Plasma_Cells" = "turquoise2", "Portal_endothelial_Cells" = "yellowgreen")
# plot estimated cell type proportions
# jitter plot
prop_plot <- ggplot(prop_all, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,size=0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values=cols)

prop_plot
# save jitter plot
ggsave("Non_Diabetic_Normal VS Non_Diabetic_Non_Advanced_20cell_rep5.tiff", plot = prop_plot, dpi = 500, width = 35, height = 8, units = "in", device = "tiff")
# violin plot
violin_plot <- ggplot(prop_all, aes(x = celltype, y = proportion, fill = celltype)) +
  geom_violin(width = 0.8, alpha = 0.8) +
  ylab('Cell Type Proportions') +
  theme_bw() +
  facet_grid(. ~ group) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  ) +
  scale_fill_manual(values = cols)
violin_plot
# save violin plot
ggsave("Non_Diabetic_Normal VS Non_Diabetic_Non_Advanced_20cell_rep5_violin_plot.tiff", plot = violin_plot, dpi = 500, width = 35, height = 8, units = "in", device = "tiff")
