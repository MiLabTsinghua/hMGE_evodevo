## Figure 4

#### Figure 4_A

DSP.seurat <- AddModuleScore(DSP.seurat,
                             features = Lineage_Score_list,
                             name = "Lineage_Score")
colnames(DSP.seurat@meta.data)[grep("Lineage_Score", colnames(DSP.seurat@meta.data))] = names(Lineage_Score_list)
DSP.seurat = subset(DSP.seurat, celltype == "ROI_DCX")
scatterpie_df = DSP.seurat@meta.data[, c("ROICoordinateX", "ROICoordinateY", names(Lineage_Score_list))]

write.xlsx(Lineage_Score_list, file = "../../results/tables/Figure4_H_ScoreGene.xlsx")

#### Figure 4_F

Express_average = (AverageExpression(ref, assays = "SCT", 
                                     features = union(gene_selected, rownames(ref@assays$SCT@data)), 
                                     group.by = "Revision_Lineage_New"))[["SCT"]]
df = as.data.frame(Express_average)
df$GeneType = "Common";df$GeneType[which(rownames(df) %in% gene_selected)] = "Variable"
df$SumUp = df$CRABP1ANGPT2 + df$EPHA5MEF2C + df$LHX6NFIA;df[, c(1,2,3)] = df[, c(1,2,3)]/df$SumUp
pannel_loosed = rownames(df)[which(MatrixGenerics::rowMins(df[, c(1,2,3)] %>% as.matrix()) <= 0.2)]
df$gene = rownames(df);df$lab = NA;df$lab[which(df$GeneType == "Variable")] = df$gene[which(df$GeneType == "Variable")]
df$GeneType[which(df$gene %in% pannel_loosed)] = "Variable"

save(df, file = "../../results/Figure4_F.RData")

#### Figure 4_G

df <- melt(scatterpie_df[, c(3:7)], id.vars=c("Node", "type"),
           measure.vars = c("CRABP1ANGPT2", "LHX6NFIA", "EPHA5MEF2C"), 
           # variable.name = "type", 
           # value.name = "x"
)
df$variable = factor(df$variable, levels = c("EPHA5MEF2C", "LHX6NFIA", "CRABP1ANGPT2"));df$type = as.factor(df$type)

save(df, file = "../../results/Figure4_G.RData")

## Figure 5

#### Figure 5_E

setname = "DynamicGene"
Seurat_test_all <- AddModuleScore(Seurat_test_all,
                                  features = gene_order,
                                  name = c("EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score"))
colnames(Seurat_test_all@meta.data)[grep(pattern = "EPHA5MEF2CScore|LHX6NFIAScore|CRABP1ANGPT2Score", x = colnames(Seurat_test_all@meta.data))] = c("EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score")

Seurat_test_all <- AddModuleScore(Seurat_test_all,
                                  features = gene_order_2,
                                  name = c("SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score"))
colnames(Seurat_test_all@meta.data)[grep(pattern = "SVZFateScore|NR2F1NR2F2Score|LHX8ISL1Score", x = colnames(Seurat_test_all@meta.data))] = c("SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score")

Seurat_test_all$test_cluster = Seurat_test_all$Revision_LineageType_New
df = Seurat_test_all@meta.data[, c("test_cluster", "EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score", "SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score")]

