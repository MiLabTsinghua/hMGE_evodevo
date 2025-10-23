## Figure 2

#### Figure 2_G

Seurat_test_all = subset(Seurat_test_all, DLF_judge == "Singlet");Seurat_test_all = subset(Seurat_test_all, Revision_CellSubtype %in% c("VZ_RGC", "SVZ_RGC", "IPC"))
Idents(Seurat_test_all) = "Revision_CellSubtype"

Seurat_stereo_all[["RNA"]] = Seurat_stereo_all[["Spatial"]] 
celltype_NV = Mine_Metaneibor(Seurat_test_all, Seurat_stereo_all, "SCSn", "ST", cal_type = 2)

save(celltype_NV, file = "../../results/processed/Figure2_G.RData")

#### Figure 2_I

Seurat_test_HumanRGC = subset(Seurat_test, Species == "Homo sapiens")
Seurat_test_MacaqueRGC = subset(Seurat_test, Species == "Macaca mulatta")
Seurat_test_MouseRGC = subset(Seurat_test, Species == "Mus musculus")

Idents(Seurat_test_HumanRGC) = "Revision_CellSubtype"
Seurat_test_MouseRGC_Petros = subset(Seurat_test_MouseRGC, batch %in% c("Petros Lab Mm_E12.5_MGE", "Petros Lab Mm_E15.5_MGE"));Idents(Seurat_test_MouseRGC_Petros) = "Anno"
S1_DEG <- Seurat_test_HumanRGC %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
S2_DEG <- Seurat_test_MouseRGC_Petros %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
gene_out = intersect(S1_DEG$gene[which(S1_DEG$cluster %in% c("SVZ_RGC"))], S2_DEG$gene)
gene_input = intersect(S1_DEG$gene, S2_DEG$gene) %>% setdiff(gene_out)
celltype_NV = Mine_Metaneibor(Seurat_test_HumanRGC, Seurat_test_MouseRGC_Petros, "humanFetal", "mouse_Petros", cal_type = 4, gene_input = gene_input)
save(celltype_NV, file = "../../results/processed/Figure2_I_Petros.RData")

Idents(Seurat_test_HumanRGC) = "Revision_CellSubtype"
Seurat_test_MouseRGC_Harwell = subset(Seurat_test_MouseRGC, batch %in% c("Harwell Lab Mm_E12.5_MGE", "Harwell Lab Mm_E14.5_MGE"));Idents(Seurat_test_MouseRGC_Harwell) = "Anno"
S1_DEG <- Seurat_test_HumanRGC %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
S2_DEG <- Seurat_test_MouseRGC_Harwell %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
gene_input = (intersect(S1_DEG$gene, S2_DEG$gene))
celltype_NV = Mine_Metaneibor(Seurat_test_HumanRGC, Seurat_test_MouseRGC_Harwell, "humanFetal", "mouse_Harwell", cal_type = 4, gene_input = gene_input)
save(celltype_NV, file = "../../results/processed/Figure2_I_Harwell.RData")

#### Figure 2_J

Idents(Seurat_test_HumanRGC) = "Revision_CellSubtype"
Seurat_test_MacaqueRGC_Alex = subset(Seurat_test_MacaqueRGC, batch %in% c("Alex A. Pollen Lab Macaque_E40_MGE", "Alex A. Pollen Lab Macaque_E50_MGE", "Alex A. Pollen Lab Macaque_E65-2019A_MGE", 
                                                                     "Alex A. Pollen Lab Macaque_E65-2019B_MGE", "Alex A. Pollen Lab Macaque_E80-2019_MGE", "Alex A. Pollen Lab Macaque_E80MGE"));
Idents(Seurat_test_MacaqueRGC_Alex) = "cell_subtype"
S1_DEG <- Seurat_test_HumanRGC %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
S2_DEG <- Seurat_test_MacaqueRGC_Alex %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
gene_input = intersect(S1_DEG$gene, S2_DEG$gene) 
celltype_NV = Mine_Metaneibor(Seurat_test_HumanRGC, Seurat_test_MacaqueRGC_Alex, "humanFetal", "macaque_Alex", cal_type = 4, gene_input = gene_input)
save(celltype_NV, file = "../../results/processed/Figure2_J_Alex.RData")

Idents(Seurat_test_HumanRGC) = "Revision_CellSubtype"
Seurat_test_MacaqueRGC_Nenad = subset(Seurat_test_MacaqueRGC, batch %in% c("Nenad Sestan Lab Macaque_E110_MGE", "Nenad Sestan Lab Macaque_E93_MGE"));Seurat_test_MacaqueRGC_Nenad = subset(Seurat_test_MacaqueRGC_Nenad, cell_subtype %in% c("GE RG HOPX NRG1", "GE RG CRYAB", "GE RG STMN2 SOX5"))
Idents(Seurat_test_MacaqueRGC_Nenad) = "cell_subtype"
S1_DEG <- Seurat_test_HumanRGC %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
S2_DEG <- Seurat_test_MacaqueRGC_Nenad %>% FindAllMarkers(logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
gene_input = intersect(S1_DEG$gene, S2_DEG$gene)
celltype_NV = Mine_Metaneibor(Seurat_test_HumanRGC, Seurat_test_MacaqueRGC_Nenad, "humanFetal", "macaque_Nenad", cal_type = 4, gene_input = gene_input)
save(celltype_NV, file = "../../results/processed/Figure2_J_Nenad.RData")
