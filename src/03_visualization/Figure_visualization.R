source("../utils/Required_packages.R")
source("../utils/Plot_parameters.R")
source("../utils/Function_tools.R")

setwd("../../results/figures")

#### Figure 1
i = 3

fig = DimPlot(subset(Seurat_test, DLF_judge == "Singlet"), group.by = c("Revision_CellClass"), reduction = umap_selected[i], shuffle = T, label = F, cols = MGE_cluster_color
              ) +
  # ggsci::scale_color_npg() +
  ggtitle("Human MGE GW9~GW39", "Cell Class") +
  theme_void();fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure1_B", nwid = 14, nhei = 12)

fig = DimPlot(subset(Seurat_test, DLF_judge == "Singlet"), group.by = c("batch"), reduction = umap_selected[i], shuffle = T, label = F, cols = Batch_color
) +
  # ggsci::scale_color_npg() +
  ggtitle("Human MGE GW9~GW39", "Batch") +
  theme_void();fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure1_B2", nwid = 14, nhei = 12)

Seurat_test_plot = subset(Seurat_test, DLF_judge == "Singlet")
fig = AverageHeatmap(object = Seurat_test_plot, assays = "RNA", slot = "data", 
                     markerGene = gene_selected, annoCol = TRUE, htCol = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 0.5)(3),
                     myanCol = MGE_cluster_color[levels(Idents(Seurat_test))],
                     width = 12, height = 20,
                     showRowNames = F, markGenes = Anno);fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure1_C", nwid = 12, nhei = 20, dpi_in = 320)

mat_out = fig@matrix %>% as.data.frame.matrix()
mat_out$Gene = rownames(mat_out)
write.xlsx(mat_out, file = "../../results/tables/Figure1_C.xlsx",)

#### Figure 2
fig = AverageHeatmap(object = Seurat_test, assays = "SCT", slot = "data", column_names_rot = 45,
                     markerGene = gene_selected, annoCol = TRUE, htCol = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 0.7)(4), htRange = c(-2, -0.3, 0, 2),
                     myanCol = MGE_progenitor_color[levels(Idents(Seurat_test))],
                     width = 12, height = 20,
                     showRowNames = F, markGenes = unlist(Anno));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2_C", nwid = 12, nhei = 20, dpi_in = 320)
mat_out = fig@matrix %>% as.data.frame.matrix()
mat_out$Gene = rownames(mat_out)
write.xlsx(mat_out, file = "../../results/tables/Figure2_C.xlsx",)

Seurat_test_sub = subset(Seurat_test, Revision_CellClass == "Progenitor")
Idents(Seurat_test_sub) = factor(Seurat_test_sub$Revision_CellSubtype, levels = c("VZ_RGC", "SVZ_RGC", "IPC"));new_idents = c("RGC_1", "RGC_2", "IPC");names(new_idents) = c("VZ_RGC", "SVZ_RGC", "IPC");Seurat_test_sub = RenameIdents(Seurat_test_sub, new_idents)
fig = Seurat_Bar(Seurat_test_sub, Seurat_test_sub$batch, color_used = MGE_cluster_color, bar_width = 0.618);fig
Single_ggplot_light(fig, "../../results/figures", if_egg = F,
                    ggfigtag = "Figure2_D", nwid = 24, nhei = 14)

#### Figure 2_F
fig = AverageHeatmap(object = Seurat_stereo_sum, 
                     assays = "Spatial", slot = "data", htCol = colorRampPalette(bias = 0.8, colors = c('#176BA0','black','#EE9A3A'))(5), 
                     htRange = c(-2,-1,0,1,2),
                     markerGene = gene_selected, annoCol = TRUE, 
                     myanCol = MGE_cluster_color[levels(Idents(Seurat_stereo_sum))],
                     width = 12, height = 20,
                     showRowNames = T, markGenes = unlist(Anno));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2_F", nwid = 12, nhei = 20, dpi_in = 320)
mat_out = fig@matrix %>% as.data.frame.matrix()
mat_out$Gene = rownames(mat_out)
write.xlsx(mat_out, file = "../../results/tables/Figure2_F.xlsx",)
