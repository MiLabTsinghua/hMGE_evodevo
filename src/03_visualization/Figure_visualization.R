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

#### Figure 2_G

fig = pheatmap::pheatmap(celltype_NV,
                         cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA, #cluster_rows = F, cluster_cols = F,
                         color = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 0.6)(1000));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2_G", nwid = 20, nhei = 20, dpi_in = 320)

#### Figure 2_H

Age_loop = "";load("") # Any specific StereoSeq dataset
pixel_x = max(Seurat_stereo$x) - min(Seurat_stereo$x);pixel_x
pixel_y = max(Seurat_stereo$y) - min(Seurat_stereo$y);pixel_y

fig = DimPlot(Seurat_stereo, reduction = "RealLocation_reduction",
              group.by = "plot_type", cols = MGE_cluster_color,
              pt.size = ptz) + 
  # scale_color_manual(values = rev(c(color_selected, "grey100"))) + 
  ggtitle(paste0(Age_loop, "_MGE")) + 
  theme_void() + NoLegend();fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure2_H_", Age_loop, "_MGE"),
                    nwid = round(pixel_x/500), nhei = round(pixel_y/500), if_egg = T)

#### Figure 2_I

fig = pheatmap::pheatmap(celltype_NV,
                          cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA, cluster_rows = F, cluster_cols = F,
                          color = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 1.8)(1000));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2I_Petros", nwid = 20, nhei = 20, dpi_in = 320)

fig = pheatmap::pheatmap(celltype_NV,
                         cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA, #cluster_rows = F, cluster_cols = F,
                         color = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 1.5)(1000));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2I_Harwell", nwid = 20, nhei = 20, dpi_in = 320)

#### Figure 2_J

fig = pheatmap::pheatmap(celltype_NV,
                         cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA, cluster_rows = T, cluster_cols = T,
                         color = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 1.8)(1000));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2J_Alex", nwid = 20, nhei = 20, dpi_in = 320)

fig = pheatmap::pheatmap(celltype_NV,
                         cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA, cluster_rows = F, cluster_cols = F,
                         color = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 1.8)(1000));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure2J_Nenad", nwid = 20, nhei = 20, dpi_in = 320)

## Figure 4

#### Figure 4_A

fig = AverageHeatmap(object = Seurat_test_plot, assays = "RNA", slot = "data", column_names_rot = 45,
                     markerGene = gene_selected, annoCol = TRUE, htCol = colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), bias = 1)(3), htRange = c(-2, 0, 2),
                     myanCol = MGE_branch_color[levels(Idents(Seurat_test))],
                     width = 12, height = 20,
                     showRowNames = F, markGenes = unlist(Anno));fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure4_A", nwid = 12, nhei = 20, dpi_in = 320)
mat_out = fig@matrix %>% as.data.frame.matrix()
mat_out$Gene = rownames(mat_out)
write.xlsx(mat_out, file = "../../results/tables/Figure4_A.xlsx",)

#### Figure 4_B

Location_Color[c("Cortical GABAergic inhibitory neurons", "Subpallial inhibitory neurons")] = c("grey85", "grey70")
fig <- ggplot(umap, aes(x= umap_1 , y = umap_2 , color = Resource)) +  
  geom_point(size = 0.7 , alpha = 0.7, stroke = 0) + 
  scale_color_manual(values = c(Location_Color, MGE_branch_color)) +
  ggtitle("Fetal Neuron Mapping to Adult Neuron") +
  xlim(c(-18,18)) + ylim(c(-15,15)) +  
  # theme(panel.grid.major = element_blank(), #主网格线
  #       panel.grid.minor = element_blank(), #次网格线
  #       panel.border = element_blank(), #边框
  #       axis.title = element_blank(),  #轴标题
  #       axis.text = element_blank(), # 文本
  #       axis.ticks = element_blank(),
  #       panel.background = element_rect(fill = 'white'), #背景色
  #       plot.background=element_rect(fill="white")) + 
  theme_void();fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure4_B1"),
                    nwid = 14, nhei = 13, if_egg = T)

fig <- ggplot(subset(umap, Resource %in% c("Subpallial inhibitory neurons", "Cortical GABAergic inhibitory neurons")), aes(x = umap_1 , y = umap_2 , color = Resource)) +  
  geom_point(size = 0.7 , alpha = 0.7, stroke = 0) + 
  scale_color_manual(values = Location_Color) +
  ggtitle("Adult Neuron (Reference)") +
  xlim(c(-18,18)) + ylim(c(-15,15)) +  
  # theme(panel.grid.major = element_blank(), #主网格线
  #       panel.grid.minor = element_blank(), #次网格线
  #       panel.border = element_blank(), #边框
  #       axis.title = element_blank(),  #轴标题
  #       axis.text = element_blank(), # 文本
  #       axis.ticks = element_blank(),
  #       panel.background = element_rect(fill = 'white'), #背景色
  #       plot.background=element_rect(fill="white")) + 
  theme_void();fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure4_B2"),
                    nwid = 14, nhei = 13, if_egg = T)

fig <- ggplot(subset(umap, !Resource %in% c("Subpallial inhibitory neurons", "Cortical GABAergic inhibitory neurons")), aes(x= umap_1 , y = umap_2 , color = Resource)) +  
  geom_point(size = 0.7 , alpha = 0.9, stroke = 0) + 
  scale_color_manual(values = MGE_branch_color) +
  ggtitle("Fetal Neuron") +
  xlim(c(-18,18)) + ylim(c(-15,15)) +  
  # theme(panel.grid.major = element_blank(), #主网格线
  #       panel.grid.minor = element_blank(), #次网格线
  #       panel.border = element_blank(), #边框
  #       axis.title = element_blank(),  #轴标题
  #       axis.text = element_blank(), # 文本
  #       axis.ticks = element_blank(),
  #       panel.background = element_rect(fill = 'white'), #背景色
  #       plot.background=element_rect(fill="white")) + 
  theme_void();fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure4_B3"),
                    nwid = 14, nhei = 13, if_egg = T)

#### Figure 4_D

fig = Seurat_Bar(Seurat_test, Seurat_test$batch, color_used = MGE_branch_color, bar_width = 0.618);fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure4_D", nwid = 22, nhei = 14)

#### Figure 4_F

fig <- ggtern(data = df, aes(CRABP1ANGPT2, LHX6NFIA, EPHA5MEF2C)) +
  # stat_density_tern(geom='polygon', bdl.val = NA,
  #   aes(fill=..level.., alpha = ..level..)) +
  # scale_fill_gradient(low = "blue", high = "red") +
  # geom_point(aes(color = GeneType), size = 1.5) + 
  geom_point(aes(color = GeneType), size = 1.5) + 
  geom_text(aes(label = lab),check_overlap = F, vjust=1) +
  scale_color_manual(values = c("red", "gray"#, "green"
  )) +
  theme(legend.position = "bottom", legend.key = element_rect(fill = NA)) +
  theme_rgbw();fig
Normal_ggplot_light(fig, f_path = "../../results/figures",
                    normtag = "Figure4_F_ggtern_label", nwid = 10, nhei = 9, dpi_in = 320)

#### Figure 4_H

pixel_x = max(DSP.seurat$ROICoordinateX) - min(DSP.seurat$ROICoordinateX);pixel_x
pixel_y = max(DSP.seurat$ROICoordinateY) - min(DSP.seurat$ROICoordinateY);pixel_y
fig = ggplot() + scatterpie::geom_scatterpie(aes(x=ROICoordinateX, y=ROICoordinateY, r = 200), 
                                             data=scatterpie_df, 
                                             cols=colnames(scatterpie_df)[3:5], color = "white") +
  theme_void()+
  scale_fill_manual(values = MGE_branch_color);fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure4_H", "_SCGW1819MGE_on_DSPGW18MGE"), 
                    nwid = round(pixel_x/300), nhei = round(pixel_y/300), if_egg = T)

#### Figure 4_G

fig = ggplot(data = df, aes(x = variable, y = value, group = Node, color = type))+
  geom_line() +
  geom_point() +
  scale_color_manual(values = (c("#FF0000FF", "#800080"))) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red') +
  # geom_hline(yintercept = 0.15, linetype = 'dotted', col = 'red') +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("GW18 DSP DCX+", "Lineage Score") + labs(x = "", y = "Score") +
  NoLegend();fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure4_G", ""), 
                    nwid = 11, nhei = 18, if_egg = F)

## Figure 5

#### Figure 5_A

i = 3
fig = DimPlot(Seurat_test, group.by = c("Revision_NeuronLineageType"), reduction = umap_selected[i], shuffle = T, label = T, cols = MGE_branch_color) +
  # ggsci::scale_color_npg() +
  ggtitle("Human MGE GW9~GW39", "Neuron Type");fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure5_A1", nwid = 14, nhei = 12)

fig = DimPlot(Seurat_test, group.by = c("predicted.LineageType"), reduction = "ref.umap", shuffle = T, label = F, cols = MGE_branch_color) +
  # ggsci::scale_color_npg() +
  ggtitle("Human MGE GW9~GW39", "IPC Lineage Type");fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure5_A2", nwid = 14, nhei = 12)

fig = DimPlot(subset(Seurat_test, cells = cells_1, invert = T), group.by = c("Revision_LineageType_New"), reduction = "ref.umap", shuffle = T, label = F, cols = MGE_branch_color, pt.size = 1.6) +
  # ggsci::scale_color_npg() +
  ggtitle("Human MGE GW9~GW39", "RGC Lineage Type") #+ 
  # theme_void();
fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure5_A3", nwid = 14, nhei = 12)

#### Figure 5_B

fig = Seurat_Bar(Seurat_test_plot, Seurat_test_plot$batch, color_used = MGE_branch_color, bar_width = 0.618);fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure5_B_Version2", nwid = 22, nhei = 14)

#### Figure 5_D

Idents(Seurat_test) = factor(Seurat_test$predicted.LineageType, levels = rev(c(
  "IPC_NR2F1NR2F2", "IPC_LHX8ISL1",
  "IPC_CRABP1ANGPT2", "IPC_LHX6NFIA", "IPC_EPHA5MEF2C"
)));Seurat_test$batch = factor(Seurat_test$batch, levels = c("GW9_XiaoqunW", "GW10_MGE", "GW13_MGE", "GW13_MGE_rep2", "GW18_MGE", "GW19_MGE", "GW26_MGE", "GW39_MGE"))
fig = Seurat_Bar(Seurat_test, Seurat_test$batch, color_used = MGE_branch_color, bar_width = 0.618);fig
Single_ggplot_light(fig, "../../results/figures", if_egg = T,
                    ggfigtag = "Figure5_D", nwid = 22, nhei = 14)

#### Figure 5_E

fig = DimPlot(Seurat_test, reduction = "Pseudotime_reduction", group.by = "Revision_CellSubtype", 
              cols = MGE_cluster_color, pt.size = 1.4
) + theme_void();fig
Single_ggplot_light(fig, "../../results/figures",if_egg = T,
                    ggfigtag = "Figure5_E1", nwid = 12, nhei = 17)

fig = DimPlot(Seurat_test, reduction = "Pseudotime_reduction", group.by = "Revision_Lineage_New", 
              cols = MGE_branch_color, pt.size = 1.4
) + theme_void();fig
Single_ggplot_light(fig, "../../results/figures",if_egg = T,
                    ggfigtag = "Figure5_E2", nwid = 12, nhei = 17)

fig = DimPlot(Seurat_test, reduction = "Pseudotime_reduction", group.by = "batch", 
              cols = Batch_color, pt.size = 1.4, 
) + theme_void();fig
Single_ggplot_light(fig, "../../results/figures",if_egg = T,
                    ggfigtag = "Figure5_E3", nwid = 12, nhei = 17)

#### Figure 5_F

split_df <- split(df, df$max_type)
score_types <- c("EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score", "NR2F1NR2F2Score", "LHX8ISL1Score");colors <- c("#FF0000FF", "#0000FFFF", "#00EDEDFF", "#009404FF", "#FFD600FF")
fig <- ggplot()
for (i in seq_along(score_types)) {
  layer_data <- split_df[[score_types[i]]]
  
  if (nrow(layer_data) > 0) {
    # 添加新的颜色标度
    if (i > 1) fig <- fig + new_scale_color()
    
    # 添加点图层
    fig <- fig +
      geom_point(
        data = layer_data,
        aes(x = x, y = y, color = max_value),
        size = 1.5,  # 点大小
        alpha = 0.7,  # 透明度
        shape = 19  # 实心圆点
      ) +
      scale_color_gradient(
        name = score_types[i],
        na.value = "white",
        low = "grey95",  # 低值颜色
        high = colors[i],  # 高值颜色
        limits = c(0, 1)  # 颜色范围
      )
  }
}
# 添加通用参数
fig <- fig +
  labs(title = "Points Colored by Highest Score Type",
       x = "X Coordinate",
       y = "Pseudotime Order") +
  # theme_minimal(base_size = 12) +  # 主题和字体大小
  theme_void(base_size = 12) +  # 主题和字体大小
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    # legend.key.height = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1);fig  # 保持坐标轴比例
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure5_F_"),
                    nwid = 12, nhei = 17, if_egg = T)

#### Figure 5_G

### the matrix for the scaled expression 
exp_mat <- df_mat %>%   
  dplyr::select(-pct.exp, -avg.exp) %>%    
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%   
  as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
### the matrix for the percentage of cells express a gene
percent_mat <- df_mat %>%   
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%    
  pivot_wider(names_from = id, values_from = pct.exp) %>%   
  as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
cell_fun = function(j, i, x, y, w, h, fill){  
  grid.rect(x = x, y = y, width = w, height = h,             
            gp = gpar(col = NA, fill = NA))  
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),              
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
# col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
col_fun = circlize::colorRamp2(c(-1, 0, 2), c('#176BA0','black','#EE9A3A'))
# expand_index = 1.25
expand_index = 1.3
cell_fun = function(j, i, x, y, w, h, fill){  
  grid.rect(x = x, y = y, width = w, height = h,             
            gp = gpar(col = NA, fill = NA))  
  grid.circle(x=x,y=y,r= expand_index*percent_mat[i, j]/100 * min(unit.c(w, h)),              
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

# celltype
colnames(exp_mat)
NeuronLineage_col = c("#009404FF", "#FFD600FF", "#00EDEDFF", "#0000FFFF", "#FF0000FF", "#00597f", "#808080");
names(NeuronLineage_col) = c("SubpalliumNR2F1+NR2F2+", "SubpalliumLHX8+ISL1+", "SubpalliumCRABP1+", "CortexCRABP1+", "CortexLHX6+EPHA5+", "VZ_RGC(toSVZ)", "VZ_RGC(Proliferative)")
cluster_anno <- colnames(exp_mat)
column_ha<- HeatmapAnnotation(
  cluster_anno = cluster_anno,
  col = list(cluster_anno = NeuronLineage_col),
  na_col = "grey"
)
layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= expand_index*pindex(percent_mat, i, j)/100 * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "Percent Expressed",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = expand_index*0 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = expand_index*0.25 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = expand_index*0.5 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = expand_index*0.75 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = expand_index*1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))
fig = Heatmap(exp_mat,
              cluster_rows = F, cluster_columns = F, show_heatmap_legend = F,
              heatmap_legend_param=list(title="Average Scaled Expression Level"),
              column_title = "Specific Dynamic genes in each lineage",
              col = col_fun,
              rect_gp = gpar(type = "none"),        
              cell_fun = cell_fun,        
              row_names_gp = gpar(fontsize = 7),
              # row_km = 4,
              top_annotation = column_ha,
              border = "black");fig
fig = draw(fig, annotation_legend_list = lgd_list);fig
Normal_ggplot_light(fig, f_path = "../../results/figures", 
                    normtag = "Figure5_G", nwid = 6, nhei = 9, dpi_in = 320)

## Figure 6

Age_loop = "";load("")
setname = paste0(Age_loop, "_MGE_Bin20");setname
pixel_x = max(Seurat_stereo$x) - min(Seurat_stereo$x);pixel_x
pixel_y = max(Seurat_stereo$y) - min(Seurat_stereo$y);pixel_y

Type_selected = "EPHA5MEF2C";col_used = "#FF0000FF";gene_remain = EPHA5_gene
fig <- Plot_Density_Joint_Only(seurat_object = Seurat_stereo, reduction = "RealLocation_reduction",
                               # features = intersect(gene_order_2[[2]], rownames(Seurat_stereo)),
                               features = gene_remain,
                               custom_palette = colorRampPalette(c('#176BA0','black','#EE9A3A'), 
                                                                 bias = 5
                               )(40), pt.size = 1) + theme_void() +
  # scale_fill_gradientn(colours = c('#176BA0','black','#EE9A3A'), values = c()) +
  ggtitle(Age_loop, paste0(Type_selected, " Module"));fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure_6_", Age_loop, "_", Type_selected, "_SpaceSmoothModule"),
                    nwid = round(pixel_x/600)*1, nhei = round(pixel_y/500)*1, if_egg = T)

Type_selected_name = Type_selected
fig <- ggplot(subset(Seurat_stereo_SVZNeuronalRGC, cells = WhichCells(Seurat_stereo_SVZNeuronalRGC, idents = Type_selected_name))@meta.data, aes(x, y)) + 
  # geom_point(alpha = 0.2, fill = "grey80", stroke = 0) + 
  # geom_bin2d(bins = 10, binwidth = c(20, 20)) +
  # geom_density_2d() +
  geom_density_2d_filled(data = subset(Seurat_stereo_SVZNeuronalRGC, cells = WhichCells(Seurat_stereo_SVZNeuronalRGC, idents = Type_selected_name))@meta.data, aes(x, y), contour_var = c("density", "ndensity", "count")[2], lineend = c("round", "butt", "square")[3] ) +
  theme_void() + #NoLegend() +
  ggtitle(paste0(Age_loop, " SVZ_RGC ",Type_selected, " in MGE")) +
  scale_fill_manual(values = colorRampPalette(colors = c("#FFFFFF00", col_used), bias = 0.1, alpha = T)(10));fig
Single_ggplot_light(fig, "../../results/figures",
                    ggfigtag = paste0("Figure6_", Age_loop, "_", Type_selected, "_SpaceDistribution"),
                    nwid = round(pixel_x/500)*1, nhei = round(pixel_y/500)*1, if_egg = T)


