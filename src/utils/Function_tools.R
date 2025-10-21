# -*- coding: utf-8 -*-
###############################################################################
# File: 
# Project: Human GE
# Author: Yiming Yan
# Email: yym21@mails.tsinghua.edu.cn
# Institution: Tsinghua Univ IDG/McGovern Institute for Brain Research
# Created: 2024-05-01
# Last modified: 2025-03-19
# License: MIT
# Version: 0.1.1
#-------------------------------------------------------------------------------

###############################################################################

# set_list = list(
#   "GW10_MGE" = rownames(GW10_MGE), 
#   "GW13_MGE" = rownames(GW13_MGE), 
#   "GW13_MGE_rep2" = rownames(GW13_MGE_rep2), 
#   "GW18_MGE" = rownames(GW18_MGE), 
#   "GW19_MGE" = rownames(GW19_MGE), 
#   "GW26_MGE" = rownames(GW26_MGE), 
#   "GW39_MGE" = rownames(GW39_MGE)
# )

gene_expression_stats <- function(seurat_obj, genes, celltype_col, 
                                  expr_slot = "data", min_threshold = 0) {
  # 验证输入
  if (!inherits(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("celltype_col not found in metadata")
  }
  
  # 提取RNA assay
  if (!"RNA" %in% names(seurat_obj@assays)) stop("RNA assay not found")
  rna_assay <- seurat_obj@assays$RNA
  
  # 获取表达矩阵（指定slot）
  if (!expr_slot %in% slotNames(rna_assay)) {
    stop("Specified slot not found in RNA assay")
  }
  expr_mat <- slot(rna_assay, expr_slot)
  
  # 获取counts矩阵（用于计算表达比例）
  counts_mat <- rna_assay@counts
  
  # 检查基因是否存在
  valid_genes <- genes[genes %in% rownames(expr_mat)]
  if (length(valid_genes) == 0) stop("No valid genes found")
  if (length(valid_genes) < length(genes)) {
    warning("Some genes not found: ", paste(setdiff(genes, valid_genes), collapse = ", "))
  }
  
  # 提取细胞类型信息
  cell_types <- seurat_obj@meta.data[[celltype_col]]
  
  # 准备结果数据框
  results <- data.frame(
    gene = character(),
    celltype = character(),
    avg_expr = numeric(),
    pct_expr = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 计算每个细胞类型的统计量
  unique_types <- unique(cell_types)
  
  for (ct in unique_types) {
    # 获取当前细胞类型的细胞索引
    ct_cells <- which(cell_types == ct)
    
    # 提取当前细胞类型的表达矩阵
    ct_expr <- expr_mat[valid_genes, ct_cells, drop = FALSE]
    ct_counts <- counts_mat[valid_genes, ct_cells, drop = FALSE]
    
    # 计算平均表达量
    avg_expr <- Matrix::rowMeans(ct_expr)
    
    # 计算表达比例（使用counts矩阵）
    pct_expr <- Matrix::rowMeans(ct_counts > min_threshold) * 100
    
    # 添加到结果
    results <- rbind(results, data.frame(
      gene = valid_genes,
      celltype = ct,
      avg_expr = avg_expr,
      pct_expr = pct_expr,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}


Mine_Metaneibor <- function(Sobj1, Sobj2, name_1 = "Set1", name_2 = "Set2", cal_type = 1, gene_input = NULL){
  Sobj1[['comparison_cluster_label']] <- Idents(Sobj1)
  # Sobj1[['comparison_cluster_label']] <- Sobj1$SCT_snn_res.0.95
  Sobj2[['comparison_cluster_label']] <- Idents(Sobj2)
  
  Loop_obj = Sobj1
  Sobj1 = CreateSeuratObject(counts = Loop_obj@assays$RNA, meta.data = Loop_obj@meta.data)
  Loop_obj = Sobj2
  Sobj2 = CreateSeuratObject(counts = Loop_obj@assays$RNA, meta.data = Loop_obj@meta.data)
  
  plan("multisession", workers = 1)
  Sobj1 <- NormalizeData(Sobj1, normalization.method = "LogNormalize", scale.factor = 10000) 
  Sobj2 <- NormalizeData(Sobj2, normalization.method = "LogNormalize", scale.factor = 10000) 
  plan("multisession", workers = 30)
  
  # table(Sobj2$cluster_label)
  
  Sobj1$experiment <- name_1
  Sobj2$experiment <- name_2
  
  Idents(Sobj1) = Sobj1$comparison_cluster_label
  Idents(Sobj2) = Sobj2$comparison_cluster_label
  
  merged <- merge(Sobj1, Sobj2, add.cell.ids = c("Set1", "Set2")
                  #, merge.data = FALSE, merge.dr = TRUE
  )
  DefaultAssay(merged) <- 'RNA'
  mergedSCE <- Seurat::as.SingleCellExperiment(merged, assay = 'RNA')
  
  if(cal_type == 4){
    varGenes = gene_input
  }else{
    # Sobj1_DEG <- FindAllMarkers(Sobj1,logfc.threshold = log(1.15),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
    # Sobj2_DEG <- FindAllMarkers(Sobj2,logfc.threshold = log(1.5),only.pos = TRUE, base=exp(1), return.thresh = 0.01)
    
    # Sobj1_DEG <- FindAllMarkers(Sobj1,logfc.threshold = log(1.33),only.pos = TRUE, base=exp(1), return.thresh = 0.01)
    # Sobj2_DEG <- FindAllMarkers(Sobj2,logfc.threshold = log(1.33),only.pos = TRUE, base=exp(1), return.thresh = 0.01)
    
    Sobj1_DEG <- FindAllMarkers(Sobj1,logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
    Sobj2_DEG <- FindAllMarkers(Sobj2,logfc.threshold = log(1.25),only.pos = TRUE, base=exp(1), return.thresh = 0.05)
    
    # varGenes<- union(Sobj1_DEG$gene, Sobj2_DEG$gene)
    if(cal_type == 1){
      varGenes<- intersect(Sobj1_DEG$gene, Sobj2_DEG$gene)
    }else if(cal_type == 2){
      varGenes<- union(Sobj1_DEG$gene, Sobj2_DEG$gene)
    }else{
      varGenes<- intersect(Sobj1_DEG$gene, rownames(Sobj2))
    }
  }
  
  # varGenes<- union(Sobj1_DEG$gene, Sobj2_DEG$gene)
  # varGenes = Sobj2_DEG$gene
  
  celltype_NV = MetaNeighborUS(var_genes = varGenes,
                               dat = mergedSCE,
                               study_id = mergedSCE$experiment,
                               cell_type = mergedSCE$comparison_cluster_label,
                               fast_version = TRUE)
  
  HeatCol <- colorRampPalette(colors = c('#176BA0','black','#EE9A3A'))(1000)
  fig1 = pheatmap::pheatmap(celltype_NV,
                            cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA,
                            color = HeatCol)
  
  threshMatrix<-as.matrix(celltype_NV>0.8)
  threshMatrix<-1*threshMatrix
  fig2 = pheatmap::pheatmap(threshMatrix[fig1$tree_row$order, fig1$tree_row$order], cluster_rows = F, cluster_cols = F,
                            cellwidth = 20, cellheight = 20, angle_col = "45", border_color = NA,
                            color = HeatCol)
  # return(list("Heat"=fig1, "Heat_cutoff"=fig2))
  return(celltype_NV)
}



Stereo_Dens_plot = function(Seurat_stereo, Type_selected, col_used, x_range, y_range, ref_add = FALSE){
  # Type_selected = "RGC"
  fig <- ggplot(subset(Seurat_stereo, cells = WhichCells(Seurat_stereo, idents = Type_selected))@meta.data, aes(x, y)) + 
    # geom_point(alpha = 0.2, fill = "grey80", stroke = 0) + 
    # geom_bin2d(bins = 10, binwidth = c(20, 20)) +
    # geom_density_2d() +
    xlim(x_range) + ylim(y_range) +
    geom_density_2d_filled(data = subset(Seurat_stereo, cells = WhichCells(Seurat_stereo, idents = Type_selected))@meta.data, aes(x, y), contour_var = c("density", "ndensity", "count")[2], lineend = c("round", "butt", "square")[3] ) +
    theme_void() + #NoLegend() +
    # ggtitle(paste0(Age_loop, " ",Type_selected, " in MGE")) +
    scale_fill_manual(values = colorRampPalette(colors = c("#FFFFFF00", col_used), bias = 0.1, alpha = T)(10));fig
  
  fig2 = DimPlot(subset(Seurat_stereo, cells = WhichCells(Seurat_stereo, idents = Type_selected)), reduction = "RealLocation_reduction", shuffle = T, label = F, 
                 cols = col_used) + 
    xlim(x_range) + ylim(y_range) +
    theme_bw();fig2
  
  
  fig3 = DimPlot(Seurat_stereo, group.by = "orig.ident", reduction = "RealLocation_reduction", shuffle = T, label = F,  raster=FALSE,
                 cols = "grey85") + 
    xlim(x_range) + ylim(y_range) +
    theme_bw();fig3
  
  if(ref_add){
    fig = wrap_plots(list(fig2, fig, fig3), ncol = 3)
  }else{
    fig = wrap_plots(list(fig2, fig), ncol = 2)
  }
  return(fig)
}

MultiSetVisual = function(set_list){
  library(UpSetR)
  library(openxlsx)
  
  data_UR_form <- fromList(set_list)
  fig = upset(data_UR_form,             # 数据
              nsets = length(set_list),          # 六个集
              number.angles = 10, # 柱状图上数字旋转的角度
              order.by = "freq",  # 按频次排序
              mainbar.y.label = "Intersection sizes", # 修改上方柱状图y轴名
              sets.x.label = "Numbers of symptoms",   # 左下sets的x轴名
              point.size = 3,     # 点的大小
              line.size = 0.5,    # 线的粗细
              shade.color = "grey",          # 点的背景阴影颜色
              matrix.color = "grey60",       # 点的颜色
              mb.ratio = c(0.55, 0.45),      # 上、下两部分图的占比
              main.bar.color = "cadetblue4", # 上方图柱的颜色
              sets.bar.color = "grey60",     # 左下图柱的颜色
              text.scale = 1.1) 
  
  
  # 生成所有可能的组合名称（例如 "SetA&SetB"）
  all_combinations <- unlist(
    lapply(1:length(set_list), function(n) {
      combn(names(set_list), n, FUN = paste, collapse = "&", simplify = FALSE)
    }
    )
  )
  # 计算每个组合的交集基因和大小
  intersection_list <- lapply(all_combinations, function(comb) {
    sets <- unlist(strsplit(comb, "&"))  # 拆分组合名为集合名称
    genes <- Reduce(intersect, set_list[sets])  # 计算交集基因
    if (length(genes) > 0) {  # 过滤空交集
      list(
        combination = comb,
        genes = paste(genes, collapse = ", "),  # 将基因转为逗号分隔的字符串
        size = length(genes)
      )
    } else {
      NULL  # 跳过空交集
    }
  })
  # 移除空值（由空交集产生）
  intersection_list <- Filter(Negate(is.null), intersection_list)
  # 转换为数据框（每行对应一个交集组合）
  intersection_df <- do.call(rbind, lapply(intersection_list, data.frame))
  # 按 size 降序排序
  intersection_df <- intersection_df[order(-intersection_df$size), ]
  # 重置行名（可选）
  rownames(intersection_df) <- NULL
  
  # # 生成所有可能的非空交集组合（包括单个集合）
  # all_combinations <- unlist(
  #   lapply(1:length(set_list), function(n) {
  #     combn(names(set_list), n, FUN = paste, collapse = "&", simplify = FALSE)
  #   })
  # )
  # # 计算每个组合的交集基因和大小
  # intersection_list <- lapply(all_combinations, function(comb) {
  #   sets <- unlist(strsplit(comb, "&"))
  #   list(
  #     combination = comb,
  #     genes = Reduce(intersect, set_list[sets]),
  #     size = length(Reduce(intersect, set_list[sets]))
  #   )
  # })
  # names(intersection_list) = all_combinations
  
  single_sets = names(set_list)
  
  # 转换为数据框并按大小排序
  intersection_df <- do.call(rbind, lapply(intersection_list, data.frame))
  intersection_df <- intersection_df[order(-intersection_df$size), ]
  intersection_df = intersection_df[-which(intersection_df$combination %in% single_sets), ]
  
  # 设定与UpSet图相同的nintersects参数（例如20）
  top_n <- 20
  displayed_combinations <- head(intersection_df, top_n)$combination
  # 提取这些组合对应的基因
  displayed_genes <- lapply(displayed_combinations, function(comb) {
    sets <- unlist(strsplit(comb, "&"))
    Reduce(intersect, set_list[sets])
  })
  names(displayed_genes) <- displayed_combinations
  
  # # 从元数据中提取显示的交叉组合名称
  # shown_intersects <- names(fig$New_data)[colSums(fig$New_data) > 0]
  # # 获取对应基因列表
  # final_genes <- lapply(shown_intersects, function(intersect_name) {
  #   sets <- unlist(strsplit(intersect_name, "&"))
  #   Reduce(intersect, set_list[sets])
  # })
  # names(final_genes) <- shown_intersects
  # 
  # # 创建带颜色标记的Excel表格
  # output_df <- data.frame(
  #   Intersection = names(final_genes),
  #   Gene_Count = sapply(final_genes, length),
  #   Genes = sapply(final_genes, paste, collapse = ", ")
  # )
  # 
  # # 添加样式：高亮高频交集
  # hs <- createStyle(textDecoration = "BOLD", fgFill = "#FFFF00")
  # wb <- write.xlsx(output_df, "Displayed_Intersections.xlsx",
  #                  headerStyle = hs, borders = "columns")
  # 
  # # 保存文件
  # saveWorkbook(wb, "Displayed_Intersections.xlsx", overwrite = TRUE)
  
  return(fig)
}



Seurat_reset = function(SeuratObj_in, new_idents = "reset", meta_remain = F){
  if(meta_remain){
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents, meta.data = SeuratObj_in@meta.data)
  }else{
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents)
  }
  return(SeuratObj_out)
}

Add_DimReduc_Seurat = function(Seurat_test, emb_new, Reduc_name = "new", assay_loc = "RNA"){
  # emb_new = Seurat_test@meta.data[,c("x", "y")]
  colnames(emb_new) = paste0(Reduc_name, "_", c(1, 2))
  Seurat_test[[paste0(Reduc_name, "_reduction")]] <- CreateDimReducObject(embeddings = as.matrix(emb_new), assay = assay_loc, key = paste0(Reduc_name, "_"))
  return(Seurat_test)
}

Idents_merged = function(idents.before, idents.after){
  idents.out = as.character(idents.before)
  names(idents.out) = names(idents.before)
  idents.out[names(idents.after)] = as.character(idents.after)
  return(idents.out)
}


Seurat_Bar = function(Seurat_test_sub, group_idents, color_used = NULL, plot_all = T, bar_width = 0.8){
  Seurat_test_sub$Sample = group_idents
  Seurat_test_sub$T_Type = Idents(Seurat_test_sub)
  plotC <- table(Seurat_test_sub@meta.data$Sample, Seurat_test_sub@meta.data$T_Type) %>% melt()
  colnames(plotC) <- c("Sample", "CellType", "Number")
  if(is.null(color_used)){
    pC1 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="stack")+
      # scale_fill_manual(values=celltype_colors) +
      scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="", y="Cell number")+
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC1
    pC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="fill")+
      # scale_fill_manual(values=celltype_colors) +
      scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="",y="Cell proportion")+
      scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC2
  }else{
    pC1 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="stack")+
      scale_fill_manual(values=color_used) +
      # scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="", y="Cell number")+
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC1
    pC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="fill")+
      scale_fill_manual(values=color_used) +
      # scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="",y="Cell proportion")+
      scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC2
  }
  if(plot_all){
    pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect');pC
  }else{
    pC = pC2
  }
  return(pC)
}


Seurat_Bar_dodge = function(Seurat_test_sub, group_idents, color_used = NULL, plot_all = T, bar_width = 0.8){
  Seurat_test_sub$Sample = group_idents
  Seurat_test_sub$T_Type = Idents(Seurat_test_sub)
  plotC <- table(Seurat_test_sub@meta.data$Sample, Seurat_test_sub@meta.data$T_Type) %>% melt()
  colnames(plotC) <- c("Sample", "CellType", "Number")
  # plotC
  
  # 按Sample分组，计算组内CellType的比例
  plotC <- plotC %>%
    group_by(Sample) %>%               # 按Sample分组
    mutate(                           # 添加新列
      Percent = Number / sum(Number)  # 计算比例
    ) %>%
    ungroup()                         # 解除分组
  
  if(is.null(color_used)){
    pC1 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="stack")+
      # scale_fill_manual(values=celltype_colors) +
      scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="", y="Cell number")+
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC1
    pC2 <- ggplot(data = plotC, aes(x = Sample, y = Percent, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="dodge")+
      # scale_fill_manual(values=celltype_colors) +
      scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="",y="Cell proportion")+
      scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC2
  }else{
    pC1 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="stack")+
      scale_fill_manual(values=color_used) +
      # scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="", y="Cell number")+
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC1
    pC2 <- ggplot(data = plotC, aes(x = Sample, y = Percent, fill = CellType)) +
      geom_bar(stat = "identity", width=bar_width,aes(group=CellType),position="dodge")+
      scale_fill_manual(values=color_used) +
      facet_grid(~CellType) +
      # scale_fill_npg() +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="",y="Cell proportion")+
      scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
      theme(axis.text = element_text(size=12, colour = "black"))+
      theme(axis.title.y = element_text(size=12, colour = "black"))+
      theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
      theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6));pC2
  }
  if(plot_all){
    pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect');pC
  }else{
    pC = pC2
  }
  return(pC)
}


elegant_violin=function(object, features, ident){
  VlnPlot(
    object = object,
    features = features,
    # idents = ident,
    pt.size = 0,  # 隐藏原始点
    cols = palette,
    same.y.lims = TRUE,
    stack = TRUE,
    flip = TRUE,
    fill.by = "ident",
    add.noise = FALSE,
    ncol = 1
  ) + 
    geom_boxplot(
      width = 0.02,
      fill = "white",
      color = "#404040",
      outlier.shape = NA,
      lwd = 0.5,
      position = position_dodge(width = 0.9))
}


prepare_dual_prop_data <- function(object, features, idents) {
  FetchData(object, vars = c(features, "ident")) %>%
    filter(ident %in% idents) %>%
    tidyr::pivot_longer(cols = -ident, names_to = "gene") %>%
    group_by(gene, ident) %>%
    summarise(
      detected = mean(value > 0),
      undetected = 1 - detected,
      .groups = 'drop'
    ) %>%
    tidyr::pivot_longer(
      cols = c(detected, undetected),
      names_to = "status",
      values_to = "proportion"
    ) %>%
    mutate(
      gene = factor(gene, levels = features),
      ident = factor(ident, levels = idents),
      status = factor(status, levels = c("undetected", "detected")))
}

aes_config <- list(
  detected_color = "#2c7bb6",  # 科学蓝
  undetected_color = "white",  # 纯白
  palette_background = "#f8f9fa",  # 浅灰背景
  font_family = "Arial",
  title_size = 14,
  axis_title_size = 12,
  axis_text_size = 10,
  facet_text_size = 11,
  bar_edge_color = "grey30",
  bar_width = 0.7,
  label_size = 3.8,
  label_threshold = 0.1
)


create_dual_stacked_plot <- function(data) {
  ggplot(data, aes(x = ident, y = proportion, fill = status)) +
    geom_col(
      width = aes_config$bar_width,
      color = aes_config$bar_edge_color,
      linewidth = 0.3,
      position = position_fill(reverse = TRUE)
    ) +
    geom_text(
      aes(label = ifelse(status == "detected" & proportion >= aes_config$label_threshold,
                         percent(proportion, accuracy = 1), "")),
      position = position_fill(vjust = 0.5, reverse = TRUE),
      color = "white",
      size = aes_config$label_size,
      family = aes_config$font_family,
      fontface = "bold"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      breaks = seq(0, 1, 0.25)
    ) +
    scale_fill_manual(
      values = c(
        "detected" = aes_config$detected_color,
        "undetected" = aes_config$undetected_color
      ),
      guide = "none"
    ) +
    facet_wrap(
      ~gene,
      nrow = 1,
      labeller = labeller(gene = label_wrap_gen(10)) +  # 自动换行基因名
        labs(
          x = "Cell Type",
          y = "Expression Proportion",
          title = "Gene Detection Profile Matrix",
          subtitle = "White areas represent non-detected proportions"
        ) +
        theme_minimal(base_family = aes_config$font_family) +
        theme(
          plot.title = element_text(size = aes_config$title_size, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = aes_config$title_size-2, hjust = 0.5, color = "grey40"),
          axis.title = element_text(size = aes_config$axis_title_size),
          axis.text.x = element_text(
            size = aes_config$axis_text_size,
            angle = 45,
            hjust = 1,
            color = "grey30"
          ),
          axis.text.y = element_text(size = aes_config$axis_text_size, color = "grey30"),
          panel.spacing = unit(1.2, "lines"),  # 分面间距
          strip.background = element_rect(fill = "#e9ecef", color = NA),  # 分面标题背景
          strip.text = element_text(
            size = aes_config$facet_text_size,
            face = "italic",
            margin = margin(b = 5)
          ),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
          plot.background = element_rect(fill = aes_config$palette_background, color = NA),
          panel.background = element_rect(fill = aes_config$palette_background, color = NA)
        ))
}


RCTD_test = function(ref, Seurat_test){
  ref$ref_anno = Idents(ref)
  gc()
  library(spacexr)
  ref <- UpdateSeuratObject(ref)
  Idents(ref) <- "ref_anno"
  table(Idents(ref))
  # table((ref$CellSubtype))
  
  # extract information to pass to the RCTD Reference function
  counts <- ref[["RNA"]]@counts
  cluster <- as.factor(ref$ref_anno)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)
  # set up query with the RCTD function SpatialRNA
  counts <- Seurat_test[["Spatial"]]@counts
  coords <- GetTissueCoordinates(Seurat_test)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  #
  RCTD <- create.RCTD(query, reference, max_cores = 32, CELL_MIN_INSTANCE = 15, UMI_min_sigma = c(300,5)[2])
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  Seurat_test <- AddMetaData(Seurat_test, metadata = RCTD@results$results_df)
  table(Seurat_test$first_type)
  
  return(Seurat_test)
}



DimHeatmap_Pro = function(Seurat_test, idents_selected, reduct_name = "integrated.cca", reduct_id = 1){
  pca_data <- Seurat::GetAssayData(Seurat_test, slot = "scale.data")
  top_genes <- union(rownames(top_n(x = as.data.frame(Seurat_test@reductions$integrated.cca@feature.loadings[, reduct_id]), n = 30)),
                     rownames(top_n(x = as.data.frame(-Seurat_test@reductions$integrated.cca@feature.loadings[, reduct_id]), n = 30)))
  pca_matrix <- pca_data[top_genes, ]
  annotation_col <- data.frame(
    Cluster = Seurat_test@meta.data[, idents_selected]
  )
  rownames(annotation_col) <- colnames(Seurat_test)
  fig = pheatmap(
    mat = pca_matrix,
    scale = "none",
    cluster_rows = F, cluster_cols = T,
    annotation_col = annotation_col,
    breaks = seq(-2, 2, length.out = 1000),
    show_colnames = FALSE,  # 隐藏列名（细胞名）
    color = HeatCol, main = paste0(reduct_name, "_", reduct_id)
  )
  return(fig)
}


Single_ggplot_light = function(ggfig, f_path, ggfigtag = "NoTag", nwid = 22, nhei = 21, dpi_in = 320, if_egg = FALSE){
  if(if_egg){
    nwid_out = nwid * 1.2;nhei_out = nhei * 1.2;
    
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".pdf"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
            width = nwid_out,  height = nhei_out, limitsize = FALSE,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".png"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
            width = nwid_out,  height = nhei_out, limitsize = FALSE,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".svg"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
            width = nwid_out,  height = nhei_out,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".bmp"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
            width = nwid_out,  height = nhei_out,
            units = "cm", dpi = dpi_in)
    # ggsave2(filename =  paste0(f_path,"/",ggfigtag,".ps"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
    #         width = nwid,  height = nhei,   
    #         units = "cm", dpi = dpi_in)
    # ggsave2(filename =  paste0(f_path,"/",ggfigtag,".eps"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
    #         width = nwid_out,  height = nhei_out,
    #         units = "cm", dpi = dpi_in)
    
  }else{
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".pdf"), plot = ggfig, device = cairo_pdf,
            width = nwid,  height = nhei, limitsize = FALSE,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".png"), plot = ggfig,
            width = nwid,  height = nhei, limitsize = FALSE,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".svg"), plot = ggfig,
            width = nwid,  height = nhei,
            units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".bmp"), plot = ggfig,
            width = nwid,  height = nhei,
            units = "cm", dpi = dpi_in)
    # ggsave2(filename =  paste0(f_path,"/",ggfigtag,".ps"), plot = ggfig,
    #         width = nwid,  height = nhei,   
    #         units = "cm", dpi = dpi_in)
    ggsave2(filename =  paste0(f_path,"/",ggfigtag,".eps"), plot = ggfig,
            width = nwid,  height = nhei,
            units = "cm", dpi = dpi_in)
  }
  return(paste0("Plot done in ", f_path))
}

Normal_ggplot_light = function(normal_fig, f_path, normtag = "NoTag", nwid = 22, nhei = 21, dpi_in = 320){
  # if(if_egg){
  #   nwid_out = nwid * 1.2;nhei_out = nhei * 1.2;
  #   
  #   ggsave2(filename =  paste0(f_path,"/",ggfigtag,".pdf"), plot = egg::set_panel_size(ggfig, width = unit(nwid, "cm"), height = unit(nhei, "cm")),
  #           width = nwid_out,  height = nhei_out, limitsize = FALSE,
  #           units = "cm", dpi = dpi_in)
  # }else{
  # ggsave2(filename =  paste0(f_path,"/",ggfigtag,".pdf"), plot = ggfig, device = cairo_pdf,
  #         width = nwid,  height = nhei, limitsize = FALSE,
  #         units = "cm", dpi = dpi_in)
  # }
  
  if(TRUE){
    pdf(file = paste0(f_path, "/", normtag, ".pdf"),
        width = nwid, height = nhei
    )
    print(normal_fig)
    dev.off()
    
    jpeg(filename = paste0(f_path, "/", normtag, ".jpeg"),
         width = nwid, height = nhei, res = dpi_in,
         units = "in", quality = 85)
    print(normal_fig)
    dev.off()
    
    png(filename = paste0(f_path, "/", normtag, ".png"),
        width = nwid, height = nhei, res = dpi_in,
        units = "in")
    print(normal_fig)
    dev.off()
    
    bmp(filename = paste0(f_path, "/", normtag, ".bmp"),
        width = nwid, height = nhei, res = dpi_in,
        units = "in")
    print(normal_fig)
    dev.off()
    
    # tiff(filename = "Rplot%03d.tiff",
    #      width = 480, height = 480, units = "px", pointsize = 12,
    #      compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
    #      bg = "white", res = NA,  ...,
    #      type = c("cairo", "Xlib", "quartz"), antialias)
    
    svg(filename = paste0(f_path, "/", normtag, ".svg"),
        width = nwid, height = nhei
    )
    print(normal_fig)
    dev.off()
    
    cairo_ps(filename = paste0(f_path, "/", normtag, ".svg"),
             width = nwid, height = nhei,
             fallback_resolution = dpi_in)
    print(normal_fig)
    dev.off()
  }
  
  return(paste0("Plot done in ", f_path))
}


ent_ID2Symbol = function(vector_in){
  # vector_in = data.frame()
  vector_out = lapply(vector_in, function(y){
    entID = unlist(strsplit(y, split = "/"))
    y=bitr(entID, fromType="ENTREZID",
           toType="SYMBOL",
           OrgDb = "org.Hs.eg.db")
    y=y$SYMBOL
    y=paste(y, collapse = "/")
  })
  return(as.vector(vector_out %>% unlist()))
}

com_gsea_kegg_human <- function(table_list, pro){
  library(org.Hs.eg.db) # human的OrgDB
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  
  for (i in 1:length(table_list)) {
    sub_name = names(table_list)[i]
    data_loop = table_list[[i]]
    
    # gene <- bitr(data_loop$gene, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
    # gene <- dplyr::distinct(gene, ENSEMBL, .keep_all=TRUE)
    # data_used <- data_loop %>%
    #   inner_join(gene, by = c("gene" = "SYMBOL"))
    
    data_used = data_loop
    
    data_sort <- data_used %>% 
      dplyr::distinct(gene, .keep_all=T) %>%
      arrange(desc(avg_log2FC))
    
    gene_list <- data_sort$avg_log2FC
    # names(gene_list) <- data_sort$ENSEMBL
    names(gene_list) <- data_sort$gene
    
    res <- gseGO(
      gene_list,    # 根据logFC排序的基因集
      ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
      OrgDb = org.Hs.eg.db,    # 使用人的OrgDb
      # keyType = "ENSEMBL",    # 基因id类型
      keyType = "SYMBOL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",    # p值校正方法
    )
    pdf(paste0(pro, "_", sub_name, '_GO_GSEA_BP.pdf') ,width = 15, height = 12)
    print(dotplot(
      res, 
      showCategory=10, title = paste0(sub_name, "_GO_GSEA"),
      split=".sign") + facet_grid(.~.sign))
    dev.off() 
    
    write.xlsx(res@result, 
               file = paste0(pro, "_", sub_name, '_GO_GSEA_BP.xlsx'))
    
    gene <- bitr(data_loop$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
    gene <- dplyr::distinct(gene, ENTREZID, .keep_all=TRUE)
    data_used <- data_loop %>%
      inner_join(gene, by = c("gene" = "SYMBOL"))
    
    data_sort <- data_used %>% 
      dplyr::distinct(gene, .keep_all=T) %>%
      arrange(desc(avg_log2FC))
    
    gene_list <- data_sort$avg_log2FC
    # names(gene_list) <- data_sort$ENSEMBL
    names(gene_list) <- data_sort$ENTREZID
    
    res <- gseKEGG(
      gene_list,    # 根据logFC排序的基因集
      organism = "hsa",    # 人的拉丁名缩写
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    
    pdf(paste0(pro, "_", sub_name, '_KEGG_GSEA.pdf') ,width = 15, height = 12)
    print(dotplot(
      res,
      showCategory=10, title = paste0(sub_name, "_GO_GSEA"),
      split=".sign") + facet_grid(.~.sign))
    dev.off()
    
    # go_gmt <- read.gmt("go.gmt")
    # res <- GSEA(
    #   gene_list,
    #   TERM2GENE = go_gmt)
    
    write.xlsx(res@result, 
               file = paste0(pro, "_", sub_name, '_KEGG_GSEA.xlsx'))
  }
}

# symbols_list = gene_list
com_go_kegg_ReactomePA_human <- function(symbols_list, pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  library(topGO)
  
  # 首先全部的symbol 需要转为 entrezID
  
  gcSample = lapply(symbols_list, function(y){ 
    # y=as.character(na.omit(select(org.Hs.eg.db,
    #                               keys = y,
    #                               columns = 'ENTREZID',
    #                               keytype = 'SYMBOL')[,2])
    # )
    y=bitr(y, fromType="SYMBOL",
           toType="ENTREZID",
           OrgDb = "org.Hs.eg.db")
    y=y$ENTREZID
    y
  })
  gcSample
  
  # 第1个注释是 KEGG 
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  +
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  ggsave(paste0(pro, '_kegg.pdf'), width = 12,height = 10)
  
  # xx2 <- compareCluster(gcSample, fun="enrichKEGG",
  #                      organism="hsa", pvalueCutoff=0.2)
  # View(xx2@compareClusterResult)
  
  # terms_selected = xx@compareClusterResult$Description[c(3,11,13, 18,52,55, 74,81,92)]
  # terms_selected = c(terms_selected, "Huntington disease")
  # terms_selected = terms_selected[-5]
  # terms_selected = terms_selected[c(1:3, 6:8, 5,4,9)]
  
  # xx@compareClusterResult = xx@compareClusterResult[which(xx@compareClusterResult$Description %in% terms_selected), ]
  # xx@compareClusterResult = xx@compareClusterResult[c(3,11,13, 66,52,55, 74,81,84), ]
  # xx@compareClusterResult = xx@compareClusterResult[-8, ]
  # xx@compareClusterResult = xx@compareClusterResult[-5, ]
  
  # dotplot(xx, color = "pvalue")  +
  #   scale_size (range=c (10, 20)) +
  #   scale_color_gradientn(colors = c("purple", "red"), breaks = c(1e-4, 1e-3, 2e-3, 3e-3, 6e-3)) +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
  #     panel.grid.major=element_line(colour=NA),
  #         panel.background = element_rect(fill = "transparent",colour = NA),
  #         plot.background = element_rect(fill = "transparent",colour = NA),
  #         panel.grid.minor = element_blank())
  # ggsave(paste0("~/NAS/NAS2/Public/SeqShare/10XDatasets/sc_sn_RNAseq/human_GE/10X_single_cell/GW13&GW18&GW26/Final_Figures_BI/",
  #               'FigureS4_PannelB_kegg.pdf'), width = 7, height = 10)
  
  
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro,'_kegg.xlsx'))
  
  # xx <- compareCluster(gcSample, fun="enrichDO",
  #                      organism="hsa", pvalueCutoff=0.05)
  # dotplot(xx)  +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  # ggsave(paste0(pro,'_DO.pdf'),width = 10,height = 8)
  # xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  # write.xlsx(xx@compareClusterResult, file = paste0(pro,'_DO.xlsx'))
  
  # xx <- compareCluster(gcSample, fun="groupGO",
  #                      organism="hsa", pvalueCutoff=0.05)
  # dotplot(xx)  +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  # ggsave(paste0(pro,'_groupGO.pdf'),width = 10,height = 8)
  # xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  # write.xlsx(xx@compareClusterResult, file = paste0(pro,'_groupGO.xlsx'))
  
  # 第2个注释是 ReactomePA 
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro, '_ReactomePA.pdf'), width = 12, height = 10)
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro, '_ReactomePA.xlsx'))
  
  # 然后是GO数据库的BP,CC,MF的独立注释
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont  = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15, height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  # write.xlsx(lineage1_ego@compareClusterResult, 
  #           file=paste0(pro,'_GO_BP_cluster_simplified.xlsx'))
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_BP_cluster_simplified.xlsx'))
  
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_CC_cluster_simplified.xlsx'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_MF_cluster_simplified.xlsx'))
  
  for (i in 1:length(gcSample)) {
    sub_name = names(gcSample)[i]
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_BP_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 20))
    dev.off()
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_CC_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 20))
    dev.off()
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_MF_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 10))
    dev.off()
  }
}

# symbols_list = gene_list
# pro = "~/NAS/NAS2/Public/SeqShare/10XDatasets/sc_sn_RNAseq/human_GE/10X_single_cell/GW13&GW18&GW26/20230628/GeneToSVZ_20241119"
com_go_kegg_ReactomePA_human_single <- function(symbols_list, pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  library(topGO)
  
  # 首先全部的symbol 需要转为 entrezID
  
  gcSample = lapply(symbols_list, function(y){ 
    # y=as.character(na.omit(select(org.Hs.eg.db,
    #                               keys = y,
    #                               columns = 'ENTREZID',
    #                               keytype = 'SYMBOL')[,2])
    # )
    y=bitr(y, fromType="SYMBOL",
           toType="ENTREZID",
           OrgDb = "org.Hs.eg.db")
    y=y$ENTREZID
    y
  })
  gcSample
  
  # 第1个注释是 KEGG 
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  +
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  ggsave(paste0(pro, '_kegg.pdf'), width = 12,height = 10)
  
  # xx2 <- compareCluster(gcSample, fun="enrichKEGG",
  #                      organism="hsa", pvalueCutoff=0.2)
  # View(xx2@compareClusterResult)
  
  # terms_selected = xx@compareClusterResult$Description[c(3,11,13, 18,52,55, 74,81,92)]
  # terms_selected = c(terms_selected, "Huntington disease")
  # terms_selected = terms_selected[-5]
  # terms_selected = terms_selected[c(1:3, 6:8, 5,4,9)]
  
  # xx@compareClusterResult = xx@compareClusterResult[which(xx@compareClusterResult$Description %in% terms_selected), ]
  # xx@compareClusterResult = xx@compareClusterResult[c(3,11,13, 66,52,55, 74,81,84), ]
  # xx@compareClusterResult = xx@compareClusterResult[-8, ]
  # xx@compareClusterResult = xx@compareClusterResult[-5, ]
  
  # dotplot(xx, color = "pvalue")  +
  #   scale_size (range=c (10, 20)) +
  #   scale_color_gradientn(colors = c("purple", "red"), breaks = c(1e-4, 1e-3, 2e-3, 3e-3, 6e-3)) +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
  #     panel.grid.major=element_line(colour=NA),
  #         panel.background = element_rect(fill = "transparent",colour = NA),
  #         plot.background = element_rect(fill = "transparent",colour = NA),
  #         panel.grid.minor = element_blank())
  # ggsave(paste0("~/NAS/NAS2/Public/SeqShare/10XDatasets/sc_sn_RNAseq/human_GE/10X_single_cell/GW13&GW18&GW26/Final_Figures_BI/",
  #               'FigureS4_PannelB_kegg.pdf'), width = 7, height = 10)
  
  
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro,'_kegg.xlsx'))
  
  # xx <- compareCluster(gcSample, fun="enrichDO",
  #                      organism="hsa", pvalueCutoff=0.05)
  # dotplot(xx)  +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  # ggsave(paste0(pro,'_DO.pdf'),width = 10,height = 8)
  # xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  # write.xlsx(xx@compareClusterResult, file = paste0(pro,'_DO.xlsx'))
  
  # xx <- compareCluster(gcSample, fun="groupGO",
  #                      organism="hsa", pvalueCutoff=0.05)
  # dotplot(xx)  +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  # ggsave(paste0(pro,'_groupGO.pdf'),width = 10,height = 8)
  # xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  # write.xlsx(xx@compareClusterResult, file = paste0(pro,'_groupGO.xlsx'))
  
  # 第2个注释是 ReactomePA 
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro, '_ReactomePA.pdf'), width = 12, height = 10)
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro, '_ReactomePA.xlsx'))
  
  # 然后是GO数据库的BP,CC,MF的独立注释
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont  = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15, height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  # write.xlsx(lineage1_ego@compareClusterResult, 
  #           file=paste0(pro,'_GO_BP_cluster_simplified.xlsx'))
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_BP_cluster_simplified.xlsx'))
  
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_CC_cluster_simplified.xlsx'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 12)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  lineage1_ego@compareClusterResult$Symbol = ent_ID2Symbol(lineage1_ego@compareClusterResult$geneID)
  write.xlsx(lineage1_ego@compareClusterResult, 
             file=paste0(pro,'_GO_MF_cluster_simplified.xlsx'))
  
  for (i in 1:length(gcSample)) {
    sub_name = names(gcSample)[i]
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_BP_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 20))
    dev.off()
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_CC_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 20))
    dev.off()
    
    GO_hyper <- enrichGO(
      gene          = gcSample[[i]],
      keyType = "ENTREZID",
      OrgDb         = org.Hs.eg.db,
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    GO_hyper_simp <- simplify(
      GO_hyper, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    pdf(paste0(pro, "_", sub_name, '_GO_MF_cluster_simplified_net.pdf') ,width = 20, height = 18)
    print(goplot(GO_hyper_simp, showCategory = 10))
    dev.off()
  }
}


##############
# 加载必要包
library(Seurat)
library(randomForest)
library(FNN)

# 步骤1：定义边缘细胞识别函数（基于k近邻距离）
find_border_cells <- function(seurat_obj, 
                              reduction = "pca", 
                              dims = 1:30, 
                              k = 50, 
                              threshold = 0.95) {
  # 获取低维空间坐标
  embed <- Embeddings(seurat_obj[[reduction]])[, dims]
  
  # 计算每个细胞到k近邻的平均距离
  knn_dist <- get.knn(embed, k = k)$nn.dist
  cell_scores <- rowMeans(knn_dist)
  
  # 识别距离分布前5%的细胞为边缘细胞
  border_idx <- which(cell_scores > quantile(cell_scores, threshold))
  
  return(list(
    border_cells = colnames(seurat_obj)[border_idx],
    non_border = colnames(seurat_obj)[-border_idx]
  ))
}

# 步骤2：构建分类训练集（非边缘细胞）
train_data <- function(seurat_obj, 
                       border_info,
                       marker_genes) {
  # 提取非边缘细胞表达矩阵
  non_border_exp <- GetAssayData(seurat_obj, slot = "data")[
    marker_genes, 
    border_info$non_border
  ]
  
  # 创建训练数据框
  train_df <- data.frame(
    t(as.matrix(non_border_exp)),
    cell_type = seurat_obj$seurat_clusters[border_info$non_border]
  )
  
  return(train_df)
}

# 步骤3：训练随机森林分类器
train_classifier <- function(train_df) {
  set.seed(123)
  rf_model <- randomForest(
    x = train_df[, -ncol(train_df)],
    y = as.factor(train_df$cell_type),
    ntree = 500,
    importance = TRUE
  )
  return(rf_model)
}

# 步骤4：预测边缘细胞类型
predict_border_cells <- function(seurat_obj, 
                                 border_info, 
                                 model, 
                                 marker_genes) {
  # 提取边缘细胞表达矩阵
  border_exp <- GetAssayData(seurat_obj, slot = "data")[
    marker_genes, 
    border_info$border_cells
  ]
  mat_extract = t(as.matrix(border_exp))
  colnames(mat_extract) = str_replace(colnames(mat_extract), "-", ".")
  
  # 预测概率
  pred_probs <- predict(model, 
                        newdata = mat_extract, 
                        type = "prob")
  
  # 获取最大概率类型
  final_pred <- colnames(pred_probs)[apply(pred_probs, 1, which.max)]
  
  return(data.frame(
    cell = border_info$border_cells,
    original_cluster = seurat_obj$seurat_clusters[border_info$border_cells],
    predicted_type = final_pred,
    confidence = apply(pred_probs, 1, max)
  ))
}

# 加载包（新增entropy包用于计算混合度）
library(entropy)
# 步骤1：重新定义边缘细胞（基于k近邻类型混合度）
find_border_cells_v2 <- function(seurat_obj, 
                                 reduction = "pca",
                                 dims = 1:30, 
                                 k = 30,
                                 mix_threshold = 0.5) {
  # 获取低维嵌入和细胞类型标签
  embed <- Embeddings(seurat_obj[[reduction]])[, dims]
  cell_types <- seurat_obj$seurat_clusters
  
  # 计算每个细胞的k近邻类型混合度
  knn_idx <- get.knn(embed, k = k)$nn.index
  mix_scores <- apply(knn_idx, 1, function(idx) {
    neighbor_types <- cell_types[idx]
    type_counts <- table(neighbor_types)
    # 使用熵值衡量混合度（可选替换为Gini系数等）
    ent <- entropy(type_counts, method = "ML")
    max_frac <- max(type_counts) / k
    c(entropy = ent, max_frac = max_frac)
  })
  
  # 转换为数据框
  mix_df <- data.frame(
    cell = colnames(seurat_obj),
    entropy = mix_scores["entropy", ],
    max_frac = mix_scores["max_frac", ]
  )
  
  # 定义混合细胞：最大类型占比<阈值 且 熵值>阈值
  border_cells <- mix_df$cell[mix_df$max_frac < mix_threshold & 
                                mix_df$entropy > quantile(mix_df$entropy, 0.75)]
  
  return(list(
    border_cells = border_cells,
    non_border = setdiff(colnames(seurat_obj), border_cells),
    mix_stats = mix_df
  ))
}

# 步骤2：构建概率分类模型（使用Soft Voting）
train_prob_classifier <- function(train_df) {
  library(caret)
  # 启用概率输出
  # ctrl <- trainControl(method = "cv", 
  #                      number = 5,
  #                      classProbs = TRUE)
  # 在train_prob_classifier中增加特征选择
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = TRUE
  )
  
  # 训练随机森林概率模型
  # rf_prob <- train(
  #   x = train_df[, -ncol(train_df)],
  #   y = as.factor(train_df$cell_type),
  #   method = "rf",
  #   ntree = 500,
  #   trControl = ctrl
  # )
  # 使用ranger加速随机森林
  rf_prob <- train(
    x = train_df[, -ncol(train_df)],
    y = as.factor(train_df$cell_type),
    method = "ranger",
    importance = "impurity",
    trControl = ctrl
  )
  
  return(rf_prob)
}

# # 强制统一基因名称格式（示例：转小写）
# rownames(seurat_obj) <- tolower(rownames(seurat_obj))
# train_genes <- tolower(train_genes)


# 步骤3：预测混合细胞类型（保留多可能性）
predict_mixed_cells <- function(seurat_obj, #assay_name = "RNA",
                                border_info, 
                                model, 
                                marker_genes,
                                prob_threshold = 0.6) {
  # 提取边缘细胞表达矩阵
  border_exp <- GetAssayData(seurat_obj, slot = "data",)[
    marker_genes, 
    border_info$border_cells
  ]
  
  # 检查训练数据和预测数据基因是否完全匹配
  train_genes <- colnames(train_df)[-ncol(train_df)]  # 训练用基因
  pred_genes <- rownames(seurat_obj)                  # 对象中的基因
  
  # 查找缺失基因
  missing_genes <- setdiff(train_genes, pred_genes)
  print(paste("缺失基因数量:", length(missing_genes)))
  print(head(missing_genes, 10))
  
  mat_used = t(as.matrix(border_exp))
  colnames(mat_used) = make.names(colnames(mat_used))
  
  # 预测概率矩阵
  pred_probs <- predict(model, 
                        newdata = mat_used,
                        type = "prob")
  
  ori_pred_labels <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # 若前两类型概率接近且未达阈值，标记为混合型
      names(top2)[1]
    } else {
      names(top2)[1]
    }
  })
  
  ori_pred_labels_2 <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # 若前两类型概率接近且未达阈值，标记为混合型
      names(top2)[2]
    } else {
      names(top2)[2]
    }
  })
  
  # 判断多可能性：
  pred_labels <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # 若前两类型概率接近且未达阈值，标记为混合型
      paste(names(top2)[1:2], collapse = "/")
    } else {
      names(top2)[1]
    }
  })
  
  # 构建结果表
  results <- data.frame(
    cell = border_info$border_cells,
    original_cluster = seurat_obj$seurat_clusters[border_info$border_cells],
    predicted_type = pred_labels,
    top_prob = apply(pred_probs, 1, max),
    top_prob_type = ori_pred_labels,
    second_prob_type = ori_pred_labels_2
    # ,
    # mix_entropy = border_info$mix_stats$entropy[match(border_info$border_cells,
    #                                                   border_info$mix_stats$cell)]
  )
  return(results)
}

Seurat_Single_Batch_pipeline = function(Seurat_test, CC_out = FALSE){
  dim_cut = 30
  Seurat_test = Seurat_test %>% 
    Seurat_reset(meta_remain = T)
  plan("multisession", workers = 1)
  Seurat_test <- Seurat_test %>% NormalizeData() %>% FindVariableFeatures()
  if(CC_out){
    Seurat_test <- CellCycleScoring(Seurat_test, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
    Seurat_test$CC.Difference <- Seurat_test$S.Score - Seurat_test$G2M.Score
    print("CC out")
  }else{
    print("No CC out")
  }
  
  Seurat_test <- ScaleData(Seurat_test, vars.to.regress = "CC.Difference", features = rownames(Seurat_test)) %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:dim_cut, reduction = "pca") %>%
    RunUMAP(reduction = "pca", dims = 1:dim_cut, n.neighbors = 30L)
  
  for (minD in c(0.05, 0.1, 0.3, 0.5)) {
    gc()
    Seurat_test = Seurat_test %>%
      RunUMAP(reduction = "pca", dims = 1:dim_cut, n.neighbors = 30L, min.dist = minD, 
              reduction.name = paste0("mD", minD, "umap"))
  }
  
  for (res in seq(0.2, 1, 0.2)){
    plan("multisession", workers = 1)
    Seurat_test <- FindClusters(Seurat_test, resolution = res, cluster.name = paste0("clusters_res_", res), algorithm = 1)
    plan("multisession", workers = 16)
    gc()
  }
  return(Seurat_test)
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}
unregister_dopar()
