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

Seurat_reset = function(SeuratObj_in, new_idents = "reset", meta_remain = F){
  if(meta_remain){
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents, meta.data = SeuratObj_in@meta.data)
  }else{
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents)
  }
  return(SeuratObj_out)
}

Basic_Seurat_QC <- function(loop_seurat, if_sct = FALSE, if_mouse = FALSE, CC_out = TRUE){
  if(if_mouse){
    # grep("^MT-", rownames(loop_seurat))
    loop_seurat[["percent.mt"]] <- PercentageFeatureSet(loop_seurat, pattern = "^mt-")
    # grep("^RP[SL]", rownames(loop_seurat))
    loop_seurat[["percent.rp"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Rp[sl]")
    loop_seurat[["percent.ercc"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Ercc")
    
    # HB.genes <- c(,"Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y""HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    # HB_m <- intersect(rownames(loop_seurat), HB.genes)
    loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Hb[abdegmqz]")
    
  }else{
    # grep("^MT-", rownames(loop_seurat))
    loop_seurat[["percent.mt"]] <- PercentageFeatureSet(loop_seurat, pattern = "^MT-")
    # grep("^RP[SL]", rownames(loop_seurat))
    loop_seurat[["percent.rp"]] <- PercentageFeatureSet(loop_seurat, pattern = "^RP[SL]")
    loop_seurat[["percent.ercc"]] <- PercentageFeatureSet(loop_seurat, pattern = "^ERCC")
    
    # HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    # HB_m <- intersect(rownames(loop_seurat), HB.genes)
    # loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, features = HB_m)
    loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, pattern = "^HB[ABDEGMQZ]")
    
  }
  
  DefaultAssay(loop_seurat) = "RNA"
  plan("multisession", workers = 1)
  loop_seurat <- loop_seurat %>% NormalizeData()
  
  if(CC_out){
    loop_seurat <- CellCycleScoring(loop_seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
    loop_seurat$CC.Difference <- loop_seurat$S.Score - loop_seurat$G2M.Score
  }else{
    print("No CC out")
  }
  
  # loop_seurat = Find_doublet(loop_seurat, sct_state = if_sct)
  return(loop_seurat)
}

gene_expression_stats <- function(seurat_obj, genes, celltype_col, 
                                  expr_slot = "data", min_threshold = 0) {
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("celltype_col not found in metadata")
  }
  
  # Extract RNA assay
  if (!"RNA" %in% names(seurat_obj@assays)) stop("RNA assay not found")
  rna_assay <- seurat_obj@assays$RNA
  
  # Get expression matrix (specified slot)
  if (!expr_slot %in% slotNames(rna_assay)) {
    stop("Specified slot not found in RNA assay")
  }
  expr_mat <- slot(rna_assay, expr_slot)
  
  # Get counts matrix (for calculating expression proportion)
  counts_mat <- rna_assay@counts
  
  # Check if genes exist
  valid_genes <- genes[genes %in% rownames(expr_mat)]
  if (length(valid_genes) == 0) stop("No valid genes found")
  if (length(valid_genes) < length(genes)) {
    warning("Some genes not found: ", paste(setdiff(genes, valid_genes), collapse = ", "))
  }
  
  # Extract cell type information
  cell_types <- seurat_obj@meta.data[[celltype_col]]
  
  # Prepare results data frame
  results <- data.frame(
    gene = character(),
    celltype = character(),
    avg_expr = numeric(),
    pct_expr = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate statistics for each cell type
  unique_types <- unique(cell_types)
  
  for (ct in unique_types) {
    # Get cell indices for current cell type
    ct_cells <- which(cell_types == ct)
    
    # Extract expression matrix for current cell type
    ct_expr <- expr_mat[valid_genes, ct_cells, drop = FALSE]
    ct_counts <- counts_mat[valid_genes, ct_cells, drop = FALSE]
    
    # Calculate average expression
    avg_expr <- Matrix::rowMeans(ct_expr)
    
    # Calculate expression proportion (using counts matrix)
    pct_expr <- Matrix::rowMeans(ct_counts > min_threshold) * 100
    
    # Add to results
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
  fig = upset(data_UR_form,             # data
              nsets = length(set_list),          # number of sets
              number.angles = 10, # angle of numbers on bar chart
              order.by = "freq",  # order by frequency
              mainbar.y.label = "Intersection sizes", # modify upper bar chart y-axis label
              sets.x.label = "Set sizes",   # lower left sets x-axis label
              point.size = 3,     # point size
              line.size = 0.5,    # line thickness
              shade.color = "grey",          # point background shadow color
              matrix.color = "grey60",       # point color
              mb.ratio = c(0.55, 0.45),      # ratio of upper/lower parts
              main.bar.color = "cadetblue4", # upper bar color
              sets.bar.color = "grey60",     # lower bar color
              text.scale = 1.1) 
  
  
  # Generate all possible combination names (e.g., "SetA&SetB")
  all_combinations <- unlist(
    lapply(1:length(set_list), function(n) {
      combn(names(set_list), n, FUN = paste, collapse = "&", simplify = FALSE)
    }
    )
  )
  # Calculate intersecting genes and size for each combination
  intersection_list <- lapply(all_combinations, function(comb) {
    sets <- unlist(strsplit(comb, "&"))  # Split combination name into set names
    genes <- Reduce(intersect, set_list[sets])  # Calculate intersecting genes
    if (length(genes) > 0) {  # Filter out empty intersections
      list(
        combination = comb,
        genes = paste(genes, collapse = ", "),  # Convert genes to comma-separated string
        size = length(genes)
      )
    } else {
      NULL  # Skip empty intersections
    }
  })
  # Remove empty values (from empty intersections)
  intersection_list <- Filter(Negate(is.null), intersection_list)
  # Convert to dataframe (each row corresponds to an intersection combination)
  intersection_df <- do.call(rbind, lapply(intersection_list, data.frame))
  # Sort by size in descending order
  intersection_df <- intersection_df[order(-intersection_df$size), ]
  # Reset row names (optional)
  rownames(intersection_df) <- NULL
  
  # # Generate all possible non-empty intersection combinations (including individual sets)
  # all_combinations <- unlist(
  #   lapply(1:length(set_list), function(n) {
  #     combn(names(set_list), n, FUN = paste, collapse = "&", simplify = FALSE)
  #   })
  # )
  # # Calculate intersecting genes and size for each combination
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
  
  # Convert to dataframe and sort by size
  intersection_df <- do.call(rbind, lapply(intersection_list, data.frame))
  intersection_df <- intersection_df[order(-intersection_df$size), ]
  intersection_df = intersection_df[-which(intersection_df$combination %in% single_sets), ]
  
  # Set nintersects parameter same as UpSet plot (e.g., 20)
  top_n <- 20
  displayed_combinations <- head(intersection_df, top_n)$combination
  # Extract genes corresponding to these combinations
  displayed_genes <- lapply(displayed_combinations, function(comb) {
    sets <- unlist(strsplit(comb, "&"))
    Reduce(intersect, set_list[sets])
  })
  names(displayed_genes) <- displayed_combinations
  
  # # Extract displayed intersection combination names from metadata
  # shown_intersects <- names(fig$New_data)[colSums(fig$New_data) > 0]
  # # Get corresponding gene lists
  # final_genes <- lapply(shown_intersects, function(intersect_name) {
  #   sets <- unlist(strsplit(intersect_name, "&"))
  #   Reduce(intersect, set_list[sets])
  # })
  # names(final_genes) <- shown_intersects
  # 
  # # Create Excel table with color markers
  # output_df <- data.frame(
  #   Intersection = names(final_genes),
  #   Gene_Count = sapply(final_genes, length),
  #   Genes = sapply(final_genes, paste, collapse = ", ")
  # )
  # 
  # # Add style: highlight high-frequency intersections
  # hs <- createStyle(textDecoration = "BOLD", fgFill = "#FFFF00")
  # wb <- write.xlsx(output_df, "Displayed_Intersections.xlsx",
  #                  headerStyle = hs, borders = "columns")
  # 
  # # Save file
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
      scale_y_continuous(labels = percent)+ ####Used to adjust y-axis position
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
      scale_y_continuous(labels = percent)+ ####Used to adjust y-axis position
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
  
  # Group by Sample and calculate the proportion of each CellType within groups
  plotC <- plotC %>%
    group_by(Sample) %>%               # Group by Sample
    mutate(                           # Add new column
      Percent = Number / sum(Number)  # Calculate proportion
    ) %>%
    ungroup()                         # Ungroup
  
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
      scale_y_continuous(labels = percent)+ ####Used to adjust y-axis position
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
      scale_y_continuous(labels = percent)+ ####Used to adjust y-axis position
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
    pt.size = 0,  # Hide original points
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
  detected_color = "#2c7bb6",  # scientific blue
  undetected_color = "white",  # pure white
  palette_background = "#f8f9fa",  # light gray background
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
      labeller = labeller(gene = label_wrap_gen(10)) +  # Automatic line wrapping for gene names
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
          panel.spacing = unit(1.2, "lines"),  # Facet spacing
          strip.background = element_rect(fill = "#e9ecef", color = NA),  # Facet title background
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
    show_colnames = FALSE,  # Hide column names (cell names)
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
  library(org.Hs.eg.db) # human OrgDB
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
      gene_list,    # Gene set sorted by logFC
      ont = "BP",    # Optional: "BP", "MF", "CC" or "ALL"
      OrgDb = org.Hs.eg.db,    # Use human OrgDb
      # keyType = "ENSEMBL",    # Gene ID type
      keyType = "SYMBOL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",    # p-value adjustment method
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
      gene_list,    # Gene set sorted by logFC
      organism = "hsa",    # Human Latin name abbreviation
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
  
  # First, convert all symbols to entrezID
  
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
  
  # First annotation: KEGG 
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
  
  # Second annotation: ReactomePA 
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro, '_ReactomePA.pdf'), width = 12, height = 10)
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro, '_ReactomePA.xlsx'))
  
  # Then perform independent annotations for GO database's BP, CC, MF categories
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
  
  # First, convert all symbols to entrezID
  
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
  
  # First annotation: KEGG 
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
  
  # Second annotation: ReactomePA 
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro, '_ReactomePA.pdf'), width = 12, height = 10)
  xx@compareClusterResult$Symbol = ent_ID2Symbol(xx@compareClusterResult$geneID)
  write.xlsx(xx@compareClusterResult, file = paste0(pro, '_ReactomePA.xlsx'))
  
  # Then perform independent annotations for GO database's BP, CC, MF categories
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
# Load required packages
library(Seurat)
library(randomForest)
library(FNN)

# Step 1: Define border cell identification function (based on k-nearest neighbor distance)
find_border_cells <- function(seurat_obj, 
                              reduction = "pca", 
                              dims = 1:30, 
                              k = 50, 
                              threshold = 0.95) {
  # Get low-dimensional space coordinates
  embed <- Embeddings(seurat_obj[[reduction]])[, dims]
  
  # Calculate average distance from each cell to k-nearest neighbors
  knn_dist <- get.knn(embed, k = k)$nn.dist
  cell_scores <- rowMeans(knn_dist)
  
  # Identify cells in the top 5% of distance distribution as border cells
  border_idx <- which(cell_scores > quantile(cell_scores, threshold))
  
  return(list(
    border_cells = colnames(seurat_obj)[border_idx],
    non_border = colnames(seurat_obj)[-border_idx]
  ))
}

# Step 2: Build classification training set (non-border cells)
train_data <- function(seurat_obj, 
                       border_info,
                       marker_genes) {
  # Extract non-border cell expression matrix
  non_border_exp <- GetAssayData(seurat_obj, slot = "data")[
    marker_genes, 
    border_info$non_border
  ]
  
  # Create training dataframe
  train_df <- data.frame(
    t(as.matrix(non_border_exp)),
    cell_type = seurat_obj$seurat_clusters[border_info$non_border]
  )
  
  return(train_df)
}

# Step 3: Train random forest classifier
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

# Step 4: Predict border cell types
predict_border_cells <- function(seurat_obj, 
                                 border_info, 
                                 model, 
                                 marker_genes) {
  # Extract border cell expression matrix
  border_exp <- GetAssayData(seurat_obj, slot = "data")[
    marker_genes, 
    border_info$border_cells
  ]
  mat_extract = t(as.matrix(border_exp))
  colnames(mat_extract) = str_replace(colnames(mat_extract), "-", ".")
  
  # Predict probabilities
  pred_probs <- predict(model, 
                        newdata = mat_extract, 
                        type = "prob")
  
  # Get highest probability type
  final_pred <- colnames(pred_probs)[apply(pred_probs, 1, which.max)]
  
  return(data.frame(
    cell = border_info$border_cells,
    original_cluster = seurat_obj$seurat_clusters[border_info$border_cells],
    predicted_type = final_pred,
    confidence = apply(pred_probs, 1, max)
  ))
}

# Load packages (added entropy package for calculating mixing degree)
library(entropy)
# Step 1: Redefine border cells (based on k-nearest neighbor type mixing degree)
find_border_cells_v2 <- function(seurat_obj, 
                                 reduction = "pca",
                                 dims = 1:30, 
                                 k = 30,
                                 mix_threshold = 0.5) {
  # Get low-dimensional embeddings and cell type labels
  embed <- Embeddings(seurat_obj[[reduction]])[, dims]
  cell_types <- seurat_obj$seurat_clusters
  
  # Calculate k-nearest neighbor type mixing degree for each cell
  knn_idx <- get.knn(embed, k = k)$nn.index
  mix_scores <- apply(knn_idx, 1, function(idx) {
    neighbor_types <- cell_types[idx]
    type_counts <- table(neighbor_types)
    # Use entropy to measure mixing degree (can be replaced with Gini coefficient, etc.)
    ent <- entropy(type_counts, method = "ML")
    max_frac <- max(type_counts) / k
    c(entropy = ent, max_frac = max_frac)
  })
  
  # Convert to dataframe
  mix_df <- data.frame(
    cell = colnames(seurat_obj),
    entropy = mix_scores["entropy", ],
    max_frac = mix_scores["max_frac", ]
  )
  
  # Define mixed cells: maximum type proportion < threshold AND entropy > threshold
  border_cells <- mix_df$cell[mix_df$max_frac < mix_threshold & 
                                mix_df$entropy > quantile(mix_df$entropy, 0.75)]
  
  return(list(
    border_cells = border_cells,
    non_border = setdiff(colnames(seurat_obj), border_cells),
    mix_stats = mix_df
  ))
}

# Step 2: Build probability classification model (using Soft Voting)
train_prob_classifier <- function(train_df) {
  library(caret)
  # Enable probability output
  # ctrl <- trainControl(method = "cv", 
  #                      number = 5,
  #                      classProbs = TRUE)
  # Add feature selection in train_prob_classifier
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = TRUE
  )
  
  # Train random forest probability model
  # rf_prob <- train(
  #   x = train_df[, -ncol(train_df)],
  #   y = as.factor(train_df$cell_type),
  #   method = "rf",
  #   ntree = 500,
  #   trControl = ctrl
  # )
  # Use ranger to accelerate random forest
  rf_prob <- train(
    x = train_df[, -ncol(train_df)],
    y = as.factor(train_df$cell_type),
    method = "ranger",
    importance = "impurity",
    trControl = ctrl
  )
  
  return(rf_prob)
}

# # Force uniform gene name format (example: convert to lowercase)
# rownames(seurat_obj) <- tolower(rownames(seurat_obj))
# train_genes <- tolower(train_genes)


# Step 3: Predict mixed cell types (retain multiple possibilities)
predict_mixed_cells <- function(seurat_obj, #assay_name = "RNA",
                                border_info, 
                                model, 
                                marker_genes,
                                prob_threshold = 0.6) {
  # Extract border cell expression matrix
  border_exp <- GetAssayData(seurat_obj, slot = "data",)[
    marker_genes, 
    border_info$border_cells
  ]
  
  # Check if genes in training data and prediction data match exactly
  train_genes <- colnames(train_df)[-ncol(train_df)]  # Training genes
  pred_genes <- rownames(seurat_obj)                  # Genes in the object
  
  # Find missing genes
  missing_genes <- setdiff(train_genes, pred_genes)
  print(paste("Number of missing genes:", length(missing_genes)))
  print(head(missing_genes, 10))
  
  mat_used = t(as.matrix(border_exp))
  colnames(mat_used) = make.names(colnames(mat_used))
  
  # Prediction probability matrix
  pred_probs <- predict(model, 
                        newdata = mat_used,
                        type = "prob")
  
  ori_pred_labels <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # If top two type probabilities are close and below threshold, mark as mixed type
      names(top2)[1]
    } else {
      names(top2)[1]
    }
  })
  
  ori_pred_labels_2 <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # If top two type probabilities are close and below threshold, mark as mixed type
      names(top2)[2]
    } else {
      names(top2)[2]
    }
  })
  
  # Determine multiple possibilities:
  pred_labels <- apply(pred_probs, 1, function(probs) {
    top2 <- sort(probs, decreasing = TRUE)[1:2]
    if (top2[1] < prob_threshold & (top2[1] - top2[2]) < 0.2) {
      # If top two type probabilities are close and below threshold, mark as mixed type
      paste(names(top2)[1:2], collapse = "/")
    } else {
      names(top2)[1]
    }
  })
  
  # Build result table
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
    print("Cell cycle effect regression applied")
  }else{
    print("No cell cycle effect regression")
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

Seurat_reset = function(SeuratObj_in, new_idents = "reset", meta_remain = F){
  if(meta_remain){
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents, meta.data = SeuratObj_in@meta.data)
  }else{
    SeuratObj_out = CreateSeuratObject(SeuratObj_in@assays$RNA@counts, project = new_idents)
  }
  return(SeuratObj_out)
}

Basic_Seurat_QC <- function(loop_seurat, if_sct = FALSE, if_mouse = FALSE, CC_out = TRUE){
  if(if_mouse){
    # grep("^MT-", rownames(loop_seurat))
    loop_seurat[["percent.mt"]] <- PercentageFeatureSet(loop_seurat, pattern = "^mt-")
    # grep("^RP[SL]", rownames(loop_seurat))
    loop_seurat[["percent.rp"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Rp[sl]")
    loop_seurat[["percent.ercc"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Ercc")
    
    # HB.genes <- c(,"Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y""HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    # HB_m <- intersect(rownames(loop_seurat), HB.genes)
    loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, pattern = "^Hb[abdegmqz]")
    
  }else{
    # grep("^MT-", rownames(loop_seurat))
    loop_seurat[["percent.mt"]] <- PercentageFeatureSet(loop_seurat, pattern = "^MT-")
    # grep("^RP[SL]", rownames(loop_seurat))
    loop_seurat[["percent.rp"]] <- PercentageFeatureSet(loop_seurat, pattern = "^RP[SL]")
    loop_seurat[["percent.ercc"]] <- PercentageFeatureSet(loop_seurat, pattern = "^ERCC")
    
    # HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    # HB_m <- intersect(rownames(loop_seurat), HB.genes)
    # loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, features = HB_m)
    loop_seurat[["percent.hb"]] <- PercentageFeatureSet(loop_seurat, pattern = "^HB[ABDEGMQZ]")
    
  }
  
  DefaultAssay(loop_seurat) = "RNA"
  plan("multisession", workers = 1)
  loop_seurat <- loop_seurat %>% NormalizeData()
  
  if(CC_out){
    loop_seurat <- CellCycleScoring(loop_seurat, s.features = intersect(rownames(loop_seurat), cc.genes$s.genes), g2m.features = intersect(rownames(loop_seurat), cc.genes$g2m.genes), set.ident = FALSE)
    loop_seurat$CC.Difference <- loop_seurat$S.Score - loop_seurat$G2M.Score
  }else{
    print("No CC out")
  }
  
  # loop_seurat = Find_doublet(loop_seurat, sct_state = if_sct)
  return(loop_seurat)
}

library(DoubletFinder)

mine_paramSweep = function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  require(parallel)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

mine_doubletF = function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                          sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, 
                          nExp, sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets = neighbors[neighbors > 
                                                    n_real.cells]
          if (length(neighbors_that_are_doublets) > 
              0) {
            neighbor_types[i, ] <- table(doublet_types1[neighbors_that_are_doublets - 
                                                          n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets - 
                                                                                                  n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, 
            ]/sum(neighbor_types[i, ])
          }
          else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, 
                          nExp, sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}

DF_finder = function(Seurat_test){
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  # Seurat_test <- NormalizeData(Seurat_test)
  Seurat_test <- FindVariableFeatures(Seurat_test, selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:10)
  
  Seurat_test = Seurat_test %>% FindNeighbors(dims = 1:10) %>% 
    FindClusters(resolution = 0.3)
  ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
  # Seurat_test <- CreateSeuratObject(kidney.data)
  # Seurat_test <- SCTransform(Seurat_test)
  # Seurat_test <- RunPCA(Seurat_test)
  # Seurat_test <- RunUMAP(Seurat_test, dims = 1:10)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- mine_paramSweep(Seurat_test, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
  # sweep.res.list_kidney <- paramSweep(Seurat_test, PCs = 1:10, sct = FALSE)
  # gt.calls <- Seurat_test@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
  # sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
  # bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations = Idents(Seurat_test)
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- Seurat_test@meta.data$ClusteringResults
  nExp_poi <- round(0.025*nrow(Seurat_test@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  Seurat_test <- mine_doubletF(Seurat_test, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  Seurat_test <- mine_doubletF(Seurat_test, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = colnames(Seurat_test@meta.data)[length(colnames(Seurat_test@meta.data)) - 1], sct = FALSE)
  colnames(Seurat_test@meta.data)[(length(colnames(Seurat_test@meta.data)) - 2):length(colnames(Seurat_test@meta.data))] = c("DLF_score", "DLF_step1", "DLF_step2")
  Seurat_test = Seurat_reset(Seurat_test, meta_remain = T)
  return(Seurat_test)
}





unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}
unregister_dopar()
