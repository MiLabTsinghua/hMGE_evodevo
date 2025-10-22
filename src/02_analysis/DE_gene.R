
## Figure 1

#### Figure 1_C
markers_test = FindAllMarkers(Seurat_test, slot = "data", only.pos = T, logfc.threshold = log2(1.4), return.thresh = 1e-7)
test.markers = markers_test;test.markers_filtered = test.markers %>% group_by(cluster) %>%
  # slice_max(n = 10, order_by = avg_log2FC) %>% as.data.frame()
  slice_min(n = 10, order_by = p_val_adj, with_ties = F) %>% as.data.frame()

Anno = c(
  # Progenitor 
  # "VIM", "HES1", "FABP7", "NES","PCNA",
  # "ASCL1","PTTG1",
  # MGE-derive "PLS3", "SOX6",
  "LHX6", "NKX2-1",
  # OPC "OLIG1", "EGFR",
  "PDGFRA", "PCDH15", 
  # APC  "S100B", "SPARCL1","AQP4",  "ALDOC", 
  "GFAP", "SOX9", 
  # Endothelium "MFSD2A",  "CD34",
  "TIE1",  "PECAM1",
  # Pericyte "KCNJ8", "PDGFRB",
  "ANPEP", "RGS5", 
  # Microglia "P2RY12", "ABI3",
  "TMEM119",  "CX3CR1"
);
gene_selected = c(union(test.markers_filtered$gene[1:10],  Anno[1:3]),
                  union(test.markers_filtered$gene[11:20], Anno[4:5]),
                  union(test.markers_filtered$gene[21:29], Anno[8:9]),
                  union(test.markers_filtered$gene[31:40], Anno[6:7]),
                  union(test.markers_filtered$gene[41:50], Anno[10:11]),
                  union(test.markers_filtered$gene[51:60], Anno[12:13]),
                  union(test.markers_filtered$gene[61:70], Anno[14:15])
)

write.xlsx(markers_test, file = "../../results/figures/Figure1_C_FindMarkerResults.xlsx")
save(test.markers_filtered, Anno, gene_selected, file = "../../data/processed/Figure1_C.RData")

## Figure 2

#### Figure 2_C

Seurat_test$Revision_CellSubtype = factor(Seurat_test$Revision_CellSubtype, levels = c("VZ_RGC", "SVZ_RGC", "IPC"));Idents(Seurat_test) = "Revision_CellSubtype"
test_markers = FindAllMarkers(Seurat_test, only.pos = T)
test.markers_filtered = test_markers %>%
  group_by(cluster) %>%
  # slice_max(n = 10, order_by = avg_log2FC) %>% as.data.frame()
  slice_min(n = 10, order_by = p_val_adj, with_ties = F) %>% as.data.frame()
Anno = list(# RGC
  c("HES1", "CLU", "FABP7", "LIX1", "PTN"),
  c("VIM", "SOX5", "NES", "SOX2"), #"", #"OLIG2", # "HES5"
  c("DACH1", "CACNA1E", "FBLN7"),
  c("ASCL1", "DLX2"));

write.xlsx(markers_test, file = "../../results/figures/Figure2_C_FindMarkerResults.xlsx")
save(test.markers_filtered, Anno, file = "../../data/processed/Figure2_C.RData")

#### Figure 2_F

Idents(Seurat_stereo_sum) = factor(Seurat_stereo_sum$plot_type, levels = c("Progenitor", "MGE_derived_Neuron", 
                                                                           "OPC", "Astrocyte", 
                                                                           "Endothelium", "Pericyte", "Microglia"
))
markers_test = FindAllMarkers(Seurat_stereo_sum, slot = "data", only.pos = T, logfc.threshold = log2(1.15), return.thresh = 1e-3);
test.markers = markers_test;table(test.markers$cluster);test.markers_filtered = test.markers %>%
  group_by(cluster) %>%
  # slice_max(n = 10, order_by = avg_log2FC) %>% as.data.frame()
  slice_min(n = 6, order_by = p_val_adj, with_ties = F) %>% as.data.frame()
Anno = list(# RGC
  # c(
  #   "SOX2", "NES", "HES5"
  #   # , "HES1", "VIM"
  #   # ,"HOPX", "CLU", "RGS6"
  # ),
  # Progenitor
  c("MKI67", "TOP2A"#, "DLX1", "DLX2"
  ),
  # IPC
  # c(
  #   "DLX1", "DLX2"#, "ASCL1", "DLX6"
  #   # ,"FABP7"#, "EGFR"
  # ),
  # Progenitor "VIM", "HES1", "FABP7", "NES","PCNA","ASCL1",
  # "TOP2A", "PTTG1",  "MKI67", "DLX1",
  # MGE-derive
  c(
    "LHX6"#, "PLS3", "SOX6", "NKX2-1"
  ),
  # OPC 
  c(
    "PDGFRA"#,"OLIG1",  "PCDH15", "EGFR"
  ),
  # APC  "S100B", "SPARCL1",
  c(
    "AQP4", "SOX9"#,"GFAP",   "ALDOC"
  ), 
  # Endothelium
  c(
    "TIE1", "PECAM1"#"MFSD2A", "CD34", 
  ),
  # Pericyte
  c(
    "RGS5"#, "PDGFRB","ANPEP", "KCNJ8", 
  ),
  # Microglia
  c(
    "CX3CR1"#"ABI3", "TMEM119", "P2RY12", 
  )
);

write.xlsx(markers_test, file = "../../results/figures/Figure2_F_FindMarkerResults.xlsx")
save(test.markers_filtered, Anno, file = "../../data/processed/Figure2_F.RData")
