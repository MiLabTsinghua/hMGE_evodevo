load("../../data/external/YingchaoS_MGE_GW9.RData")
Seurat_test = Seurat_reset(Seurat_test, new_idents = "GW9_MGE", meta_remain = F)
Seurat_test$Age = "GW9";Seurat_test$batch = "GW9_XiaoqunW";Seurat_test$Region = "GE";Seurat_test$Sample = "GW9_GE"
Seurat_test = Basic_Seurat_QC(Seurat_test);Seurat_test = DF_finder(Seurat_test)

# you can get them at https://ngdc.cncb.ac.cn/omix/ 
GW10_MGE = Read10X_h5("../../data/processed/Mida_GE_GW10_20230103_GW10-M/cellbender_filtered.h5")
GW13_MGE = Read10X_h5("../../data/processed/Mida_GE_GW13_20230214_gw13-mge/cellbender_filtered.h5")
GW13_MGE_rep2 = scCustomize::Read_CellBender_h5_Mat("../../data/processed/GW13-MGE_20250504/cellbender_filtered.h5")
GW18_MGE = Read10X_h5("../../data/processed/Mida_GE_GW18_20230214_gw18-mge/cellbender_filtered.h5")
GW19_MGE = scCustomize::Read_CellBender_h5_Mat("../../data/processed/GW19-MGE_20250430/cellbender_filtered.h5")
GW26_MGE = Read10X_h5("../../data/processed/Mida_GW26_20230305_GW26-MGE/cellbender_filtered.h5")
GW39_MGE = Read10X_h5("../../data/processed/MGE_GW39_20240122/forcesample_filtered_feature_bc_matrix.h5")

Age = c("GW10", "GW13", "GW13", "GW18", "GW19", "GW26", "GW39");setname_all = c("GW10_MGE", "GW13_MGE", "GW13_MGE_rep2", "GW18_MGE", "GW19_MGE", "GW26_MGE", "GW39_MGE")
dgCM_list = list(GW10_MGE, GW13_MGE, GW13_MGE_rep2, GW18_MGE, GW19_MGE, GW26_MGE, GW39_MGE);names(dgCM_list) = setname_all
DataIntegrated.list = c()
for (i in setname_all[1:6]) {
  # loop_seurat = CreateSeuratObject(counts = dgCM_list[[i]], project = i, min.features = 10)
  loop_seurat = CreateSeuratObject(counts = dgCM_list[[i]], project = i, min.features = 200)
  loop_seurat = RenameCells(loop_seurat, i);loop_seurat$Sample = i; 
  loop_seurat$Age = substr(i, 1,4);loop_seurat$Region = "MGE";
  loop_seurat = loop_seurat %>% 
    Basic_Seurat_QC() %>% 
    DF_finder()
  loop_seurat = subset(loop_seurat, (percent.mt <= 5) & (percent.rp <= 30))
  DataIntegrated.list = c(DataIntegrated.list, list(loop_seurat))
}

i = setname_all[7];i
loop_seurat = CreateSeuratObject(counts = dgCM_list[[i]], project = i, min.features = 200)
loop_seurat = RenameCells(loop_seurat, i);loop_seurat$Sample = i; loop_seurat$Age = substr(i, 1,4);loop_seurat$Region = "MGE";
loop_seurat = loop_seurat %>% 
  Basic_Seurat_QC() %>% 
  DF_finder()
loop_seurat = subset(loop_seurat, (percent.mt <= 10) )
dim(loop_seurat)
DataIntegrated.list = c(DataIntegrated.list, list(loop_seurat))
names(DataIntegrated.list) = setname_all

Seurat_test <- merge(x = Seurat_test, y = DataIntegrated.list)
Seurat_test$orig.ident[which(Seurat_test$Age == "GW9")] = "GW9_GE"
Seurat_test$batch[which(Seurat_test$Sample != "GW9_GE")] = Seurat_test$Sample[which(Seurat_test$Sample != "GW9_GE")]

# save(Seurat_test, file = "../../data/processed/GW9101318192539_MGE_all_raw.RData")