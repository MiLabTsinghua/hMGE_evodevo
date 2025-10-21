# -*- coding: utf-8 -*-
###############################################################################
# File: 
# Project: Human GE ()
# Author: Yiming Yan
# Email: yym21@mails.tsinghua.edu.cn
# Institution: Tsinghua Univ IDG/McGovern Institute for Brain Research
# Created: 2024-05-01
# Last modified: 2025-03-19
# License: MIT
# Version: 0.1.1
#-------------------------------------------------------------------------------

###############################################################################

HeatCol <- colorRampPalette(colors = c('#176BA0','black','#EE9A3A'))(1000)

MajorType_Col = c("#F0BCCA", "#00BFFF", "gray70");names(MajorType_Col) = c(
  "Progenitor", "Neuron", "Others"
)

# MGE_cluster_color <- c("red", "#176BA0", "#FF8C00",  "#F0BCCA", "#F0BCCA",
#                        "#00BFFF", "#D4F2E7", "#FBC241", "bisque4", "#f056b4", "grey70",
#                        "red", "#9370db", "red", "#9370db", "#176BA0","#9370db"
#                        )
# names(MGE_cluster_color) = c("RGC", "APC(Astrocyte)", "IPC", "Progenitor", "Neuroblast", 
#                              "MGE_derived_Neuron", "OPC", "Microglia", "Pericyte", "Endothelium", "Others",
#                              "VZ_RGC", "non-VZ_RGC", "VZRGC", "non-VZRGC", "Astrocyte", "SVZ_RGC"
#                              )

MGE_cluster_color <- c("red", "#176BA0", "#FF8C00",  "#F0BCCA", "#F0BCCA",
                       "#00BFFF", "#D4F2E7", "#FBC241", "bisque4", "#f056b4", "grey70",
                       "red", "#9370db", "red", "#9370db", "#176BA0","#9370db",
                       "brown", "darkorchid3",
                       "red", "#9370db"
)
names(MGE_cluster_color) = c("RGC", "APC(Astrocyte)", "IPC", "Progenitor", "Neuroblast", 
                             "MGE_derived_Neuron", "OPC", "Microglia", "Pericyte", "Endothelium", "Others",
                             "VZ_RGC", "non-VZ_RGC", "VZRGC", "non-VZRGC", "Astrocyte", "SVZ_RGC",
                             "VZ_Progenitor", "SVZ_Progenitor",
                             "RGC_1", "RGC_2"
)


# MGE_progenitor_color <- c("red", "#9370db", "#FF8C00",  "#F0BCCA", "#00BFFF")
# names(MGE_progenitor_color) = c("VZ_RGC", "SVZ_RGC", "IPC", "Neuroblast", "Neuron")

MGE_progenitor_color <- c("red", "#9370db", "#FFCB0A",  "#F0BCCA", "#00BFFF")
names(MGE_progenitor_color) = c("VZ_RGC", "SVZ_RGC", "IPC", "Neuroblast", "Neuron")


MGE_progenitor_color2 <- c(
  "#FF8C00", "#6943CD", "#33cc33", "#808080",
  "#FF9FED","#C700C9", "#B4FBFF")
names(MGE_progenitor_color2) = c(
  "IPC(Neuronal)", "SVZ_RGC(Neuronal)", "VZ_RGC(Heterogeneous)", "VZ_RGC(Homogeneous)",
  "VZ_RGC(Astrocyte)", "SVZ_RGC(Astrocyte)", "SVZ_RGC(OPC)"
)

method_selected = c("integrated.cca", "integrated.rpca", "harmony", "integrated.mnn");umap_selected = c("umap.unintegrated", "umap.integrated.cca", "umap.integrated.rpca", "umap.harmony", "umap.integrated.mnn")


MGE_branch_color = c('#FF0000FF', "#0000FFFF", "#00EDEDFF", "#009404FF", "#FFD600FF",
                     '#FF0000FF', "#0000FFFF", "#00EDEDFF", "#009404FF", "#FFD600FF",
                     '#FF0000FF', "#0000FFFF", "#00EDEDFF", "#009404FF", "#FFD600FF",
                     '#FF0000FF', "#0000FFFF", "#00EDEDFF", "#009404FF", "#FFD600FF",
                     "#00597F", "#808080")
names(MGE_branch_color) = c('Neuron_EPHA5MEF2C', 'Neuron_LHX6NFIA', 'Neuron_CRABP1ANGPT2',
                            'Neuron_NR2F1NR2F2', 'Neuron_LHX8ISL1',
                            'IPC_EPHA5MEF2C', 'IPC_LHX6NFIA', 'IPC_CRABP1ANGPT2',
                            'IPC_NR2F1NR2F2', 'IPC_LHX8ISL1',
                            'EPHA5MEF2C', 'LHX6NFIA', 'CRABP1ANGPT2',
                            'NR2F1NR2F2', 'LHX8ISL1',
                            'SVZRGC_EPHA5MEF2C', 'SVZRGC_LHX6NFIA', 'SVZRGC_CRABP1ANGPT2',
                            'VZRGC_NR2F1NR2F2', 'VZRGC_LHX8ISL1',
                            "VZRGC_SVZFate", "VZ_RGC(Homogeneous)")

# MGE_cluster_color <- c("red", "#176BA0", "#FF8C00",  "#F0BCCA", "#F0BCCA",
#                        "#00BFFF", "#D4F2E7", "#FBC241", "bisque4", "#f056b4", "grey70",
#                        "red", "#9370db", "red", "#9370db"
# )
# names(MGE_cluster_color) = c("RGC", "APC(Astrocyte)", "IPC", "Progenitor", "Neuroblast",
#                              "MGE_derived_Neuron", "OPC", "Microglia", "Pericyte", "Endothelium", "Others",
#                              "VZ_RGC", "non-VZ_RGC", "VZRGC", "non-VZRGC"
# )

Batch_color <- c("#FBC241", "#E64B35FF", "#4DBBD5FF", "#4DBBD8FF", "#00A087FF", "#00A023FF", "#3C5488FF", "#F39B7FFF");names(Batch_color) = c("GW9_XiaoqunW", "GW10_MGE", "GW13_MGE", "GW13_MGE_rep2", "GW18_MGE", "GW19_MGE", "GW26_MGE", "GW39_MGE")

Speices_color = c(
  "Mouse_RGP_1" = "#D4F2E7", "Mouse_RGP_2" = "#FBC241", "Mouse_RGP_3" = "#00BFFF",
  "Human_VZ_RGC" = "#FF0000", "Human_SVZ_RGC" = "#9370db",
  "Macaque_RG_HOPX_NRG1" = "#F0BCCA", "Macaque_RG_STMN2_SOX5" = "#29AAE1", "Macaque_RG_CRYAB" = "#662D90")
