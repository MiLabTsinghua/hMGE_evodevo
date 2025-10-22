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

# ------ Initialize environment ------
# rm(list = ls())
options(stringsAsFactors = FALSE)
# set.seed(2024)

# Load required packages
# require_package <- function(pkg) {
#   if (!require(pkg, character.only = TRUE)) {
#     stop("Required dependency package missing: ", pkg, "\nInstall using install.packages('", pkg, "')")
#   }
# }
# require_package("tidyverse")
# require_package("here")  # Use here package for path management

require(Seurat)
require(clustree)

require(pheatmap)

require(dbplyr)
require(dplyr)
require(reshape2)
require(scales)
require(tidyverse)

require(ggplot2)
require(cowplot)
require(patchwork)
require(ggsci)
require(ggnewscale)
require(gghalves)

require(scRNAtoolVis)
require(scCustomize)

require(SeuratDisk)


require(rhdf5)
require(Matrix)

require(RColorBrewer)
#palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(3)  # Advanced color scheme

library(org.Hs.eg.db) # human OrgDB
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)



library(hdWGCNA)
library(openxlsx)
library(ggthemes)
library(scPred)

require(MetaNeighbor)


require(future)
availableCores()
plan("multisession", workers = 64)
nbrOfWorkers()
options(future.globals.maxSize= (10000*1024^2) )



















