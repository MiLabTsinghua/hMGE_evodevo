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

# ------ 初始化环境 ------
# rm(list = ls())
options(stringsAsFactors = FALSE)
# set.seed(2024)

# 加载必要包
# require_package <- function(pkg) {
#   if (!require(pkg, character.only = TRUE)) {
#     stop("缺少必需依赖包: ", pkg, "\n使用install.packages('", pkg, "')安装")
#   }
# }
# require_package("tidyverse")
# require_package("here")  # 使用here包管理路径

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
#palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(3)  # 高级配色方案

library(org.Hs.eg.db) # human的OrgDB
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



















