## Figure 4

#### Figure 4_A

DSP.seurat <- AddModuleScore(DSP.seurat,
                             features = Lineage_Score_list,
                             name = "Lineage_Score")
colnames(DSP.seurat@meta.data)[grep("Lineage_Score", colnames(DSP.seurat@meta.data))] = names(Lineage_Score_list)
DSP.seurat = subset(DSP.seurat, celltype == "ROI_DCX")
scatterpie_df = DSP.seurat@meta.data[, c("ROICoordinateX", "ROICoordinateY", names(Lineage_Score_list))]

write.xlsx(Lineage_Score_list, file = "../../results/tables/Figure4_H_ScoreGene.xlsx")

#### Figure 4_F

Express_average = (AverageExpression(ref, assays = "SCT", 
                                     features = union(gene_selected, rownames(ref@assays$SCT@data)), 
                                     group.by = "Revision_Lineage_New"))[["SCT"]]
df = as.data.frame(Express_average)
df$GeneType = "Common";df$GeneType[which(rownames(df) %in% gene_selected)] = "Variable"
df$SumUp = df$CRABP1ANGPT2 + df$EPHA5MEF2C + df$LHX6NFIA;df[, c(1,2,3)] = df[, c(1,2,3)]/df$SumUp
pannel_loosed = rownames(df)[which(MatrixGenerics::rowMins(df[, c(1,2,3)] %>% as.matrix()) <= 0.2)]
df$gene = rownames(df);df$lab = NA;df$lab[which(df$GeneType == "Variable")] = df$gene[which(df$GeneType == "Variable")]
df$GeneType[which(df$gene %in% pannel_loosed)] = "Variable"

save(df, file = "../../results/Figure4_F.RData")

#### Figure 4_G

df <- melt(scatterpie_df[, c(3:7)], id.vars=c("Node", "type"),
           measure.vars = c("CRABP1ANGPT2", "LHX6NFIA", "EPHA5MEF2C"), 
           # variable.name = "type", 
           # value.name = "x"
)
df$variable = factor(df$variable, levels = c("EPHA5MEF2C", "LHX6NFIA", "CRABP1ANGPT2"));df$type = as.factor(df$type)

save(df, file = "../../results/Figure4_G.RData")

## Figure 5

#### Figure 5_E~G

setname = "DynamicGene"
Seurat_test_all <- AddModuleScore(Seurat_test_all,
                                  features = gene_SET1,
                                  name = c("EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score"))
colnames(Seurat_test_all@meta.data)[grep(pattern = "EPHA5MEF2CScore|LHX6NFIAScore|CRABP1ANGPT2Score", x = colnames(Seurat_test_all@meta.data))] = c("EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score")

Seurat_test_all <- AddModuleScore(Seurat_test_all,
                                  features = gene_SET2,
                                  name = c("SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score"))
colnames(Seurat_test_all@meta.data)[grep(pattern = "SVZFateScore|NR2F1NR2F2Score|LHX8ISL1Score", x = colnames(Seurat_test_all@meta.data))] = c("SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score")

Seurat_test_all$test_cluster = Seurat_test_all$Revision_LineageType_New
df = Seurat_test_all@meta.data[, c("test_cluster", "EPHA5MEF2CScore", "LHX6NFIAScore", "CRABP1ANGPT2Score", "SVZFateScore", "NR2F1NR2F2Score", "LHX8ISL1Score")]

# library(URD)
DefaultAssay(Seurat_test) = "RNA"
Seurat_test = NormalizeData(Seurat_test) %>% FindVariableFeatures() %>% ScaleData()
Seurat_test <- Seurat_test %>% SCTransform(verbose = T, vst.flavor ='v2')
Seurat_test$DotType = Seurat_test$predicted.LineageType
table(Seurat_test$DotType)

setname = "VZ_RGC_Fate" # setname = "SVZ_RGC_Fate" 
table(Seurat_test$DotType)
Seurat_test$test_cluster = Seurat_test$DotType
Seurat_test_loop = Seurat_reset(Seurat_test, meta_remain = T)

Idents(Seurat_test_loop) = "test_cluster"
Seurat_test_loop$URD_Type = Idents(Seurat_test_loop)

testsample <- createURD(count.data = Seurat_test_loop@assays$RNA@counts, meta = Seurat_test_loop@meta.data, min.cells=3, min.counts=3)
table(testsample@meta$URD_Type)

testsample@group.ids$stage <- as.character(testsample@meta[rownames(testsample@group.ids), "Age"])
stages <- paste0("GW", c(9, 10, 13, 18, 19, 26, 39));stages

var.by.stage <- lapply(seq(3,6,1), function(n) {
  findVariableGenes(testsample, cells.fit = cellsInCluster(testsample, "stage", stages[(n-2):n]),
                    set.object.var.genes=F, diffCV.cutoff=0.3,
                    mean.min=.005, mean.max=100,
                    main.use=paste0("Stages ", stages[n-2], " to ", stages[n]), do.plot=T)
})
var.genes <- sort(unique(unlist(var.by.stage)));testsample@var.genes <- var.genes;testsample <- calcPCA(testsample, mp.factor = 2)
pcSDPlot(testsample)
testsample <- calcTsne(object = testsample)
plotDim(testsample, "Age", plot.title = "tSNE: Stage");plotDim(testsample, "URD_Type", plot.title = "tSNE: Stage")

knn_selected = round((dim(Seurat_test_loop@assays$RNA)[2])^0.5);knn_selected
testsample <- calcDM(testsample, knn = knn_selected, sigma = 16)
plotDimArray(testsample, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 16, NNs): Stage", label="stage", plot.title="", legend=F)
plotDim(testsample, "URD_Type", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

# Here we use all cells from the first stage as the root
root.cells <- rownames(testsample@meta)[which( (testsample@meta$URD_Type == "VZ_RGC(Homogeneous)") & (testsample@meta$Age == "GW9") )]

# Then we run 'flood' simulations
testsample.floods <- floodPseudotime(testsample, root.cells = root.cells, n = 40, minimum.cells.flooded = 2, verbose=T)
# The we process the simulations into a pseudotime
testsample <- floodPseudotimeProcess(testsample, testsample.floods, floods.name="pseudotime")
testsample@group.ids$URD_Type = testsample@meta$URD_Type
cells_outs = rownames(testsample@meta)[is.na(testsample@pseudotime$pseudotime)]
testsample <- urdSubset(testsample, cells.keep = setdiff(colnames(testsample@count.data), cells_outs))
pseudotimePlotStabilityOverall(testsample)
plotDim(testsample, "pseudotime");plotDists(testsample, "pseudotime", "stage", plot.title="Pseudotime by stage")

testsample@group.ids$URD_num = as.numeric(as.factor(testsample@group.ids$URD_Type))
# Create a subsetted object of just those cells from the final stage
id2 = cellsInCluster(testsample, "URD_Type", c("VZRGC_SVZFate", "VZRGC_LHX8ISL1", "VZRGC_NR2F1NR2F2"))
testsample.final <- urdSubset(testsample, cells.keep = id2)
table(testsample.final@meta$URD_Type)
# Use the variable genes that were calculated only on the final group of stages (which
# contain the last stage).
# testsample.final@var.genes <- var.by.stage[[4]]
testsample.final@var.genes <- var.by.stage[[2]]
# Calculate PCA and tSNE
testsample.final <- calcPCA(testsample.final, mp.factor = 1.5)
pcSDPlot(testsample.final)
testsample.final <- calcTsne(testsample.final)
# Calculate graph clustering of these cells
testsample.final <- graphClustering(testsample.final, num.nn = 50, do.jaccard=T, method="Louvain")
# Copy cluster identities from axial.6somite object to a new clustering ("tip.clusters") in the full axial object.
# testsample@group.ids$tip.clusters=NULL
testsample@group.ids[rownames(testsample.final@group.ids), "tip.clusters"] <- testsample.final@group.ids$URD_num
table(testsample@group.ids$tip.clusters)
# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
testsample.ptlogistic <- pseudotimeDetermineLogistic(testsample, "pseudotime", optimal.cells.forward=10, max.cells.back = 10, do.plot = T)
# testsample.ptlogistic <- pseudotimeDetermineLogistic(testsample, "pseudotime", optimal.cells.forward = 200, max.cells.back = 100, do.plot = T)
# Bias the transition matrix acording to pseudotime
testsample.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(testsample, "pseudotime", logistic.params=testsample.ptlogistic))
dim(testsample.biased.tm)
# Simulate the biased random walks from each tip
testsample.walks <- simulateRandomWalksFromTips(testsample,
                                                tip.group.id = "tip.clusters",
                                                root.cells = root.cells,
                                                transition.matrix = testsample.biased.tm,
                                                # n.per.tip = 1000, root.visits = 1, verbose = T
                                                max.steps = ncol(testsample@logupx.data) * 1
)
# Process the biased random walks into visitation frequencies
testsample <- processRandomWalksFromTips(testsample, testsample.walks, verbose = F)
# Load the cells used for each tip into the URD object
testsample.tree <- loadTipCells(testsample, "tip.clusters")
# Build the tree
testsample.tree <- buildTree(testsample.tree, pseudotime = "pseudotime",
                             divergence.method = "preference", cells.per.pseudotime.bin = 20,#5,
                             bins.per.pseudotime.window = 10,#2,
                             save.all.breakpoint.info = T, p.thresh=0.05)
pseudotime_df = testsample.tree@pseudotime
plotTree(testsample.tree, "Age", title=paste0("Developmental Stage of ", setname))
plotTree(testsample.tree, "URD_Type", title=paste0("Developmental Stage of ", setname))

#
change_table = c("SVZRGC_EPHA5MEF2C", "SVZRGC_CRABP1ANGPT2", "SVZRGC_LHX6NFIA")
change_table = c("VZRGC_SVZFate", "VZRGC_LHX8ISL1", "VZRGC_NR2F1NR2F2")

Seurat_test_loop = AddMetaData(Seurat_test_loop, metadata = testsample.tree@group.ids)
meta_add = testsample.tree@pseudotime
# colnames(meta_add)[2:4] = c("SVZRGC_EPHA5MEF2C", "SVZRGC_CRABP1ANGPT2", "SVZRGC_LHX6NFIA")
colnames(meta_add)[2:4] = c("VZRGC_SVZFate", "VZRGC_LHX8ISL1", "VZRGC_NR2F1NR2F2")
# colnames(meta_add)[c(2:5)] = change_table[colnames(meta_add)[c(2:5)]]
colnames(meta_add) = paste0("URD_", colnames(meta_add))
Seurat_test_loop = AddMetaData(Seurat_test_loop, metadata = meta_add, col.name = colnames(meta_add))

Idents(Seurat_test_loop) = "URD_Type"

gc()
# subline = change_table[1];subline
# subline = change_table[2];subline
subline = change_table[3];subline

Seurat_test_Line = subset(Seurat_test_loop, URD_Type %in% setdiff(change_table, subline), invert = T)
# cells_out = names(which(is.na(Seurat_test_Line$URD_SVZRGC_EPHA5MEF2C)))
# cells_out = names(which(is.na(Seurat_test_Line$URD_SVZRGC_CRABP1ANGPT2)))
# cells_out = names(which(is.na(Seurat_test_Line$URD_SVZRGC_LHX6NFIA)))
# cells_out = names(which(is.na(Seurat_test_Line$URD_VZRGC_SVZFate)))
# cells_out = names(which(is.na(Seurat_test_Line$URD_VZRGC_LHX8ISL1)))
cells_out = names(which(is.na(Seurat_test_Line$URD_VZRGC_NR2F1NR2F2)))
Seurat_test_Line = subset(Seurat_test_Line, cells = cells_out, invert = T)
dim(Seurat_test_Line)
Seurat_test_Line$URD_Ptime = Seurat_test_Line$URD_VZRGC_NR2F1NR2F2

quantile(Seurat_test_Line$URD_Ptime, seq(0, 1, length.out = 20))
for (i in unique(Seurat_test_Line$URD_Type)) {
  print(i)
  print(quantile(Seurat_test_Line$URD_Ptime[which(Seurat_test_Line$URD_Type == i)], seq(0, 1, length.out = 20)))
}
Seurat_test_Line$URD_Ptime = max(Seurat_test_Line$URD_Ptime) - Seurat_test_Line$URD_Ptime

Seurat_test_Line = FindVariableFeatures(Seurat_test_Line, nfeatures = 3000);top_hvg <- VariableFeatures(Seurat_test_Line)
# Prepare data for random forest
dat_use <- t(Seurat_test_Line@assays$RNA@data[intersect(top_hvg, rownames(Seurat_test_Line@assays$RNA@data)), ]);dim(dat_use)
dat_use_df <- cbind(Seurat_test_Line$URD_Ptime, dat_use);colnames(dat_use_df)[1] <- "URD_Ptime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[, 1]), ])

library(rsample);library(parsnip);library(tradeSeq)
dat_split <- initial_split(dat_use_df);dat_train <- training(dat_split);dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 50) %>%
  parsnip::fit(URD_Ptime ~ ., data = dat_train)
val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[, -1]) %>% pull()) %>% 
  select(truth = URD_Ptime, estimate)
model_outs = yardstick::metrics(data = val_results, truth, estimate)
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
importance_table = as.data.frame(var_imp);importance_table$gene = rownames(importance_table)

pseudotime_fake = matrix(data = Seurat_test_Line$URD_Ptime, nrow = length(Seurat_test_Line$URD_Ptime), ncol = 1, 
                         dimnames = list(rownames(Seurat_test_Line@meta.data), "curve1"))
# pseudotime_fake[which(is.na(pseudotime_fake))] = 0
cellWeights_fake = matrix(data = rep(1, length(Seurat_test_Line$URD_Ptime)), nrow = length(Seurat_test_Line$URD_Ptime), ncol = 1, 
                          dimnames = list(rownames(Seurat_test_Line@meta.data), "curve1"))
sce <- fitGAM(counts = Seurat_test_Line@assays$RNA@data[importance_table$gene[1:2500], ], pseudotime = pseudotime_fake, parallel = T,
              cellWeights = cellWeights_fake,
              nknots = 6, verbose = FALSE)
ATres <- associationTest(sce);topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:500]
pst.ord <- order(sce$crv$pseudotime.curve1, na.last = NA)

DefaultAssay(Seurat_test_Line) = "RNA"
Seurat_test_Line = Seurat_test_Line %>% FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
cds_subset <- as.CellDataSet(Seurat_test_Line)
newdata <- data.frame(URD_Ptime = seq(
  min(Seurat_test_Line$URD_Ptime),
  max(Seurat_test_Line$URD_Ptime),
  length.out = 100))
library(monocle)
diff_test_result <- monocle::differentialGeneTest(cds_subset[rownames(ATres[order(ATres$pvalue), ])[1:1500]],
                                                  fullModelFormulaStr = "~sm.ns(URD_Ptime, df=3)",
                                                  reducedModelFormulaStr = "~1", relative_expr = F, cores = 16,
                                                  verbose = FALSE)

