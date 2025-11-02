## Spatial dataform pipeline
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(SeuratDisk))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(RColorBrewer))

gem_to_seuratObject <- function(data, prefix = 'sample', binsize = 1){
    #' group counts into bins
    data$x <- trunc(data$x / binsize) * binsize
    data$y <- trunc(data$y / binsize) * binsize
    #trunc取整函数，只取整数部分
    #将横纵坐标根据binsize取整，比如，如果binszie=10,x=2401,y=4522,
    #那么处理后，x=2400,y=4520
    
    if ('MIDCount' %in% colnames(data)) {
      # aa %in% bb，判断aa是否存在于bb中
      data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]
    } else {
      data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
    }
    
    #' create sparse matrix from stereo
    data$cell <- paste0(prefix, ':', data$x, '_', data$y)
    data$geneIdx <- match(data$geneID, unique(data$geneID))
    data$cellIdx <- match(data$cell, unique(data$cell))
    #match(aa,bb)返回bb在aa中的位置
    
    # if (! is.null(opts$binsize)){
    #   #is.null(aa) 判断aa是否为空，如果空，返回Ture；
    #   #! is.null(aa)，相反地，判断aa是否非空，如果非空，返回True；
    #   write.table(data, file = paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '.tsv'), 
    #               quote = FALSE, sep = '\t', row.names = FALSE)
    # }
    #write.table(aa,file=xx,sep =" ",quote=TRUE,row.names =TRUE, col.names =TRUE)
    #把内容aa输出到文件xx中
    
    mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                        dimnames = list(unique(data$geneID), unique(data$cell)))
    #sparseMatrix稀疏矩阵函数
    
    cell_coords <- unique(data[, c('cell', 'x', 'y')])
    #unique去除重复函数，删除cell,x和y都一样的行
    
    rownames(cell_coords) <- cell_coords$cell
    
    #cell_coords$cell <- NULL
    
    seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial', 
                                            names.delim = ':', meta.data = cell_coords)
    
    
    #' create pseudo image
    cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
    cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
    
    tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
    #matrix(aa,x,y)以aa为输入向量，创建一个x行y列的矩阵
    #构造一个seruat image
    
    tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                        tissue = 1,
                                        row = cell_coords$y, col = cell_coords$x,
                                        imagerow = cell_coords$y, imagecol = cell_coords$x)
    
    
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                     tissue_hires_scalef = 1,
                                     tissue_lowres_scalef = 1))
    #toJSON: 把json格式 转换成 list格式
    
    
    #' function to create image object
    generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
      if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
      }
      
      unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
      
      spot.radius <- unnormalized.radius / max(dim(x = image))
      
      return(new(Class = 'VisiumV1', 
                 image = image, 
                 scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                              fiducial = scale.factors$fiducial_diameter_fullres, 
                                              hires = scale.factors$tissue_hires_scalef, 
                                              lowres = scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius))
    }
    
    spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                      scale.factors = fromJSON(scalefactors_json), 
                                      tissue.positions = tissue_positions_list)
    #可以理解为构建一个spatial背景
    
    #' import image into seurat object
    spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
    DefaultAssay(spatialObj) <- 'Spatial'
    
    seurat_spatialObj[['slice1']] <- spatialObj
    rm("spatialObj")
    rm("data")
    rm("mat")
    
    return(seurat_spatialObj)
    }

# Example
RDS_output_path = "../../data/processed/Spatial_Data/"
gem_name = "D01873D3.bin20.GW18.MGE"
data = fread('D01873D3.lasso.gem.gz') # download from https://ngdc.cncb.ac.cn/omix/
Seurat_stereo <- gem_to_seuratObject(data, binsize = 20, prefix = gem_name)
Seurat_stereo = UpdateSeuratObject(Seurat_stereo)
save(Seurat_stereo, file = paste0(RDS_output_path, gem_name, ".RData"))


## Spatial mapping pipeline

Seurat_test # SC dataset
Seurat_stereo # Stereoseq dataset

Loc = Seurat_stereo@meta.data[, c("x", "y")];colnames(Loc) = paste0("Loc_", c(1, 2))
Seurat_stereo = Add_DimReduc_Seurat(Seurat_stereo, Loc, Reduc_name = "RealLocation", assay_loc = "Spatial")

Age_loop = "GW26";Age_loop2 = "GW26"

filepath_out = "../../data/processed/Spatial_Mapping"
if(!dir.exists(filepath_out)){dir.create(filepath_out);print("created")}else{print("exsited")}
filepath_add = paste0(filepath_out, "/", Age_loop)
if(!dir.exists(filepath_add)){dir.create(filepath_add);print("created")}else{print("exsited")}

setname = paste0(Age_loop, "_All_Bin20");setname

Seurat_test$CellSubtype = Idents(Seurat_test)
ref = subset(Seurat_test, Age == Age_loop2);Idents(ref) = "CellSubtype";ref = subset(ref, idents = names(which(table(Idents(ref)) > 15)))
names(which(table(Idents(ref)) > 15))
Seurat_stereo = RCTD_test(ref, Seurat_stereo)

Idents(Seurat_stereo) = "first_type"
Seurat_stereo$SolidType = "";Seurat_stereo$SolidType[which(Seurat_stereo$spot_class == "singlet")] = as.character(Seurat_stereo$first_type)[which(Seurat_stereo$spot_class == "singlet")]
Idents(Seurat_stereo) = 'SolidType'
table(Idents(Seurat_stereo))

# For example VZ RGC subtype
Seurat_test  # SC dataset RGC subtype idents
Seurat_test = Seurat_Object_Reset(Seurat_test, if_SCT = F)

ref = subset(Seurat_test, idents = "VZ_RGC")
ref = subset(ref, Age %in% Age_loop2)
Idents(ref) = "RGCSubtype"

setname_sub = paste0(Age_loop, "_VZRGC_Bin20");setname_sub
Seurat_stereo_sub = subset(Seurat_stereo, SolidSubtype_RGC == "VZ_RGC")
Seurat_stereo_sub = RCTD_test(ref, Seurat_stereo_sub)

Seurat_stereo_sub$SolidSubtype_VZRGC = "";Seurat_stereo_sub$SolidSubtype_VZRGC[which(Seurat_stereo_sub$spot_class == "singlet")] = as.character(Seurat_stereo_sub$first_type)[which(Seurat_stereo_sub$spot_class == "singlet")]
Seurat_stereo_VZRGC = Seurat_stereo_sub
Seurat_stereo$SolidSubtype_VZRGC = "";Seurat_stereo$SolidSubtype_VZRGC = Idents_merged(Seurat_stereo$SolidSubtype_VZRGC, Seurat_stereo_VZRGC$SolidSubtype_VZRGC)

# For SVZ RGC IPC Neuron etc. Same as VZ RGC
