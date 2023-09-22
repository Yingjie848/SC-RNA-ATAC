# calculate gene signature score by taking average log2-CPM (counts per million) for each cell

calc_cpm <- function(counts){

    a = colSums(counts) # get sums from each column
    mat = as.matrix(counts)
    c = sweep(mat, 2, a, '/') # divide the matrix by the vector
    cpm = c * 1e6 # multiple 1 million to form cpm

}

calc_sigscore <- function(counts,genes=NULL){

    cpm <- calc_cpm(counts)/100

    cpm_siggenes <- cpm[rownames(cpm) %in% genes,]

    sigscore <- colMeans(cpm_siggenes,na.rm=T)

}

main <- function(){

    file_seurat="seurat_cellRangerData_emptyDrops_scDblFinder_cells_results/reclustered_rmC4C12/reclustered_tumors_C0_1_2_4_5_7/seurat_rna.rds"
    obj      <- readRDS(file_seurat)
    counts   <- obj$RNA@data
    metaData <- obj@meta.data

    metaData$chondroid_sigscore <- calc_sigscore(counts,genes)

    library(ggpubr)
    ggdensity(metaData,x='chondroid_sigscore',color='seurat_clusters')

}