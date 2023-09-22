# Goal: In silico CD45 sorting
# Methods: 1) define CD45+/- identify highly variable genes; 

# rm(list=ls())

library(dplyr)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)

# make heatmap for one gene and multiple peaks
make_heatmap_rna_and_atac <- function(obj,gene_name=NULL,peaks=NULL){

    prepare_rna_data <- function(obj,gene_name=NULL){

        rna_data  <- as.matrix(obj$RNA@data)
        rna_data  <- rna_data[gene_name,,drop=F]; print(dim(rna_data))
        rna_data

    }

    prepare_atac_data <- function(obj,peaks=NULL){

        atac_data  <- as.matrix(obj$ATAC@data)
        peaks      <- peaks[peaks %in% rownames(atac_data)]
        atac_data  <- atac_data[peaks,,drop=F]; print(dim(atac_data))
        atac_data

    }

    rna_data  <- prepare_rna_data(obj,gene_name)
    atac_data <- prepare_atac_data(obj,peaks)

    # clustering cells
    c <- hclust(dist(t(rna_data)))
    rna_data <- rna_data[,c$order,drop=F]
    atac_data <- atac_data[,c$order,drop=F]

    atac_cell_mean <- t(as.matrix(apply(atac_data,2,function(x) mean(x,na.rm=T))))
    rownames(atac_cell_mean) <- 'Mean ATAC'


    cell_info <- obj@meta.data

    
    
    ########################################################################################
    # make heatmap 
    ########################################################################################
    print("make heatmap")

    # make heatmap annotations
    make_heatmap_annotations <- function(cell_info,rna_data){

        # // prepare cell clusters annotations for heatmap
        print("prepare cell clusters annotations for heatmap")
        cell_clusters     <- sort(unique(cell_info$seurat_clusters))

        cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                              "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                              "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                              "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
        cluster_cols             = cluster_cols[1:length(cell_clusters)]
        names(cluster_cols)      = cell_clusters

        cols=list('Cluster' = cluster_cols)


        col_ha = HeatmapAnnotation(
                                  'Cluster'=cell_info[match(colnames(rna_data),cell_info$cell),]$seurat_clusters,
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')


        annotation_list <- list(col_ha)

    }


    # make heatmap annotations
    annotation_list <- make_heatmap_annotations(cell_info,rna_data)
    col_ha          <- annotation_list[[1]]

    
    p1 <- Heatmap(rna_data,name="RNA",show_column_names=F,column_split = cell_info[match(colnames(rna_data),cell_info$cell),]$seurat_clusters,column_title=paste0(ncol(rna_data)," cells"),top_annotation = col_ha,cluster_columns = F,show_column_dend=F,row_names_side='left',cluster_rows=F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),row_title_side='right',width=unit(8,"inch"))
    p2 <- Heatmap(atac_cell_mean,name="Mean ATAC",show_column_names=F,column_split = cell_info[match(colnames(rna_data),cell_info$cell),]$seurat_clusters,column_title=paste0(ncol(rna_data)," cells"),top_annotation = NULL,cluster_columns = F,show_column_dend=F,row_names_side='left',cluster_rows=F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),row_title_side='right',width=unit(8,"inch"))
    p3 <- Heatmap(atac_data,name="ATAC",show_column_names=F,column_split = cell_info[match(colnames(atac_data),cell_info$cell),]$seurat_clusters,column_title=paste0(ncol(atac_data)," cells"),top_annotation = NULL,cluster_columns = F,show_column_dend=F,row_names_side='left',cluster_rows=F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),row_title_side='right',width=unit(8,"inch"))

    p1 %v% p2 %v% p3

}