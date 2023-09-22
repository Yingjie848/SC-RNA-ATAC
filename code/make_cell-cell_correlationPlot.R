# To understand relationship between cells, for example, one cluster could be mixture of cells from all other cells

library(ComplexHeatmap)

#file_rdata <- "MBC4/seurat_cellRangerData_emptyDrops_scDblFinder_cells_results/InSilicoCD45sorting_DEGs/data.Rdata"
#outDir <- "MBC4/seurat_cellRangerData_emptyDrops_scDblFinder_cells_results/InSilicoCD45sorting_DEGs/"

make_cell_cell_corr_plot <- function(file_rdata,outDir){

    rm(list=ls())

    load(file_rdata)

    ########################################################################################
    # make cell-cell correlation heatmap
    ########################################################################################

    # make heatmap column annotations
    make_heatmap_annotations <- function(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg){

        # // prepare cell clusters annotations for heatmap
        print("prepare cell clusters annotations for heatmap")
        cell_clusters     <- unique(as.character(sort(as.numeric(as.vector(cell_info$seurat_clusters)))))

        cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                              "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                              "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                              "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")

        #ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        #  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        #  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
        #}
        #cluster_cols <- ggplotColours(n=length(cell_clusters))
        cluster_cols             = cluster_cols[1:length(cell_clusters)]
        names(cluster_cols)      = cell_clusters

        library(circlize)
        CD45_cols            = colorRamp2( c(min(log10(cell_info$CD45+1)), max(log10(cell_info$CD45+1))), c("#f5f2f2", "red"))
        immune_gsea_cols     = colorRamp2( c(min(cell_info$immune_gsea), max(cell_info$immune_gsea)), c("#f5f2f2", "red"))
        CellTotalCounts_cols = colorRamp2( c(0, max(log10(cell_info$nCount_RNA+1))), c("#f5f2f2", "#ff00a6"))
        CellTotalGenes_cols  = colorRamp2( c(0, max(log10(cell_info$nFeature_RNA+1))), c("#f5f2f2", "#01c504"))
        panck_cols           = colorRamp2( c(0, max(log10(cell_info$PanCK+1))), c("white", "black"))
        dbl_cols             = colorRamp2( c(0, max(cell_info$scDblFinder.score)), c("white", "purple"))
        scDbl_cols           = c('singlet'='white','doublet'='purple')
        hybrid_cols          = colorRamp2( c(0, max(cell_info$hybrid_score)), c("white", "darkgreen"))
        doublet_cols         = c('No'='white','Yes'='purple')

        cols=list('Seurat Cluster' = cluster_cols,
                  'CD45 counts (log10)'=CD45_cols,
                  'Immune GSEA'=immune_gsea_cols,
                  'Cell Total Counts (log10)'=CellTotalCounts_cols,
                  'Cell Total Genes (log10)'=CellTotalGenes_cols,
                  'Pan-CK (log10)'=panck_cols,
                  'Doublets (scDblFinder)'=scDbl_cols,
                  'Doublets score (scDblFinder)'=dbl_cols,
                  'Doublets score (scds)'=hybrid_cols,
                  'Doublets'=doublet_cols)


        col_ha = HeatmapAnnotation(
                                  'Seurat Cluster'=cell_info[match(colnames(SCT_data),cell_info$cell),]$seurat_clusters,
                                  'CD45 counts (log10)'=log10(cell_info$CD45+1),
                                  'Immune GSEA'=cell_info$immune_gsea,
                                  'Pan-CK (log10)'=log10(cell_info$PanCK),
                                  'Cell Total Counts (log10)'=log10(cell_info$nCount_RNA),
                                  'Cell Total Genes (log10)'=log10(cell_info$nFeature_RNA),
                                  'Doublets (scDblFinder)'=cell_info$scDblFinder.doublets,
                                  'Doublets score (scDblFinder)'=cell_info$scDblFinder.score,
                                  'Doublets score (scds)'=cell_info$hybrid_score,
                                  # 'Doublets'=cell_info$Doublets,
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')

        col_ha_CD45pos = HeatmapAnnotation(
                                  'Seurat Cluster'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$seurat_clusters,
                                  'CD45 counts (log10)'=log10(cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$CD45+1),
                                  'Immune GSEA'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$immune_gsea,
                                  'Pan-CK (log10)'=log10(cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$PanCK),
                                  'Cell Total Counts (log10)'=log10(cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$nCount_RNA),
                                  'Cell Total Genes (log10)'=log10(cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$nFeature_RNA),
                                  'Doublets (scDblFinder)'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$scDblFinder.doublets,
                                  'Doublets score (scDblFinder)'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$scDblFinder.score,
                                  'Doublets score (scds)'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$hybrid_score,
                                  # 'Doublets'=cell_info[match(colnames(SCT_data_CD45pos),cell_info$cell),]$Doublets,
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')

        col_ha_CD45neg = HeatmapAnnotation(
                                  'Seurat Cluster'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$seurat_clusters,
                                  'CD45 counts (log10)'=log10(cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$CD45+1),
                                  'Immune GSEA'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$immune_gsea,
                                  'Pan-CK (log10)'=log10(cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$PanCK),
                                  'Cell Total Counts (log10)'=log10(cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$nCount_RNA),
                                  'Cell Total Genes (log10)'=log10(cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$nFeature_RNA),
                                  'Doublets (scDblFinder)'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$scDblFinder.doublets,
                                  'Doublets score (scDblFinder)'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$scDblFinder.score,
                                  'Doublets score (scds)'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$hybrid_score,
                                  # 'Doublets'=cell_info[match(colnames(SCT_data_CD45neg),cell_info$cell),]$Doublets,
                                  col=cols, 
                                  show_annotation_name=FALSE,annotation_name_side = 'right')

        annotation_list <- list(col_ha,col_ha_CD45pos,col_ha_CD45neg)

    }

    # make heatmap row annotations
    make_heatmap_row_annotation_for_cells <- function(cell_info){

        # // prepare cell clusters annotations for heatmap
        print("prepare cell clusters annotations for heatmap")
        cell_clusters     <- unique(as.character(sort(as.numeric(as.vector(cell_info$seurat_clusters)))))

        cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                              "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                              "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                              "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")

        cluster_cols             = cluster_cols[1:length(cell_clusters)]
        names(cluster_cols)      = cell_clusters

        cols=list('Seurat Cluster' = cluster_cols)


        row_ha = rowAnnotation(
                                  'Seurat Cluster'=cell_info[match(colnames(SCT_data),cell_info$cell),]$seurat_clusters,
                                  col=cols, 
                                  show_annotation_name=FALSE)

        row_ha

    }

    all(cell_info$cell==colnames(RNA_data))

    mat_cor <- cor(RNA_data)

    #saveRDS(mat_cor,file=paste0(outDir,"/","heatmap_RNAdata_correlation.rds"))

    # make heatmap annotations
    annotation_list <- make_heatmap_annotations(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg)
    col_ha          <- annotation_list[[1]]

    row_ha <- make_heatmap_row_annotation_for_cells(cell_info)

    p <- Heatmap(mat_cor,name="Correlation Coefficient",column_title=paste0(ncol(mat_cor)," cells"),left_annotation = row_ha, top_annotation = col_ha,cluster_rows = F, cluster_columns = F,show_row_dend=F,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    
    pdf(paste0(outDir,"/","heatmap_RNAdata_correlation.pdf"),14,8); draw(p); dev.off()

    print("Heatmap plotted")

}


samples=c()
for(sample in samples){

  print(sample)

  file_rdata=paste0(sample,"/seurat_cellRangerData_emptyDrops_scDblFinder_cells_results/InSilicoCD45sorting_DEGs/data.Rdata")
  outDir=paste0(sample,"/seurat_cellRangerData_emptyDrops_scDblFinder_cells_results/InSilicoCD45sorting_DEGs/")

  make_cell_cell_corr_plot(file_rdata,outDir)

}
