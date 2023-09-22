# Goal: In silico CD45 sorting
# Methods: 1) define CD45+/- identify highly variable genes; 

# rm(list=ls())

library(dplyr)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)

make_heatmap_atac <- function(obj,outDir="outdir",show_row_names=FALSE,cluster_rows=T,peak_info=NULL,peaks=NULL){

    # // prepare data

    prepare_data <- function(obj,outDir="InSilicoCD45sorting",peaks=NULL){

        dir.create(outDir,recursive = T)
        
        if(class(obj)!='Seurat' & is.character(obj)){
            if(file.exists(obj))
                obj <- readRDS(obj)
        }

        cell_info <- obj@meta.data %>% 
                    mutate(cell=rownames(obj@meta.data))

        atac_data  = as.matrix(obj$ATAC@data)

        print("extract data")
        print(peaks[!peaks %in% rownames(atac_data)])
        peaks <- peaks[peaks %in% rownames(atac_data)]
        atac_data  <- t(t(atac_data)[,peaks,drop=F]); print(dim(atac_data))
        
        # out

        write.table(cell_info,paste0(outDir,"/metadata_cells.txt"),quote=F,sep="\t",row.names=F)

        out <- list(cell_info=cell_info,atac_data=atac_data)

    }

    print(outDir)
    
    dir.create(outDir,recursive = T)

    out = prepare_data(obj=obj,outDir=outDir,peaks=peaks)

    cell_info <- out$cell_info
    atac_data  <- out$atac_data

    saveRDS(out,file=paste0(outDir,"/","data.rds"))

    
    
    ########################################################################################
    # make heatmap 
    ########################################################################################
    print("make heatmap")

    # make heatmap annotations
    make_heatmap_annotations <- function(cell_info,atac_data){

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
                                  'Cluster'=cell_info[match(colnames(atac_data),cell_info$cell),]$seurat_clusters,
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')


        annotation_list <- list(col_ha)

        if(!is.null(peak_info)){
            row_ha <- rowAnnotation(
                'Cluster'=peak_info[match(rownames(atac_data),peak_info$gene),]$cluster,
                col=cols
            )
            annotation_list <- list(col_ha,row_ha)
        }

    }


    # make heatmap annotations
    annotation_list <- make_heatmap_annotations(cell_info,atac_data)
    col_ha          <- annotation_list[[1]]
    if(length(annotation_list)==2){
        row_ha <- annotation_list[[2]]
    }

    # replace peak in atac_data to alt_id with gene name
    if('gene_name' %in% colnames(peak_info)){
        peak_info <- peak_info %>% mutate(alt_id=paste0(gene_name,' (',gene,')'))
        rownames(atac_data) <- peak_info[match(rownames(atac_data),peak_info$gene),]$alt_id
    }

    # use raw RNA data
    #p <- Heatmap(atac_data,name="Log Counts",column_title=paste0(ncol(atac_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = F,show_row_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),width=unit(8,"inch"))
    
    #pdf(paste0(outDir,"/","heatmap_atac_logcounts.pdf"),20,10); draw(p); dev.off()

    # use raw RNA data, clustering cells
    p <- Heatmap(atac_data,name="Log Counts",column_title=paste0(ncol(atac_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = T,column_split = cell_info[match(colnames(atac_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),width=unit(8,"inch"))
    
    pdf(paste0(outDir,"/","heatmap_atac_logcounts.clustered_cells.pdf"),20,10); draw(p); dev.off()

    if(length(annotation_list)==2){

        p <- Heatmap(atac_data,name="Log Counts",column_title=paste0(ncol(atac_data)," cells"),top_annotation = col_ha,left_annotation = row_ha, cluster_rows = cluster_rows, cluster_columns = T,row_split=peak_info[match(rownames(out$atac_data),peak_info$gene),]$cluster,column_split = cell_info[match(colnames(atac_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),row_title_side='right',width=unit(8,"inch"))
    
        pdf(paste0(outDir,"/","heatmap_atac_logcounts.clustered_cells.splitGenes.pdf"),20,10); draw(p); dev.off()

        p2 <- Heatmap(atac_data,name="Log Counts",column_title=paste0(ncol(atac_data)," cells"),top_annotation = col_ha,left_annotation = row_ha, cluster_rows = cluster_rows, cluster_columns = T,row_split=peak_info[match(rownames(out$atac_data),peak_info$gene),]$cluster,column_split = cell_info[match(colnames(atac_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = TRUE ,show_column_names = F,row_names_side='left',row_names_gp = gpar(fontsize = 8),use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10),row_title_side='right',width=unit(8,"inch"))

        page_height=20
        if(nrow(atac_data)<100)
            page_height=15
        pdf(paste0(outDir,"/","heatmap_atac_logcounts.clustered_cells.splitGenes.showGeneName.pdf"),20,page_height); draw(p2); dev.off()

        save(atac_data,row_ha,col_ha,cluster_rows,peak_info,cell_info,show_row_names,p,p2,file=paste0(outDir,"/","heatmap_RNA_logcounts.clustered_cells.splitGenes.Rdata"))
    }

}
