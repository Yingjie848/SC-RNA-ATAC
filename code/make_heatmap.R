# Goal: In silico CD45 sorting
# Methods: 1) define CD45+/- identify highly variable genes; 

# rm(list=ls())

library(dplyr)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)

make_heatmap <- function(obj,outDir="outdir",geneExXCells=1,CD45posMinReads=1,genes=NULL,show_row_names=F,normalize_to_logcounts=FALSE,cluster_rows=T,gene_info=NULL,page_height_when_show_gene_name=NULL){

    # // prepare data

    prepare_data <- function(obj,outDir="InSilicoCD45sorting",geneExXCells=1,CD45posMinReads=1,genes=NULL,normalize_to_logcounts=FALSE){


        calc_cpm <- function(counts){

            a = colSums(counts) # get sums from each column
            mat = as.matrix(counts)
            c = sweep(mat, 2, a, '/') # divide the matrix by the vector
            cpm = c * 1e6 # multiple 1 million to form cpm

        }


        dir.create(outDir,recursive = T)
        
        if(class(obj)!='Seurat' & is.character(obj)){
        if(file.exists(obj))
            obj <- readRDS(obj)
        }
        
        cat("Number of cells: ",nrow(obj@meta.data),"\n")
        cat("Number of genes: ",nrow(obj),"\n")
        cat("Number of genes in gene list: ",length(unique(genes)),"\n")
        cat("Number of genes in gene list and RNA object: ",sum(genes %in% rownames(obj$RNA@data)),"\n")

        cell_info <- obj@meta.data %>% 
                    mutate(cell=rownames(obj@meta.data))
        cell_info <- cell_info[,! colnames(cell_info) %in% c('CD45','PanCK')]

        RNA_data  = as.matrix(obj$RNA@counts)

        if(isTRUE(normalize_to_logcounts)){
            obj <- NormalizeData(obj,assay="RNA")
            logcounts <- as.matrix(obj$RNA@data)
        }

        # create Pan Cytokeratin (Pan-CK) counts 
        print("create Pan Cytokeratin (Pan-CK) counts")
        panck_genes = rownames(RNA_data)[grepl('^KRT',rownames(RNA_data))]
        panck_genes = panck_genes[!grepl('CAP',panck_genes)]

        RNA_data_panck <- RNA_data[panck_genes,]

        df_panck <- data.frame(cell=colnames(RNA_data_panck),PanCK=colSums(RNA_data_panck))
        cell_info <- merge(cell_info,df_panck,by='cell',all.x=T,sort=F)

        # divide CD45+/CD45- cells
        print("divide CD45+/CD45- cells")
        cat("cutoff:",CD45posMinReads,"\n")
        CD45_expr   = t(RNA_data)[,'PTPRC']
        CD45_counts = data.frame(cell=names(CD45_expr),CD45=CD45_expr) %>% 
                    arrange(desc(CD45)) %>% 
                    mutate(CD45_status=ifelse(CD45>=CD45posMinReads,'CD45+','CD45-'))

        # order cells by CD45 counts
        cell_info   = merge(cell_info,CD45_counts,by='cell',all.x=T,sort=F) %>% 
                    arrange(desc(CD45)) 
        print("Cells ordered by seurat_clusters,desc(CD45)")


        # filter genes
        RNA_data  = RNA_data[,cell_info$cell]
        nexpressed_cells = apply(RNA_data,1,function(x) sum(x>0))
        RNA_data  = RNA_data[nexpressed_cells >= geneExXCells,] # filter out genes expressed in less than 20 cells 
        cat("Number of genes in gene list and RNA object after filtering: ",sum(genes %in% rownames(RNA_data)),"\n")

        # extract data
        print("extract data")
        print(genes[!genes %in% rownames(RNA_data)])
        genes <- genes[genes %in% rownames(RNA_data)]
        RNA_data  <- t(t(RNA_data)[,genes,drop=F]); print(dim(RNA_data))




        # // order cells by cell clusters, clusters having more CD45+ cells will be put at the beginning

        # get frequency of CD45+ cells for each cluster
        cell_cluster_CD45_freq <- cell_info %>% 
                                group_by(seurat_clusters) %>% 
                                summarise(ncells=length(cell),CD45_frequency=sum(CD45>=CD45posMinReads)) %>% 
                                mutate(CD45_cell_pct=CD45_frequency/ncells*100) %>% 
                                arrange(desc(CD45_cell_pct)) # order cell clusters by CD45+ frequency
        cell_info.ordered_by_clusters_ranked_by_CD45 <- lapply(as.vector(cell_cluster_CD45_freq$seurat_clusters),
                                                                function(x){
                                                                tmp = cell_info[cell_info$seurat_clusters==x,] %>% arrange(desc(CD45))
                                                                })
        cell_info.ordered_by_clusters_ranked_by_CD45 <- do.call(rbind,cell_info.ordered_by_clusters_ranked_by_CD45) # reorder cells by cell clusters, clusters having more CD45+ cells will be put at the beginning
        cell_info <- cell_info.ordered_by_clusters_ranked_by_CD45


        # for RNA_data, they are ordered by cell clusters ranked by frequency of CD45+ cells 

        # now extract data
        RNA_data         = RNA_data[,as.vector(cell_info$cell)]

        if(isTRUE(normalize_to_logcounts)){
            print("Use data normalized to logcounts")
            RNA_data         = logcounts[rownames(RNA_data),colnames(RNA_data)]
        }
        
        # out

        write.table(cell_info,paste0(outDir,"/metadata_cells.txt"),quote=F,sep="\t",row.names=F)

        out <- list(cell_info=cell_info,RNA_data=RNA_data)

    }

    print(outDir)
    
    dir.create(outDir,recursive = T)

    out = prepare_data(obj=obj,outDir=outDir,geneExXCells = geneExXCells,CD45posMinReads = CD45posMinReads,genes=genes,normalize_to_logcounts=normalize_to_logcounts)

    cell_info <- out$cell_info
    RNA_data  <- out$RNA_data

    saveRDS(out,file=paste0(outDir,"/","output.rds"))

    
    ########################################################################################
    # make heatmap 
    ########################################################################################
    print("make heatmap")

    # make heatmap annotations
    make_heatmap_annotations <- function(cell_info,RNA_data){

        # // prepare cell clusters annotations for heatmap
        print("prepare cell clusters annotations for heatmap")
        cell_clusters     <- sort(unique(cell_info$seurat_clusters))

        cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                              "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                              "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                              "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
        cluster_cols             = cluster_cols[1:length(cell_clusters)]
        names(cluster_cols)      = cell_clusters

        library(circlize)
        CD45_cols            = colorRamp2( c(min(log10(cell_info$CD45+1)), max(log10(cell_info$CD45+1))), c("#f5f2f2", "red"))
        fibroblast_gsea_cols     = colorRamp2( c(min(cell_info$fibroblast_gsea), max(cell_info$fibroblast_gsea)), c("#f5f2f2", "red"))
        CellTotalCounts_cols = colorRamp2( c(0, max(log10(cell_info$nCount_RNA+1))), c("#f5f2f2", "#ff00a6"))
        CellTotalGenes_cols  = colorRamp2( c(0, max(log10(cell_info$nFeature_RNA+1))), c("#f5f2f2", "#01c504"))
        panck_cols           = colorRamp2( c(0, max(log10(cell_info$PanCK+1))), c("white", "black"))
        dbl_cols             = colorRamp2( c(0, max(cell_info$scDblFinder.score)), c("white", "purple"))
        scDbl_cols           = c('singlet'='white','doublet'='purple')
        hybrid_cols          = colorRamp2( c(0, max(cell_info$hybrid_score)), c("white", "darkgreen"))
        doublet_cols         = c('No'='white','Yes'='purple')

        cols=list('Cluster' = cluster_cols,
                  'CD45 counts (log10)'=CD45_cols,
                  'Fibroblast GSEA'=fibroblast_gsea_cols,
                  'Cell Total Counts (log10)'=CellTotalCounts_cols,
                  'Cell Total Genes (log10)'=CellTotalGenes_cols,
                  'Pan-CK (log10)'=panck_cols,
                  'Doublets (scDblFinder)'=scDbl_cols,
                  'Doublets score (scDblFinder)'=dbl_cols,
                  'Doublets score (scds)'=hybrid_cols,
                  'Doublets'=doublet_cols)


        col_ha = HeatmapAnnotation(
                                  'Cluster'=cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,
                                  'CD45 counts (log10)'=log10(cell_info$CD45+1),
                                  #'Fibroblast GSEA'=cell_info$fibroblast_gsea,
                                  'Pan-CK (log10)'=log10(cell_info$PanCK),
                                  'Cell Total Counts (log10)'=log10(cell_info$nCount_RNA),
                                  'Cell Total Genes (log10)'=log10(cell_info$nFeature_RNA),
                                  #'Doublets (scDblFinder)'=cell_info$scDblFinder.doublets,
                                  #'Doublets score (scDblFinder)'=cell_info$scDblFinder.score,
                                  #'Doublets score (scds)'=cell_info$hybrid_score,
                                  # 'Doublets'=cell_info$Doublets,
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')


        annotation_list <- list(col_ha)

        if(!is.null(gene_info)){
            row_ha <- rowAnnotation(
                'Cluster'=gene_info[match(rownames(RNA_data),gene_info$gene),]$cluster,
                col=cols
            )
            annotation_list <- list(col_ha,row_ha)
        }

    }


    # make heatmap annotations
    annotation_list <- make_heatmap_annotations(cell_info,RNA_data)
    col_ha          <- annotation_list[[1]]
    if(length(annotation_list)==2){
        row_ha <- annotation_list[[2]]
    }


    if(isTRUE(normalize_to_logcounts)){

        # use raw RNA data
        p <- Heatmap(RNA_data,name="Log Counts",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = F,show_row_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNA_logcounts.pdf"),15,10); draw(p); dev.off()

        # use raw RNA data, clustering cells
        p <- Heatmap(RNA_data,name="Log Counts",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = T,column_split = cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNA_logcounts.clustered_cells.pdf"),15,10); draw(p); dev.off()

        if(length(annotation_list)==2){

            p <- Heatmap(RNA_data,name="Log Counts",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,left_annotation = row_ha, cluster_rows = cluster_rows, cluster_columns = T,row_split=gene_info[match(rownames(RNA_data),gene_info$gene),]$cluster,column_split = cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=14),column_title_gp = gpar(fontsize=14))
        
            pdf(paste0(outDir,"/","heatmap_RNA_logcounts.clustered_cells.splitGenes.pdf"),15,10); draw(p); dev.off()

            p2 <- Heatmap(RNA_data,name="Log Counts",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,left_annotation = row_ha, cluster_rows = cluster_rows, cluster_columns = T,row_split=gene_info[match(rownames(RNA_data),gene_info$gene),]$cluster,column_split = cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = TRUE ,show_column_names = F,row_names_side='left',row_names_gp = gpar(fontsize = 10),use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=14),column_title_gp = gpar(fontsize=14))

            page_height=20
            if(nrow(RNA_data)<100)
                page_height=15
            if(!is.null(page_height_when_show_gene_name)){
                page_height=page_height_when_show_gene_name
            }
            pdf(paste0(outDir,"/","heatmap_RNA_logcounts.clustered_cells.splitGenes.showGeneName.pdf"),15,page_height); draw(p2); dev.off()

            save(RNA_data,row_ha,col_ha,cluster_rows,gene_info,cell_info,show_row_names,p,p2,file=paste0(outDir,"/","heatmap_RNA_logcounts.clustered_cells.splitGenes.Rdata"))
        }

        # use z-score
        scaled_RNA_data <- t(scale(t(RNA_data)))
        scaled_RNA_data <- scales::squish(scaled_RNA_data,range = c(-0.5,0.5))

        p <- Heatmap(scaled_RNA_data,name="Log Counts (Z-score)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = F,show_row_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNA_logcounts_zscore.pdf"),15,10); draw(p); dev.off()

    }else{
  
        # use raw RNA data
        p <- Heatmap(log10(RNA_data+1),name="RNA counts (log10)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = F,show_row_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNAcounts.pdf"),15,10); draw(p); dev.off()

        # use raw RNA data, clustering cells
        p <- Heatmap(log10(RNA_data+1),name="RNA counts (log10)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = T,column_split = cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNAcounts.clustered_cells.pdf"),15,10); draw(p); dev.off()

        # use z-score
        scaled_RNA_data <- t(scale(t(RNA_data)))
        scaled_RNA_data <- scales::squish(scaled_RNA_data,range = c(-0.5,0.5))

        p <- Heatmap(scaled_RNA_data,name="RNA counts (Z-score)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = cluster_rows, cluster_columns = F,show_row_dend=F,show_row_names = show_row_names ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
        
        pdf(paste0(outDir,"/","heatmap_RNAcounts_zscore.pdf"),15,10); draw(p); dev.off()

    }


    ########################################################################################
    # make cell-cell correlation heatmap
    ########################################################################################

    
    # make heatmap annotations
    make_heatmap_row_annotation_for_cells <- function(cell_info){

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


        row_ha = rowAnnotation(
                                  'Cluster'=cell_info[match(colnames(RNA_data),cell_info$cell),]$seurat_clusters,
                                  col=cols, 
                                  show_annotation_name=FALSE)

        row_ha

    }

    mat_cor <- cor(RNA_data)

    #saveRDS(mat_cor,file=paste0(outDir,"/","heatmap_RNAdata_correlation.rds"))

    row_ha <- make_heatmap_row_annotation_for_cells(cell_info)

    p <- Heatmap(mat_cor,name="Correlation Coefficient",column_title=paste0(ncol(mat_cor)," cells"),left_annotation = row_ha, top_annotation = col_ha,cluster_rows = F, cluster_columns = F,show_row_dend=F,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    pdf(paste0(outDir,"/","heatmap_RNAdata_correlation.pdf"),14,8); draw(p); dev.off()


}
