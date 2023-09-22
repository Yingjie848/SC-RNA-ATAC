# explore cell origin of tumors
# Question: should we exclude genes expressed in few cells before calculating gene signature score?

library(GSEABase)
library(GSVA)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)



# Do ssGSEA analysis
ssgsea_analysis <- function(obj,outDir,gene.sets,gene.sets.name="geneset"){

    print(gene.sets.name)

    library(reshape2)
    source("seurat_workflow/code/annotate_cells.R")

    dir.create(paste0(outDir,'/',gene.sets.name),recursive = T)

    obj <- annotate_cells(obj)

    # calculate gsea score
    gsva_output <- data.frame()
    #if(!file.exists(paste0(outDir,'/',gene.sets.name,'/gsva_output.txt'))){
        gsva_output <- gsva(as.matrix(obj$RNA@data),gene.sets,method="ssgsea",ssgsea.norm=TRUE)
        print("GSVA done")
    #}else{
    #    gsva_output <- read.table(paste0(outDir,'/',gene.sets.name,'/gsva_output.txt'),header=T,check.names=F)
    #    print("GSVA results loaded")
    #}
    print(gsva_output[1,1:5])

    # make heatmap

    cell_info <- obj@meta.data %>% mutate(cell=rownames(obj@meta.data)) %>% arrange(seurat_clusters)

    # set colors for cell clusters
    cell_clusters     <- sort(unique(cell_info$seurat_clusters))
    cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
    cluster_cols             = cluster_cols[1:length(cell_clusters)]
    names(cluster_cols)      = cell_clusters

    library(circlize)
    CD45_cols            = colorRamp2( c(min(log10(cell_info$CD45+1)), max(log10(cell_info$CD45+1))), c("#f5f2f2", "red"))
    panck_cols           = colorRamp2( c(0, max(log10(cell_info$PanCK+1))), c("white", "black"))
    CellTotalCounts_cols = colorRamp2( c(0, max(log10(cell_info$nCount_RNA+1))), c("#f5f2f2", "#ff00a6"))
    CellTotalGenes_cols  = colorRamp2( c(0, max(log10(cell_info$nFeature_RNA+1))), c("#f5f2f2", "#01c504"))
    
    cols=list('Seurat Cluster' = cluster_cols,
                  'CD45 counts (log10)'=CD45_cols,
                  'Pan-CK (log10)'=panck_cols,
                  'Cell Total Counts (log10)'=CellTotalCounts_cols,
                  'Cell Total Genes (log10)'=CellTotalGenes_cols)

    col_ha = HeatmapAnnotation(
                                  'Seurat Cluster'=cell_info$seurat_clusters,
                                  'CD45 counts (log10)'=log10(cell_info$CD45+1),
                                  'Pan-CK (log10)'=log10(cell_info$PanCK),
                                  'Cell Total Counts (log10)'=log10(cell_info$nCount_RNA),
                                  'Cell Total Genes (log10)'=log10(cell_info$nFeature_RNA),
                                  col=cols, 
                                  show_annotation_name=TRUE,annotation_name_side = 'left')

    # ssGSEA
    p1 <- Heatmap(gsva_output[,cell_info$cell,drop=F],name="ssGSEA",column_title=paste0(ncol(gsva_output)," cells"),top_annotation = col_ha,cluster_rows = T, cluster_columns = T,column_split=cell_info$seurat_clusters,show_row_dend=F,show_column_dend=F,show_row_names = T ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(0.5, "mm"),row_title_gp = gpar(fontsize=6),row_names_gp = gpar(fontsize = 8),width=unit(7,units="inch"))
    # Z-score
    zscore <- t(scale(t(gsva_output[,cell_info$cell,drop=F])))
    p2 <- Heatmap(zscore,name="ssGSEA",column_title=paste0(ncol(gsva_output)," cells"),top_annotation = col_ha,cluster_rows = T, cluster_columns = T,column_split=cell_info$seurat_clusters,show_column_dend=F,show_row_dend=F,show_row_names = T ,show_column_names = F,row_names_side='left',use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=6),row_names_gp = gpar(fontsize = 8),width=unit(7,units="inch"))

    # boxplot for each geneset
    
    df_gsva <- melt(gsva_output)
    colnames(df_gsva) <- c('geneset','cell','ssGSEA')
    df_gsva <- merge(df_gsva,cell_info[,c('cell','seurat_clusters')])

    plots <- lapply(rownames(gsva_output),function(gs){
        avg_score <- mean(df_gsva[df_gsva$geneset==gs,]$ssGSEA,na.rm=T)
        ggviolin(df_gsva %>% filter(geneset==gs),x='seurat_clusters',y='ssGSEA',add='boxplot',color='seurat_clusters',title=gs,xlab="") + theme(legend.position = 'none') + geom_hline(yintercept = avg_score,col='grey')
    })
    p3 <- ggarrange(plotlist=plots,ncol=1,nrow=3)

    # out
    pdf(paste0(outDir,'/',gene.sets.name,'/ssgsea_heatmap.pdf'),15,10);print(p1);dev.off()
    pdf(paste0(outDir,'/',gene.sets.name,'/ssgsea_zscore_heatmap.pdf'),15,10);print(p2);dev.off()
    pdf(paste0(outDir,'/',gene.sets.name,'/ssgsea_violoinPlot.pdf'),10,7);print(p3);dev.off()

    write.table(gsva_output,paste0(outDir,'/',gene.sets.name,'/gsva_output.txt'),,quote=F,sep="\t",row.names=T)
    write.table(cell_info,paste0(outDir,'/',gene.sets.name,'/metaData.txt'),,quote=F,sep="\t",row.names=F)

}

