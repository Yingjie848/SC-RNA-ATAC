# explore cell origin of tumors
# Question: should we exclude genes expressed in few cells before calculating gene signature score?

library(GSEABase)
library(GSVA)
library(vcd)
library(Seurat)
library(UCell)
library(ggpubr)

source("seurat_workflow/code/calc_sigscore.R")
source("seurat_workflow/code/make_heatmap.R")



# cluster cells according to given genes
clustering_cells_by_given_genes <- function(obj,genes,resolution=0.2,outFigure){

    # set default colors for ggplot, which are used for cell clusters
    cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                    "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                    "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                    "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
    options(ggplot2.discrete.colour = cluster_cols)

    genes <- genes[genes %in% rownames(obj$RNA@data)]
    #obj.2 <- NormalizeData(obj,assay='RNA')
    #obj.2 <- ScaleData(obj.2,assay='RNA')
    obj.2 <- RunPCA(obj,assay='SCT',features=genes,npcs = 30, verbose = FALSE)
    obj.2 <- RunUMAP(obj.2,assay='SCT',reduction = "pca", dims = 1:30, verbose = FALSE)
    obj.2 <- FindClusters(obj.2,resolution = resolution, verbose = FALSE)

    p <- DimPlot(obj.2, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("")
    pdf(outFigure,8,7); print(p);dev.off()
    obj.2

}

# explore cell origin using methods:
# 1) signature score: calculated from average log2CPM for signature genes
# 2) ssGSEA for signature genes
# 3) UCell score for signature genes: UCell scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity
make_ternaryPlot <- function(obj,outDir,gene.sets,gene.sets.name="geneset"){

    # set default colors for ggplot, which are used for cell clusters
    cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
    cell_clusters     <- sort(unique(obj@meta.data$seurat_clusters))
    cluster_cols             = cluster_cols[1:length(cell_clusters)]
    names(cluster_cols)      = cell_clusters

    dir.create(outDir,recursive = T)

    #genesets <- getGmt(gmt.file)

    #MaSC_genes <- geneIds(genesets[['LIM_MAMMARY_STEM_CELL_UP']])
    #LP_genes   <- geneIds(genesets[['LIM_MAMMARY_LUMINAL_PROGENITOR_UP']])
    #ML_genes   <- geneIds(genesets[['LIM_MAMMARY_LUMINAL_MATURE_UP']])
    #ene.sets <- list(MaSC=MaSC_genes,LP=LP_genes,ML=ML_genes)

    # calculate signature score using average of CPM/100
    MaSC_score <- calc_sigscore(obj$RNA@counts,gene.sets$MaSC)
    LP_score   <- calc_sigscore(obj$RNA@counts,gene.sets$LP)
    ML_score   <- calc_sigscore(obj$RNA@counts,gene.sets$ML)

    # calculate gsea score
    gsva_output <- gsva(as.matrix(obj$RNA@counts),gene.sets,method="ssgsea",ssgsea.norm=TRUE)

    # calculate UCell score
    ucell_score <- ScoreSignatures_UCell(obj$RNA@counts, features=gene.sets,maxRank=2000)

    # combine scores
    df <- data.frame(cell=names(MaSC_score),
                        MaSC_score=MaSC_score,
                        LP_score=LP_score,
                        ML_score=ML_score,
                        MaSC_gsea=t(gsva_output)[,'MaSC']+abs(min(gsva_output)),
                        LP_gsea=t(gsva_output)[,'LP']+abs(min(gsva_output)),
                        ML_gsea=t(gsva_output)[,'ML']+abs(min(gsva_output)),
                        MaSC_ucell_score=ucell_score[,1],
                        LP_ucell_score=ucell_score[,2],
                        ML_ucell_score=ucell_score[,3]
                        )
    
    
    obj@meta.data <- obj@meta.data %>% mutate(cell=rownames(obj@meta.data))
    df <- merge(obj@meta.data,df)

    write.table(df,paste0(outDir,"/",gene.sets.name,".signature_score.txt"),quote=F,sep="\t",row.names=F)

    # ternary plot: signature score
    pdf(paste0(outDir,"/",gene.sets.name,".signature_score_ternaryplot.pdf"))
    df_tmp <- as.matrix(df[,c('MaSC_score','LP_score','ML_score')])
    rownames(df_tmp) <- df$cell
    ternaryplot(df_tmp, scale = 100, dimnames_position='corner',
        pch = 19, cex = 0.3, col = cluster_cols[match(df$seurat_clusters,names(cluster_cols))],
        dimnames = c("MaSC", "LP", "ML"),
        labels = c("outside"),
        prop_size=FALSE,
        main = "")
    dev.off()

    pdf(paste0(outDir,"/",gene.sets.name,".gsea_ternaryplot.pdf"))
    df_tmp <- as.matrix(df[,c('MaSC_gsea','LP_gsea','ML_gsea')])
    rownames(df_tmp) <- df$cell
    ternaryplot(df_tmp[rowSums(df_tmp)!=0,], scale = 100, dimnames_position='corner',
        pch = 19, cex = 0.3, col = cluster_cols[match(df$seurat_clusters,names(cluster_cols))],
        dimnames = c("MaSC", "LP", "ML"),
        labels = c("outside"),
        prop_size=FALSE,
        main = "")
    dev.off()


    pdf(paste0(outDir,"/",gene.sets.name,".ucell_score_ternaryplot.pdf"))
    df_tmp <- as.matrix(df[,c('MaSC_ucell_score','LP_ucell_score','ML_ucell_score')])
    rownames(df_tmp) <- df$cell
    ternaryplot(df_tmp[rowSums(df_tmp)!=0,], scale = 100, dimnames_position='corner',
        pch = 19, cex = 0.3, col = cluster_cols[match(df$seurat_clusters,names(cluster_cols))],
        dimnames = c("MaSC", "LP", "ML"),
        labels = c("outside"),
        prop_size=FALSE,
        main = "")
    dev.off()


    # make boxplot

    p1 <- df %>% ggviolin(x='seurat_clusters',y='MaSC_score',color='seurat_clusters',title='MaSC',xlab='Clusters',add='boxplot') + scale_color_manual(values=cluster_cols) + theme(legend.position = 'none') 
    p2 <- df %>% ggviolin(x='seurat_clusters',y='LP_score',color='seurat_clusters',title='LP',xlab='Clusters',add='boxplot') + scale_color_manual(values=cluster_cols) + theme(legend.position = 'none')
    p3 <- df %>% ggviolin(x='seurat_clusters',y='ML_score',color='seurat_clusters',title='ML',xlab='Clusters',add='boxplot') + scale_color_manual(values=cluster_cols) + theme(legend.position = 'none')

    p <- ggarrange(p1,p2,p3,ncol=3)
    pdf(paste0(outDir,"/",gene.sets.name,".signature_score_boxplot.pdf"),10,5);print(p);dev.off()
}

explore_cell_origin <- function(obj,do_clustering=TRUE,resolution=0.2,outDir,gene.sets,gene.sets.name){

        print(outDir)

        dir.create(paste0(outDir,'/',gene.sets.name),recursive = T)

        genes      <- unique(unlist(gene.sets))

        # cluster cells
        if(do_clustering){
            obj.2 <- clustering_cells_by_given_genes(obj,genes,resolution = resolution,outFigure=paste0(outDir,'/',gene.sets.name,'/',gene.sets.name,'_DimPlot_umap.pdf'))
        }else{
            obj.2 <- obj
        }

        # make ternaryPlot
        make_ternaryPlot(obj.2,outDir=paste0(outDir,'/',gene.sets.name),gene.sets,gene.sets.name=gene.sets.name)

        # make heatmap
        if(length(gene.sets$MaSC) < 50){
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/MaSC_heatmap"),genes=gene.sets$MaSC,show_row_names = T,geneExXCells=20,normalize_to_logcounts=TRUE)
        }else{
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/MaSC_heatmap"),genes=gene.sets$MaSC,show_row_names = F,geneExXCells=20,normalize_to_logcounts=TRUE)
        }
        if(length(gene.sets$LP) < 50){
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/LP_heatmap"),  genes=gene.sets$LP,show_row_names = T,geneExXCells=20,normalize_to_logcounts=TRUE)
        }else{
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/LP_heatmap"),  genes=gene.sets$LP,show_row_names = F,geneExXCells=20,normalize_to_logcounts=TRUE)
        }
        if(length(gene.sets$ML) < 50){
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/ML_heatmap"),  genes=gene.sets$ML,show_row_names = T,geneExXCells=20,normalize_to_logcounts=TRUE)
        }else{
            make_heatmap(obj.2,outDir=paste0(outDir,"/",gene.sets.name,"/ML_heatmap"),  genes=gene.sets$ML,show_row_names = F,geneExXCells=20,normalize_to_logcounts=TRUE)
        }

}


