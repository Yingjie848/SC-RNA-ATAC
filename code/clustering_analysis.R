
# perform subclusering for cells
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)


clustering_analysis <- function(obj,outDir,npc=50,resolution = 1, stop_after_clustering=FALSE,stop_after_marker_genes = FALSE,regress_out_cellcycle=FALSE,clustering_assay='RNA',secondary_analysis=TRUE ){

        dir.create(outDir)

        #########################################################################################################
        # Add CD45 and PanCK counts
        print("Add CD45 and PanCK counts")
        source("seurat_workflow/code/annotate_cells.R")
        obj <- annotate_cells(obj)


        #########################################################################################################
        # scoring cell cycle
        # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
        # segregate this list into markers of G2/M phase and markers of S phase
        print("Scoring cell cycle")

        print("Set default assay to RNA")
        DefaultAssay(obj) <- "RNA"

        s.genes   <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes

        print("Normalize data")
        obj <- NormalizeData(obj,assay='RNA')
        print("Scale data")
        obj <- ScaleData(obj)
        print("Scoring cell cycle")
        obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

        # scale data and regress out variables
        if(regress_out_cellcycle){
          print("Regress out cell cycle from scaled data in RNA assay")
          obj <- ScaleData(obj,vars.to.regress=c("S.Score", "G2M.Score"))
        }

        

        #########################################################################################################
        # do clustering

        if(!regress_out_cellcycle){
            cat("Do clustering without regressing out cell cycle\n")
            obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE) %>%
                      RunPCA(npcs = npc, verbose = FALSE) %>%
                      RunUMAP(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindNeighbors(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindClusters(resolution = resolution, verbose = FALSE)
            obj <- obj %>% RunTSNE(reduction = "pca", dims = 1:npc, verbose = FALSE)
        }else{
            cat("Do clustering with regressing out cell cycle\n")
            obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE,vars.to.regress=c("S.Score", "G2M.Score")) %>%
                      RunPCA(npcs = npc, verbose = FALSE) %>%
                      RunUMAP(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindNeighbors(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindClusters(resolution = resolution, verbose = FALSE)  
            obj <- ScaleData(obj, assay='RNA', vars.to.regress = c("S.Score", "G2M.Score"))
            obj <- obj %>% RunTSNE(reduction = "pca", dims = 1:npc, verbose = FALSE)
        }

        if('ATAC' %in% clustering_assay){
            # ATAC analysis
            # We exclude the first dimension as this is typically correlated with sequencing depth
            DefaultAssay(obj) <- "ATAC"
            obj <- RunTFIDF(obj) %>% 
                    FindTopFeatures(min.cutoff = 'q0') %>%
                    RunSVD() %>%
                    RunUMAP(reduction = 'lsi', dims = 2:npc, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

            # We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
            # We use this graph for UMAP visualization and clustering
            obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:npc, 2:npc)) %>%
                    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
                    FindClusters(graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
            obj <- RunTSNE(obj, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_")

        }

        cat("Number of clusters:",length(unique(obj@meta.data$seurat_clusters)))
        print(table(obj@meta.data$seurat_clusters))

        
        #########################################################################################################
        # save obj
        DefaultAssay(obj) <- 'RNA'
        saveRDS(obj,file=paste0(outDir,"/seurat_rna.rds"))


        #########################################################################################################
        if(stop_after_clustering){
          return("Stopped after clustering")
        }

        #########################################################################################################
        # secondary analysis
        if(secondary_analysis){
            secondary_analysis(obj,outDir,stop_after_marker_genes = stop_after_marker_genes, clustering_assay=clustering_assay )
        }

        print("Done")


}


subclustering <- function(obj,outDir='subclustered',cluster=NULL,subcluster.name='subcluster',resolution=0.5,stop_after_marker_genes = FALSE,clustering_assay="RNA", secondary_analysis=TRUE,evaluate_stemness=TRUE){

    dir.create(outDir,recursive = T)

    print("Normalize data")
    obj <- NormalizeData(obj,assay='RNA')
    print("Scale data")
    obj <- ScaleData(obj)

    if('ATAC' %in% clustering_assay){
        obj <- FindSubCluster(obj,cluster=cluster,subcluster.name=subcluster.name,graph.name='wsnn',algorithm=3,resolution=resolution)
    }else{
        obj <- FindSubCluster(obj,cluster=cluster,subcluster.name=subcluster.name,graph.name='SCT_snn',resolution=resolution)
    }

    obj@meta.data$seurat_clusters.1 <- obj@meta.data$seurat_clusters
    obj@meta.data$seurat_clusters <- factor(obj@meta.data[,subcluster.name],
                                            levels=sort(unique(obj@meta.data[,subcluster.name])))
    Idents(obj) <- obj@meta.data$seurat_clusters

    DefaultAssay(obj) <- 'RNA'
    saveRDS(obj,file=paste0(outDir,"/seurat_rna.rds"))


    #########################################################################################################
    # secondary analysis
    if(secondary_analysis){
        secondary_analysis(obj,outDir,stop_after_marker_genes = stop_after_marker_genes, clustering_assay=clustering_assay,evaluate_stemness=evaluate_stemness )
    }

    obj

}



secondary_analysis <- function(obj,outDir,stop_after_marker_genes = FALSE, stop_after_DEA=FALSE,clustering_assay='RNA', evaluate_stemness=TRUE ){

        dir.create(outDir)

        #########################################################################################################
        # QC
        cat("Create QC plots\n")

        if('ATAC' %in% clustering_assay){
            p1 <- VlnPlot(obj, features = c("nCount_RNA","nCount_ATAC","nFeature_RNA","nFeature_ATAC","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 4,log = TRUE, pt.size = 0.5) + NoLegend()
            p2 <- VlnPlot(obj, features = c("nCount_RNA","nCount_ATAC","nFeature_RNA","nFeature_ATAC","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 4,log = FALSE, pt.size = 0.5) + NoLegend()

            pdf(paste0(outDir,"/qc_metrics.with_clustering.log.pdf"),20,7);print(p1);dev.off()
            pdf(paste0(outDir,"/qc_metrics.with_clustering.raw.pdf"),20,7);print(p2);dev.off()

        }else{
            p1 <- VlnPlot(obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 3,log = TRUE, pt.size = 0.5) + NoLegend()
            p2 <- VlnPlot(obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 3,log = FALSE, pt.size = 0.5) + NoLegend()

            pdf(paste0(outDir,"/qc_metrics.with_clustering.log.pdf"),15,7);print(p1);dev.off()
            pdf(paste0(outDir,"/qc_metrics.with_clustering.raw.pdf"),15,7);print(p2);dev.off()
        }

        #########################################################################################################
        # visualize clustering based on gene expression 
        cat("Create dimentional reduction plots\n")

        p <- DimPlot(obj, reduction = "pca", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("PCA (RNA)")

        pdf(paste0(outDir,"/DimPlot_pca_rna.pdf"),5,5);print(p);dev.off()

        p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("UMAP (RNA)")

        pdf(paste0(outDir,"/DimPlot_umap_rna.pdf"),5,5);print(p);dev.off()

        p <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("tSNE (RNA)")

        pdf(paste0(outDir,"/DimPlot_tsne_rna.pdf"),5,5);print(p);dev.off()

        if('ATAC' %in% clustering_assay){

            p <- DimPlot(obj, reduction = "lsi", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("lsi (ATAC)")

            pdf(paste0(outDir,"/DimPlot_lsi_atac.pdf"),5,5);print(p);dev.off()

            p <- DimPlot(obj, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("UMAP (ATAC)")

            pdf(paste0(outDir,"/DimPlot_umap_atac.pdf"),5,5);print(p);dev.off()

            p <- DimPlot(obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("UMAP (RNA + ATAC)")

            pdf(paste0(outDir,"/DimPlot_umap_wnn.pdf"),5,5);print(p);dev.off()

        }

        # meta.data
        metaData = obj@meta.data %>% mutate(cell_id=rownames(obj@meta.data))
        write.table(metaData,paste0(outDir,"/metaData.txt"),quote=F,sep="\t",row.names=F)

        # count cells in cluster
        cells_per_cluster = metaData %>% group_by(orig.ident,seurat_clusters) %>% summarise(cells=length(cell_id))
        write.table(cells_per_cluster,paste0(outDir,"/metaData.cells_per_cluster.txt"),quote=F,sep="\t",row.names=F)


        #########################################################################################################
        # visualize marker genes
        cat("Visualize marker genes\n")

        source("seurat_workflow/code/make_bubblePlot.R")

        make_bubblePlot(obj,paste0(outDir,'/marker_gene_plots'))


        #########################################################################################################
        # make FeaturePlot in UMAP

        print("make FeaturePlot ")
        DefaultAssay(obj) <- 'RNA'

        if('ATAC' %in% clustering_assay){
          reduction_method='wnn.umap'
        }else{
          reduction_method='umap'
        }

        source("seurat_workflow/code/make_featurePlot.R")

        make_featurePlots(obj,outDir,reduction_method=reduction_method,order=FALSE)

        make_featurePlots_for_chondroid_markers(obj,paste0(outDir,"/FeaturePlot"),reduction_method=reduction_method,order=FALSE)
        make_featurePlots_for_chondroid_markers(obj,paste0(outDir,"/FeaturePlot_Ordered"),reduction_method=reduction_method,order=TRUE)

        make_featurePlots_for_myoepithelial_markers(obj,paste0(outDir,"/FeaturePlot"),reduction_method='umap',order=FALSE)
        make_featurePlots_for_myoepithelial_markers(obj,paste0(outDir,"/FeaturePlot_Ordered"),reduction_method='umap',order=TRUE)


        #########################################################################################################
        if(stop_after_marker_genes){
          return("Stopped after feature plot")
        }


        #########################################################################################################
        # // Finding differentially expressed features (cluster biomarkers)

        print("Finding differentially expressed features")

        DefaultAssay(obj) <- "SCT"

        # find markers for every cluster compared to all remaining cells
        dir.create(paste0(outDir,"/DE_analysis"))

        # identify markers genes both positive and negative, and take top50 genes for each cluster
        obj <- PrepSCTFindMarkers(obj)
        markers <- FindAllMarkers(obj, assay='SCT',test.use="wilcox",only.pos = FALSE, min.pct = 0.25, logfc.threshold = log2(1.2))

        markers_top50 <-markers %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(desc(abs(avg_log2FC))) %>%
          dplyr::filter(p_val_adj < 0.001) %>%
          dplyr::slice_head(n = 50)

        markers_top10 <- markers %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(desc(abs(avg_log2FC))) %>%
          dplyr::filter(p_val_adj < 0.001) %>%
          dplyr::slice_head(n = 10)

        write.table(markers,paste0(outDir,"/DE_analysis/DE_genes.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_top50,paste0(outDir,"/DE_analysis/DE_genes_top50.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_top10,paste0(outDir,"/DE_analysis/DE_genes_top10.txt"),quote=F,sep="\t",row.names=F)

        # markers, only the positive genes
        markers_pos <- markers %>% dplyr::filter(avg_log2FC>0 & p_val_adj < 0.001)
        markers_neg <- markers %>% dplyr::filter(avg_log2FC<0 & p_val_adj < 0.001)

        markers_pos_top50 <-markers_pos %>%
          dplyr::group_by(cluster) %>%
          dplyr::top_n(n = 50, wt = avg_log2FC)

        markers_pos_top10 <-markers_pos %>%
          dplyr::group_by(cluster) %>%
          dplyr::top_n(n = 10, wt = avg_log2FC)

        markers_neg_top10 <-markers_neg %>%
          dplyr::group_by(cluster) %>%
          dplyr::filter(avg_log2FC<0) %>%
          dplyr::arrange(avg_log2FC) %>%
          dplyr::slice_head(n = 10)
          

        write.table(markers_pos,paste0(outDir,"/DE_analysis/posDE_genes.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_neg,paste0(outDir,"/DE_analysis/negDE_genes.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_pos_top50,paste0(outDir,"/DE_analysis/posDE_genes_top50.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_pos_top10,paste0(outDir,"/DE_analysis/posDE_genes_top10.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_neg_top10,paste0(outDir,"/DE_analysis/negDE_genes_top10.txt"),quote=F,sep="\t",row.names=F)


        #########################################################################################################
        # // annotate up-regulated genes for each cluster

        cat("Annotate up-regulated genes for each cluster\n")

        source("seurat_workflow/code/pathway.R")

        annotate_posDEG_to_pathway(markers_pos,outDir,topn=5)
        annotate_posDEG_to_pathway(markers_pos,outDir,topn=10)
        

        #########################################################################################################
        # make heatmap

        print("Make heatmap for top 50 up-regulated DEGs")

        source("seurat_workflow/code/make_heatmap.R")

        # top 50 up-regulated genes

        gene_info <- markers_pos_top50 %>% group_by(gene) %>% arrange(p_val_adj) %>% summarise(cluster=cluster[1])

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top50_posDEGs_logcounts"),
                 genes=unique(markers_pos_top50$gene),
                 show_row_names = F,
                 normalize_to_logcounts = TRUE,
                 gene_info=gene_info)

        # top 10 up-regulated genes

        gene_info <- markers_pos_top10 %>% group_by(gene) %>% arrange(p_val_adj) %>% summarise(cluster=cluster[1])

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top10_posDEGs_logcounts"),
                 genes=unique(markers_pos_top10$gene),
                 show_row_names = F,
                 normalize_to_logcounts = TRUE,
                 gene_info=gene_info,page_height_when_show_gene_name = 15)


        # top 10 down-regulated genes

        gene_info <- markers_neg_top10 %>% group_by(gene) %>% arrange(p_val_adj) %>% summarise(cluster=cluster[1])

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top10_negDEGs_logcounts"),
                 genes=unique(markers_neg_top10$gene),
                 show_row_names = F,
                 normalize_to_logcounts = TRUE,
                 gene_info=gene_info,page_height_when_show_gene_name = 15)

        #make_heatmap(obj=obj,
        #         outDir=paste0(outDir,"/heatmap_top50_posDEGs_counts"),
        #         genes=unique(markers_pos_top50$gene),
        #         show_row_names = F,
        #         normalize_to_logcounts = FALSE)

        # for EMT markers

        EMT_markers <- c(# mesenchymal markers
                        'VIM', # vimentin
                        'FN1', # fibronectin
                        'CDH2', # N-cadherin
                        'FOXC2',
                        'SNAI1',
                        'SNAI2',
                        'SOX10',
                        'TWIST1',
                        'ITGB6', #integrin αvβ6
                        'MMP2',
                        'MMP3',
                        'MMP9',
                        'GCS',
                        # epithelial markers
                        'CDH1', # E-cadherin
                        'DSP', # desmoplakin
                        'OCLN' # occludin
                        )

        make_heatmap(obj=obj,
                      outDir=paste0(outDir,"/heatmap_EMT_markers"),
                      genes=EMT_markers,
                      show_row_names = T,
                      cluster_rows = F,
                      normalize_to_logcounts = TRUE)


        #########################################################################################################
        if(stop_after_DEA){
          return("Stopped after differential expression analysis")
        }


        ##########################################################################################################################################################################
        # GSEA analysis according to fold change
        source("seurat_workflow/code/gsea.R")

        # MSigDB
        gmt.file       <- "seurat_workflow/data/msigdb_bp_kegg_reactome_h.symbol.gmt"
        gene.sets      <- read.gmt(gmt.file)
        gene.sets.name <- 'MSigDB'

        do_gsea(obj,paste0(outDir,'/GSEA'),gene.sets,gene.sets.name,plot_nterms = 50)

        # HALLMARK
        gmt.file       <- "seurat_workflow/data/h.all.v7.2.symbols.gmt"
        gene.sets      <- read.gmt(gmt.file)
        gene.sets.name      <- 'HALLMARK'

        do_gsea(obj,paste0(outDir,'/GSEA'),gene.sets,gene.sets.name,plot_nterms = 50)

        # cell origin
        gmt.file       <- "seurat_workflow/data/marker_genes/cell_origin/cell_origin_Pommier_LIM.gmt"
        gene.sets      <- read.gmt(gmt.file)
        gene.sets.name      <- 'Cell_Origin'

        do_gsea(obj,paste0(outDir,'/GSEA'),gene.sets,gene.sets.name,plot_nterms = 50)

        # Stemness signature
        gmt.file       <- "seurat_workflow/data/marker_genes/cell_origin/my_stemness_signatures.gmt"
        gene.sets      <- read.gmt(gmt.file)
        gene.sets.name <- 'Stemness'

        do_gsea(obj,paste0(outDir,'/GSEA'),gene.sets,gene.sets.name,plot_nterms = 50)



        ##########################################################################################################################################################################
        # cell origin analysis

        print("cell origin analysis")

        source("seurat_workflow/code/explore_cell_origin.R")

        # Pommier gene sets
        gmt.file   <- "seurat_workflow/data/marker_genes/cell_origin/Pommier.genesymbol.gmt"
        genesets   <- getGmt(gmt.file)
        MaSC_genes <- geneIds(genesets[['Pommier_MaSC1-2-3_UP']])
        LP_genes   <- geneIds(genesets[['Pommier_LP_UP']])
        ML_genes   <- geneIds(genesets[['Pommier_mL_UP']])
        gene.sets  <- list(MaSC=MaSC_genes,LP=LP_genes,ML=ML_genes)
        gene.sets.name <- "Pommier"

        explore_cell_origin(obj,do_clustering=FALSE, resolution = NULL,outDir=paste0(outDir,"/cell_origin_notClustered"),gene.sets,gene.sets.name)

        # explore_cell_origin(obj,do_clustering=TRUE, resolution = 0.2,outDir=paste0(outDir,"/cell_origin_clustered"),gene.sets,gene.sets.name)


        ##########################################################################################################################################################################
        # pseudotime

        print("Peudotime")

        source("seurat_workflow/code/destiny.R")

        destiny_pca(obj,outDir=paste0(outDir,'/destiny_pca'))
        

        ##########################################################################################################################################################################
        # ssGSEA analysis

        print("ssGSEA analysis")

        source("seurat_workflow/code/pathway_analysis.R")

        # hallmark
        gmt.file       <- "seurat_workflow/data/h.all.v7.2.symbols.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'HALLMARK'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)

        # Chondroid signature
        gmt.file       <- "seurat_workflow/data/chondroid.symbol.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'Chondroid'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)

        # Cell origin
        gmt.file       <- "seurat_workflow/data/marker_genes/cell_origin/Pommier.genesymbol.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'CellOrigin'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)

        # Pece MaSC
        gmt.file       <- "seurat_workflow/data/marker_genes/cell_origin/stemchecker_Hs_MaSC_Pece.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'Pece_MaSC'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)


        # Epithelial signature
        gmt.file       <- "seurat_workflow/data/epithelial.symbol.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'Epithelial'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)

        # Mesenchymal signature
        gmt.file       <- "seurat_workflow/data/mesenchymal.symbol.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'Mesenchymal'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)


        # Stemness signature
        gmt.file       <- "seurat_workflow/data/marker_genes/cell_origin/my_stemness_signatures.gmt"
        gene.sets      <- getGmt(gmt.file)
        gene.sets.name <- 'Stemness'
        ssgsea_analysis(obj,outDir=paste0(outDir,'/ssGSEA'),gene.sets,gene.sets.name)

        #########################################################################################################
        # evaluate stemness according to stemness sigature score
        if(evaluate_stemness){ # this may take 20min

            source("seurat_workflow/code/evaluate_stemness.R")

            # for all stemness signatures
            gmt.file   <- "seurat_workflow/data/marker_genes/cell_origin/my_stemness_signatures.gmt"
            genesets   <- getGmt(gmt.file)

            lapply(names(genesets),function(x){
                    print(x)
                    genes <- geneIds(genesets[[x]])
                    gene.sets.name <- x
                    evaluate_stemness(obj,outDir=paste0(outDir,'/stemness/',gene.sets.name),genes,gene.sets.name)
            })
        }


        #########################################################################################################
        # Run StemSC

        source("seurat_workflow/code/stemSC.R")

        obj <- Run_StemSC(obj,outDir=paste0(outDir,"/StemSC"))

}

