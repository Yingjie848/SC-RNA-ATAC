
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


clustering_analysis <- function(obj,outDir,npc=50,resolution = 1, stop_after_marker_genes = FALSE,regress_out_cellcycle=FALSE ){

        
        source("seurat_workflow/code/annotate_cells.R")
        

        dir.create(outDir)


        #########################################################################################################
        # Add CD45 and PanCK counts
        print("Add CD45 and PanCK counts")
        obj <- annotate_cells(obj)


        #########################################################################################################
        # scoring cell cycle
        # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
        # segregate this list into markers of G2/M phase and markers of S phase
        s.genes   <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes

        obj <- NormalizeData(obj,assay='RNA')
        obj <- ScaleData(obj)
        obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

        # scale data and regress out variables
        if(regress_out_cellcycle){
          obj <- ScaleData(obj,vars.to.regress=c("S.Score", "G2M.Score"))
        }

        #########################################################################################################
        # do clustering

        if(!regress_out_cellcycle){
            cat("Do clustering without regressing out cell cycle\n")
            obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE) %>%
                      RunPCA(npcs = npc, verbose = FALSE) %>%
                      RunUMAP(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      RunTSNE(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindNeighbors(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindClusters(resolution = resolution, verbose = FALSE)
        }else{
            cat("Do clustering with regressing out cell cycle\n")
            obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE,vars.to.regress=c("S.Score", "G2M.Score")) %>%
                      RunPCA(npcs = npc, verbose = FALSE) %>%
                      RunUMAP(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      RunTSNE(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindNeighbors(reduction = "pca", dims = 1:npc, verbose = FALSE) %>%
                      FindClusters(resolution = resolution, verbose = FALSE)  
            obj <- ScaleData(obj, assay='RNA', vars.to.regress = c("S.Score", "G2M.Score"))
        }
        

        saveRDS(obj,file=paste0(outDir,"/seurat_rna.rds"))


        #########################################################################################################
        # QC cells after clustering
        cat("Create QC plots\n")
        p1 <- VlnPlot(obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 3,log = TRUE, pt.size = 0.5) + NoLegend()
        p2 <- VlnPlot(obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","MALAT1","percent.ribo","percent.hb"), ncol = 3,log = FALSE, pt.size = 0.5) + NoLegend()

        pdf(paste0(outDir,"/qc_metrics.with_clustering.log.pdf"),15,7);print(p1);dev.off()
        pdf(paste0(outDir,"/qc_metrics.with_clustering.raw.pdf"),15,7);print(p2);dev.off()


        #########################################################################################################
        # visualize clustering based on gene expression 
        cat("Create dimentional reduction plots\n")

        p <- DimPlot(obj, reduction = "pca", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("PCA")

        pdf(paste0(outDir,"/DimPlot_pca.pdf"),5,5);print(p);dev.off()

        p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("UMAP")

        pdf(paste0(outDir,"/DimPlot_umap.pdf"),5,5);print(p);dev.off()

        p <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("tSNE")

        pdf(paste0(outDir,"/DimPlot_tsne.pdf"),5,5);print(p);dev.off()

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
        DefaultAssay(obj) <- 'RNA'

        pdf(paste0(outDir,'/FeaturePlot.pdf'))
        print(FeaturePlot(obj,features='nCount_RNA'))
        print(FeaturePlot(obj,features='nFeature_RNA'))
        print(FeaturePlot(obj,features='percent.mt'))
        print(FeaturePlot(obj,features='MALAT1'))
        print(FeaturePlot(obj,features='percent.ribo'))
        print(FeaturePlot(obj,features='percent.hb'))
        print(FeaturePlot(obj,features='CD45'))
        print(FeaturePlot(obj,features='PanCK'))
        print(FeaturePlot(obj,features='hybrid_score'))
        print(FeaturePlot(obj,features='cxds_score'))
        print(FeaturePlot(obj,features='bcds_score'))
        print(FeaturePlot(obj,features='scDblFinder.score'))
        print(FeaturePlot(obj,features='S.Score'))
        print(FeaturePlot(obj,features='G2M.Score'))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Mesenchymal.pdf'))
        print(FeaturePlot(obj,features='CDH2'))
        print(FeaturePlot(obj,features='VIM'))
        print(FeaturePlot(obj,features='FN1'))
        print(FeaturePlot(obj,features='ITGB6'))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Epithelial.pdf'))
        print(FeaturePlot(obj,features='CDH1'))
        print(FeaturePlot(obj,features='DSP'))
        print(FeaturePlot(obj,features='OCLN'))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.EMT.pdf'))
        print(FeaturePlot(obj,features='SNAI1')) # don't have this gene
        print(FeaturePlot(obj,features='SNAI2'))
        print(FeaturePlot(obj,features='ZEB1'))
        print(FeaturePlot(obj,features='ZEB2'))
        print(FeaturePlot(obj,features='TWIST1'))
        print(FeaturePlot(obj,features='TWIST2'))
        dev.off()


        #########################################################################################################
        if(stop_after_marker_genes){
          return("Stopped after marker genes")
        }


        #########################################################################################################
        # // Finding differentially expressed features (cluster biomarkers)

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

        write.table(markers,paste0(outDir,"/DE_analysis/DE_genes.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_top50,paste0(outDir,"/DE_analysis/DE_genes_top50.txt"),quote=F,sep="\t",row.names=F)

        # markers, only the positive genes
        markers_pos <- markers %>% dplyr::filter(avg_log2FC>0 & p_val_adj < 0.001)

        markers_pos_top50 <-markers_pos %>%
          dplyr::group_by(cluster) %>%
          dplyr::top_n(n = 50, wt = avg_log2FC)

        write.table(markers_pos,paste0(outDir,"/DE_analysis/posDE_genes.txt"),quote=F,sep="\t",row.names=F)
        write.table(markers_pos_top50,paste0(outDir,"/DE_analysis/posDE_genes_top50.txt"),quote=F,sep="\t",row.names=F)


        #########################################################################################################
        # // annotate up-regulated genes for each cluster

        cat("Annotate up-regulated genes for each cluster\n")

        source("seurat_workflow/code/pathway.R")

        # load hgnc gene id map
        hgnc = fread("seurat_workflow/data/hgnc_20220106.txt")
        hgnc = hgnc[,c('Approved symbol','NCBI Gene ID')]
        colnames(hgnc) = c('gene_symbol','ENTREZID')

        # gene sets
        gmtfile.bp        <- "seurat_workflow/data/c5.bp.v7.1.entrez.gmt"
        gmtfile.kegg      <- "seurat_workflow/data/c2.cp.kegg.v7.1.entrez.gmt"
        gmtfile.reactome  <- "seurat_workflow/data/c2.cp.reactome.v7.1.entrez.gmt"
        gmtfile.hallmark  <- "seurat_workflow/data/h.all.v7.2.entrez.gmt"

        # 
        y = table(markers_pos$cluster); print(y)
        clusters = names(y[y>=5]); print(clusters)

        annotated_clusters <- lapply(sort(clusters),function(cluster){

          cat(cluster,': ',y[cluster],"\n")
          genes = markers_pos[markers_pos$cluster==cluster,]$gene
          genes = genes[!grepl('^RP[SL]',genes)]
          genes = genes[!grepl('^MT-',genes)]
          #genemap  <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F) %>% dplyr::filter(!is.na(ENTREZID))
          genemap  <- hgnc[hgnc$gene_symbol %in% genes,]
          res_gobp <- pathwayAnnotator.gmt(genemap$ENTREZID,gmtfile=gmtfile.bp)
          res_kegg <- pathwayAnnotator.gmt(genemap$ENTREZID,gmtfile=gmtfile.kegg)
          res_reactome  <- pathwayAnnotator.gmt(genemap$ENTREZID,gmtfile=gmtfile.reactome)
          res_hallmark  <- pathwayAnnotator.gmt(genemap$ENTREZID,gmtfile=gmtfile.hallmark)

          res <- data.frame()
          if(!is.null(res_gobp)){
            res <- rbind(res,as.data.frame(res_gobp))
          }
          if(!is.null(res_kegg)){
            res <- rbind(res,as.data.frame(res_kegg))
          }
          if(!is.null(res_reactome)){
            res <- rbind(res,as.data.frame(res_reactome))
          }
          if(!is.null(res_hallmark)){
            res <- rbind(res,as.data.frame(res_hallmark))
          }
          res <- res %>% 
                  arrange(p.adjust) %>% 
                  mutate(gene_cluster=cluster)
        })

        out_annotated_clusters = do.call(rbind,annotated_clusters)

        write.table(out_annotated_clusters,paste0(outDir,"/DE_analysis/","DEG_annotation.tsv"),quote=F,sep="\t",row.names=F)


        #########################################################################################################
        print("Save obj")
        saveRDS(obj,file=paste0(outDir,"/seurat_rna.rds"))


        #########################################################################################################
        # make heatmap

        print("Make heatmap for top 50 up-regulated DEGs")

        source("seurat_workflow/code/make_heatmap.R")

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top50_posDEGs_logcounts"),
                 genes=unique(markers_pos_top50$gene),
                 show_row_names = F,
                 normalize_to_logcounts = TRUE)

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top50_posDEGs_counts"),
                 genes=unique(markers_pos_top50$gene),
                 show_row_names = F,
                 normalize_to_logcounts = FALSE)

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
        # In silico CD45 sorting

        #print("In silico CD45 sorting")

        #source("seurat_workflow/code/InSilicoCD45sorting.R")

        # make in silico CD45 sorting heatmap

        #InSilicoCD45sorting(fdata=paste0(outDir,"/seurat_rna.rds"),
        #                    outDir=paste0(outDir,"/InSilicoCD45sorting_DEGs"),
        #                    geneExXCells=20,
        #                    CD45posMinReads=2,
        #                    immune_gsea=obj@meta.data$immune_gsea,
        #                    test=FALSE,
        #                    prepare_data=TRUE,
        #                    gene_clustering_method = 'kmeans',
        #                    annotate_gene_clusters=TRUE,
        #                    genes=unique(markers_pos_top50$gene),
        #                    hgnc=hgnc)


        #InSilicoCD45sorting(fdata=paste0(outDir,"/seurat_rna.rds"),
        #                    outDir=paste0(outDir,"/InSilicoCD45sorting_DEGs"),
        #                    geneExXCells=20,
        #                    CD45posMinReads=2,
        #                    immune_gsea=obj@meta.data$immune_gsea,
        #                    test=FALSE,
        #                    prepare_data=FALSE,
        #                    gene_clustering_method = 'kmeans',
        #                    genes=unique(markers_pos_top50$gene),
        #                    annotate_gene_clusters=FALSE,
        #                    hgnc=hgnc,
        #                    order_cells_by_CD45_doublets = TRUE)

        #InSilicoCD45sorting(fdata=paste0(outDir,"/seurat_rna.rds"),
        #                    outDir=paste0(outDir,"/InSilicoCD45sorting_DEGs"),
        #                    geneExXCells=20,
        #                    CD45posMinReads=2,
        #                    immune_gsea=obj@meta.data$immune_gsea,
        #                    test=FALSE,
        #                    prepare_data=FALSE,
        #                    gene_clustering_method = 'kmeans',
        #                    genes=unique(markers_pos_top50$gene),
        #                    annotate_gene_clusters=FALSE,
        #                    hgnc=hgnc,
        #                    order_cells_by_CD45_scDblFinder = TRUE)

        #InSilicoCD45sorting(fdata=paste0(outDir,"/seurat_rna.rds"),
        #                    outDir=paste0(outDir,"/InSilicoCD45sorting_HVGs"),
        #                    geneExXCells=20,
        #                    CD45posMinReads=2,
        #                    immune_gsea=obj@meta.data$immune_gsea,
        #                    test=FALSE,
        #                    prepare_data=TRUE,
        #                    gene_clustering_method = 'kmeans',
        #                    hgnc=hgnc,
        #                    annotate_gene_clusters=TRUE)

}