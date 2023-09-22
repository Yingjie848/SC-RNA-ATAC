# Goal: In silico CD45 sorting
# Methods: 1) define CD45+/- identify highly variable genes; 

# rm(list=ls())

library(dplyr)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)

# // prepare RNA_data and SCT_data

prepare_data <- function(fdata="seurat_data.rds",outDir="InSilicoCD45sorting",geneExXCells=20,CD45posMinReads=1,immune_gsea=NULL,genes=NULL,gene_clustering_method='kmeans'){

    dir.create(outDir)
    print(fdata)
    obj=readRDS(fdata)

    # identify HVGs
    print("identify HVGs")
    DefaultAssay(obj) <- "SCT"
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    
    cat("Number of cells: ",nrow(obj@meta.data),"\n")

    cell_info <- obj@meta.data %>% 
                  mutate(cell=rownames(obj@meta.data),immune_gsea_score=immune_gsea)
    cell_info <- cell_info[,! colnames(cell_info) %in% c('CD45','PanCK')]

    RNA_data  = as.matrix(obj$RNA@data)

    # // create Pan Cytokeratin (Pan-CK) counts 
    print("create Pan Cytokeratin (Pan-CK) counts")
    panck_genes = rownames(RNA_data)[grepl('^KRT',rownames(RNA_data))]
    panck_genes = panck_genes[!grepl('CAP',panck_genes)]

    RNA_data_panck <- RNA_data[panck_genes,]

    df_panck  <- data.frame(cell=colnames(RNA_data_panck),PanCK=colSums(RNA_data_panck))
    cell_info <- merge(cell_info,df_panck,by='cell',all.x=T,sort=F)

    # // filter genes
    RNA_data         = RNA_data[,cell_info$cell]
    nexpressed_cells = apply(RNA_data,1,function(x) sum(x>0))
    RNA_data         = RNA_data[nexpressed_cells >= geneExXCells | rownames(RNA_data) == 'PTPRC',] # filter out genes expressed in less than X cells 

    # find highly variable genes
    print("Find 1000 highly variable genes")
    HVGs      <- VariableFeatures(obj)
    HVGs      <- head(HVGs[HVGs %in% rownames(RNA_data)], 1000)

    # extract data
    print("Extract data by genes")
    if(is.null(genes)){
        genes <- HVGs
    }
    if(!'PTPRC' %in% genes){
        genes <- c('PTPRC',genes)
    }
    print(genes[!genes %in% rownames(RNA_data)])
    RNA_data  <- RNA_data[rownames(RNA_data) %in% genes,]; print(dim(RNA_data));print('PTPRC' %in% rownames(RNA_data));
    SCT_data  <- as.matrix(obj$SCT@data)
    sum(rownames(SCT_data) %in% rownames(RNA_data)) # check how many SCT genes not in filtered data. Based on calculation, SCT should excluded genes mutated in less than 5 cells

    # divide CD45+/CD45- cells
    print("divide CD45+/CD45- cells")
    cat("cutoff:",CD45posMinReads,"\n")
    CD45_expr   = t(RNA_data)[,'PTPRC']; 
    CD45_counts = data.frame(cell=names(CD45_expr),CD45=CD45_expr) %>% 
                  arrange(desc(CD45)) %>% 
                  mutate(CD45_status=ifelse(CD45>=CD45posMinReads,'CD45+','CD45-'))

    # order cells by CD45 counts and immune gsea
    cell_info   = merge(cell_info,CD45_counts,by='cell',all.x=T,sort=F) %>% 
                  arrange(desc(CD45),desc(immune_gsea_score)) 
    cell_info.ordered_by_CD45_and_immune_gsea <- cell_info
    print("Cells ordered by seurat_clusters,desc(CD45),desc(immune_gsea_score)")
    write.table(cell_info,paste0(outDir,"/metadata_cells.txt"),quote=F,sep="\t",row.names=F)

    p1 <- ggplot(CD45_counts %>% mutate(order=1:length(CD45)),aes(x=CD45,y=..count../sum(..count..))) + geom_histogram(binwidth=1,fill='#c9838e') + ggtitle("CD45 Read frequency")
    p2 <- ggplot(CD45_counts %>% mutate(order=1:length(CD45)),aes(x=order,y=CD45)) + geom_point(color='#c9838e',size=0.5) + ggtitle("Ranked counts")
    p = ggarrange(p1,p2,ncol=2)
    pdf(paste0(outDir,"/CD45_counts.pdf"),7,4);print(p);dev.off()

    CD45pos_cells = cell_info[cell_info$CD45_status=='CD45+',]$cell
    CD45neg_cells = cell_info[cell_info$CD45_status=='CD45-',]$cell

    # // order cells by cell clusters, clusters having more CD45+ cells will be put at the beginning

    # get frequency of CD45+ cells for each cluster
    cell_cluster_CD45_freq <- cell_info %>% 
                              group_by(seurat_clusters) %>% 
                              summarise(ncells=length(cell),CD45_frequency=sum(CD45>=CD45posMinReads)) %>% 
                              mutate(CD45_cell_pct=CD45_frequency/ncells*100) %>% 
                              arrange(desc(CD45_cell_pct)) # order cell clusters by CD45+ frequency
    cell_info.ordered_by_clusters_ranked_by_CD45 <- lapply(as.vector(cell_cluster_CD45_freq$seurat_clusters),
      function(x){
        tmp = cell_info[cell_info$seurat_clusters==x,] %>% arrange(desc(CD45),desc(immune_gsea_score))
      })
    cell_info.ordered_by_clusters_ranked_by_CD45 <- do.call(rbind,cell_info.ordered_by_clusters_ranked_by_CD45) # reorder cells by cell clusters, clusters having more CD45+ cells will be put at the beginning
    cell_info <- cell_info.ordered_by_clusters_ranked_by_CD45

    # // order cells by CD45 counts and immune gsea
    cell_info.ordered_by_CD45_and_immune_gsea   = cell_info %>% 
                                                  arrange(desc(CD45)) 
    print("Cells ordered by seurat_clusters,desc(CD45),desc(immune_gsea_score), immune cell clusters")
    write.table(cell_info.ordered_by_CD45_and_immune_gsea,paste0(outDir,"/metadata_cells.txt"),quote=F,sep="\t",row.names=F)

    CD45pos_cells = cell_info.ordered_by_CD45_and_immune_gsea[cell_info.ordered_by_CD45_and_immune_gsea$CD45_status=='CD45+',]$cell
    CD45neg_cells = cell_info.ordered_by_CD45_and_immune_gsea[cell_info.ordered_by_CD45_and_immune_gsea$CD45_status=='CD45-',]$cell

    # thus cells CD45pos_cells and CD45neg_cells are ordered by 1) CD45 counts; 2) immune gsea; 3) immune clusters
    # for RNA_data and SCT_data, they are ordered by cell clusters ranked by frequency of CD45+ cells 

    # now extract data
    RNA_data         = RNA_data[,as.vector(cell_info$cell)]
    SCT_data         = t(t(SCT_data[,colnames(RNA_data)])[,rownames(RNA_data)])

    RNA_data_CD45pos = RNA_data[,CD45pos_cells]
    RNA_data_CD45neg = RNA_data[,CD45neg_cells]

    SCT_data_CD45pos = SCT_data[,CD45pos_cells]
    SCT_data_CD45neg = SCT_data[,CD45neg_cells]

    # // do clustering based on SCT_data
    print("do clustering")
    SCT_row_distance <- dist(SCT_data,method='euclidean')
    SCT_row_clust    <- hclust(SCT_row_distance,method='complete')

    SCT_col_distance <- dist(t(SCT_data),method='euclidean')
    SCT_col_clust    <- hclust(SCT_col_distance,method='complete')

    SCT_col_distance_CD45pos <- dist(t(SCT_data_CD45pos),method='euclidean')
    SCT_col_clust_CD45pos    <- hclust(SCT_col_distance_CD45pos,method='complete')

    SCT_col_distance_CD45neg <- dist(t(SCT_data_CD45neg),method='euclidean')
    SCT_col_clust_CD45neg    <- hclust(SCT_col_distance_CD45neg,method='complete')

    RNA_data <- RNA_data[SCT_row_clust$order,]
    SCT_data <- SCT_data[SCT_row_clust$order,]
    RNA_data_CD45pos <- RNA_data_CD45pos[SCT_row_clust$order,]
    RNA_data_CD45neg <- RNA_data_CD45neg[SCT_row_clust$order,]
    SCT_data_CD45pos <- SCT_data_CD45pos[SCT_row_clust$order,]
    SCT_data_CD45neg <- SCT_data_CD45neg[SCT_row_clust$order,]


    # // group genes by K-means
    if(gene_clustering_method=='kmeans'){
        print("Making gene groups by kmeans")
        library(factoextra)
        p1 <- fviz_nbclust(SCT_data, kmeans, method = "wss",k.max=15)
        p2 <- fviz_nbclust(t(SCT_data), kmeans, method = "wss",k.max=15)
        pdf(paste0(outDir,"/","nclusters_genes.kmeans.pdf"),15,7); print(p1); dev.off()
        pdf(paste0(outDir,"/","nclusters_cells.kmeans.pdf"),15,7); print(p2); dev.off()

        # Compute k-means
        set.seed(123)
        km.res.genes <- kmeans(SCT_data, 10, nstart = 1)
        km.res.cells <- kmeans(t(SCT_data), 12, nstart = 1)
        gene_clusters <- data.frame(gene_name=names(km.res.genes$cluster),gene_cluster=km.res.genes$cluster)
        write.table(gene_clusters,paste0(outDir,"/","kmeans_gene_clusters.tsv"),quote=F,sep="\t",row.names=F)
    }else if(gene_clustering_method=='hcut'){
        print("Making gene groups by hierarchical clustering")
        #library(factoextra)
        #p1 <- fviz_nbclust(SCT_data, hcut, method = "wss",k.max=15)
        #p2 <- fviz_nbclust(t(SCT_data), hcut, method = "wss",k.max=15)
        #pdf(paste0(outDir,"/","nclusters_genes.hcut.pdf"),15,7); print(p1); dev.off()
        #pdf(paste0(outDir,"/","nclusters_cells.hcut.pdf"),15,7); print(p2); dev.off()

        # Compute k-means
        hcut.res.genes <- cutree(SCT_row_clust,k=12)
        hcut.res.cells <- cutree(SCT_col_clust,k=12)
        gene_clusters <- data.frame(gene_name=names(hcut.res.genes),gene_cluster=hcut.res.genes)
        write.table(gene_clusters,paste0(outDir,"/","hcut_gene_clusters.tsv"),quote=F,sep="\t",row.names=F)
    }

    # save
    save(cell_info,cell_info.ordered_by_CD45_and_immune_gsea,cell_info.ordered_by_clusters_ranked_by_CD45,gene_clusters,
         CD45_counts,RNA_data,RNA_data_CD45pos,RNA_data_CD45neg,RNA_data_panck,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg,
         SCT_row_clust,SCT_col_clust,SCT_col_clust_CD45pos,SCT_col_clust_CD45neg,file=paste0(outDir,"/data.Rdata"))

}


InSilicoCD45sorting <- function(fdata="seurat_data.rds",outDir="InSilicoCD45sorting",geneExXCells=20,CD45posMinReads=1,test=FALSE,
                                immune_gsea=NULL,prepare_data=FALSE,gene_clustering_method='kmeans',genes=NULL,annotate_gene_clusters=TRUE,
                                order_cells_by_doublets=FALSE,order_cells_by_CD45_doublets=FALSE,order_cells_by_CD45_scDblFinder=FALSE,hgnc=NULL){

    if(! file.exists(paste0(outDir,"/data.Rdata")) | prepare_data){
      print("Data not exists, preparing")
      prepare_data(fdata,outDir,geneExXCells,CD45posMinReads,immune_gsea,genes,gene_clustering_method)
      load(paste0(outDir,"/data.Rdata"))
    }else if(file.exists(paste0(outDir,"/data.Rdata"))){
      print("Data exists")
      print(paste0(outDir,"/data.Rdata"))
      load(paste0(outDir,"/data.Rdata"))
    }


    ########################################################################################
    # annotate gene clusters defined by K-means
    ########################################################################################

    print("annotate gene clusters")

    func_annotate_gene_clusters <- function(){

        # gene sets
        gmtfile.bp        <- "~/my_databases/MSigDB/c5.bp.v7.1.entrez.gmt"
        gmtfile.kegg      <- "~/my_databases/MSigDB/c2.cp.kegg.v7.1.entrez.gmt"
        gmtfile.reactome  <- "~/my_databases/MSigDB/c2.cp.reactome.v7.1.entrez.gmt"
        gmtfile.hallmark  <- "/home/zhuy1/my_databases/MSigDB/h.all.v7.2.entrez.gmt"

        # 
        y = table(gene_clusters$gene_cluster)
        print(y)
        clusters = names(y[y>=2])
        print(clusters)

        annotated_clusters <- lapply(sort(clusters),function(cluster){
          print(cluster)
          gene_name = gene_clusters[gene_clusters$gene_cluster==cluster,]$gene_name
          #genemap  <- bitr(gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F) %>% dplyr::filter(!is.na(ENTREZID))
          genemap  <- hgnc[hgnc$gene_symbol %in% gene_name,]

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

        write.table(out_annotated_clusters,paste0(outDir,"/","kmeans_gene_clusters_annotation.tsv"),quote=F,sep="\t",row.names=F)

    }

    if(annotate_gene_clusters)
        func_annotate_gene_clusters()

    ########################################################################################
    # make heatmap 
    ########################################################################################
    print("make heatmap")

    # make heatmap annotations
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
        immune_gsea_cols     = colorRamp2( c(min(cell_info$immune_gsea), max(cell_info$immune_gsea)+0.1), c("#f5f2f2", "red"))
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

    
    # reorder cells in cell_info and RNA_data by doublets
    if(order_cells_by_doublets){

      print("order cells by doublets")

      cell_info <- cell_info %>% arrange(desc(scDblFinder.doublets),desc(scDblFinder.score),desc(hybrid_score))
      
      RNA_data  <- RNA_data[,cell_info$cell]
      SCT_data  <- SCT_data[,cell_info$cell]
      SCT_data_CD45pos <- SCT_data_CD45pos[,colnames(SCT_data_CD45pos) %in% cell_info$cell]
      SCT_data_CD45neg <- SCT_data_CD45neg[,colnames(SCT_data_CD45neg) %in% cell_info$cell]

      annotation_list <- make_heatmap_annotations(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg)
      col_ha          <- annotation_list[[1]]
      col_ha_CD45pos  <- annotation_list[[2]]
      col_ha_CD45neg  <- annotation_list[[3]]

      p <- Heatmap(log10(RNA_data+1),name="RNA counts (log10)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
      pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_doublets.pdf"),15,10); draw(p); dev.off()

      return("made heatmap ordered cells by doublets")
    }

    if(order_cells_by_CD45_doublets){

      print("order cells by CD45 and doublets")

      # add doublet to cell_info
      cell_info <- cell_info %>% mutate(Doublets=ifelse(scDblFinder.doublets=='doublet' | hybrid_score>=1,'Yes','No'))

      cell_info_CD45pos_doublets <- cell_info %>% 
                                          filter(CD45_status=='CD45+' & Doublets=='Yes')
      cell_info_CD45pos_singlets <- cell_info %>% 
                                          filter(CD45_status=='CD45+' & Doublets=='No')
      cell_info_CD45neg_doublets <- cell_info %>% 
                                          filter(CD45_status=='CD45-' & Doublets=='Yes')
      cell_info_CD45neg_singlets <- cell_info %>% 
                                          filter(CD45_status=='CD45-' & Doublets=='No')

      cell_info <- rbind(cell_info_CD45pos_doublets,cell_info_CD45pos_singlets,cell_info_CD45neg_doublets,cell_info_CD45neg_singlets)

      write.table(cell_info,paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_CD45_doublets.cell_info.txt"),quote=F,sep="\t",row.names=F)

      RNA_data  <- RNA_data[,cell_info$cell]
      RNA_data_CD45pos <- RNA_data_CD45pos[,cell_info[cell_info$CD45_status=='CD45+',]$cell]
      RNA_data_CD45neg <- RNA_data_CD45neg[,cell_info[cell_info$CD45_status=='CD45-',]$cell]
      SCT_data  <- SCT_data[,cell_info$cell]
      SCT_data_CD45pos <- SCT_data_CD45pos[,cell_info[cell_info$CD45_status=='CD45+',]$cell]
      SCT_data_CD45neg <- SCT_data_CD45neg[,cell_info[cell_info$CD45_status=='CD45-',]$cell]

      annotation_list <- make_heatmap_annotations(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg)
      col_ha          <- annotation_list[[1]]
      col_ha_CD45pos  <- annotation_list[[2]]
      col_ha_CD45neg  <- annotation_list[[3]]

      p5 <- Heatmap(log10(RNA_data_CD45pos+1),name="RNA counts (log10)",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
      p6 <- Heatmap(log10(RNA_data_CD45neg+1),name="RNA counts (log10)",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
      pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_CD45_doublets.pdf"),15,10); draw(p5+p6); dev.off()

      return("made heatmap ordered cells by CD45 and doublets")
    }

    if(order_cells_by_CD45_scDblFinder){

      print("order cells by CD45 and scDblFinder")

      # add doublet to cell_info
      cell_info <- cell_info %>% mutate(Doublets=ifelse(scDblFinder.doublets=='doublet','Yes','No'))

      cell_info_CD45pos_doublets <- cell_info %>% 
                                          filter(CD45_status=='CD45+' & Doublets=='Yes')
      cell_info_CD45pos_singlets <- cell_info %>% 
                                          filter(CD45_status=='CD45+' & Doublets=='No')
      cell_info_CD45neg_doublets <- cell_info %>% 
                                          filter(CD45_status=='CD45-' & Doublets=='Yes')
      cell_info_CD45neg_singlets <- cell_info %>% 
                                          filter(CD45_status=='CD45-' & Doublets=='No')

      cell_info <- rbind(cell_info_CD45pos_doublets,cell_info_CD45pos_singlets,cell_info_CD45neg_doublets,cell_info_CD45neg_singlets)

      write.table(cell_info,paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_CD45_scDblFinder.cell_info.txt"),quote=F,sep="\t",row.names=F)

      RNA_data  <- RNA_data[,cell_info$cell]
      RNA_data_CD45pos <- RNA_data_CD45pos[,cell_info[cell_info$CD45_status=='CD45+',]$cell]
      RNA_data_CD45neg <- RNA_data_CD45neg[,cell_info[cell_info$CD45_status=='CD45-',]$cell]
      SCT_data  <- SCT_data[,cell_info$cell]
      SCT_data_CD45pos <- SCT_data_CD45pos[,cell_info[cell_info$CD45_status=='CD45+',]$cell]
      SCT_data_CD45neg <- SCT_data_CD45neg[,cell_info[cell_info$CD45_status=='CD45-',]$cell]

      annotation_list <- make_heatmap_annotations(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg)
      col_ha          <- annotation_list[[1]]
      col_ha_CD45pos  <- annotation_list[[2]]
      col_ha_CD45neg  <- annotation_list[[3]]

      p5 <- Heatmap(log10(RNA_data_CD45pos+1),name="RNA counts (log10)",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
      p6 <- Heatmap(log10(RNA_data_CD45neg+1),name="RNA counts (log10)",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
      pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_CD45_scDblFinder.pdf"),15,10); draw(p5+p6); dev.off()

      return("made heatmap ordered cells by CD45 and scDblFinder")
    }

    # make heatmap annotations
    annotation_list <- make_heatmap_annotations(cell_info,SCT_data,SCT_data_CD45pos,SCT_data_CD45neg)
    col_ha          <- annotation_list[[1]]
    col_ha_CD45pos  <- annotation_list[[2]]
    col_ha_CD45neg  <- annotation_list[[3]]

    # make heatmap, ordered by cell clusters
    #p <- Heatmap(SCT_data,name="SCT",column_title=paste0(ncol(SCT_data)," cells"),top_annotation = col_ha,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    #pdf(paste0(outDir,"/","heatmap_SCTdata_kmeansByGene.ordered_by_cell_clusters.pdf"),15,10); draw(p); dev.off()

    p <- Heatmap(log10(RNA_data+1),name="RNA counts (log10)",column_title=paste0(ncol(RNA_data)," cells"),top_annotation = col_ha,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_cell_clusters.pdf"),15,10); draw(p); dev.off()


    ########################################################################################
    # make heatmap by separating CD45+/CD45- cells
    ########################################################################################

    # RNA, separate cells by CD45
    p5 <- Heatmap(log10(RNA_data_CD45pos+1),name="RNA counts (log10)",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    p6 <- Heatmap(log10(RNA_data_CD45neg+1),name="RNA counts (log10)",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.ordered_by_CD45.pdf"),15,10); draw(p5+p6); dev.off()

    p7 <- Heatmap(log10(RNA_data_CD45pos+1),name="RNA counts (log10)",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = SCT_col_clust_CD45pos,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    p8 <- Heatmap(log10(RNA_data_CD45neg+1),name="RNA counts (log10)",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = SCT_col_clust_CD45neg,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    pdf(paste0(outDir,"/","heatmap_RNAdata_kmeansByGene.clustered_by_cells.pdf"),15,10); draw(p7+p8); dev.off()


    # SCT, separate cells by CD45
    #SCT_data_CD45pos[SCT_data_CD45pos>3] = 3
    #SCT_data_CD45neg[SCT_data_CD45neg>3] = 3

    #p1 <- Heatmap(SCT_data_CD45pos,name="SCT",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    #p2 <- Heatmap(SCT_data_CD45neg,name="SCT",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = F,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    #pdf(paste0(outDir,"/","heatmap_SCTdata_kmeansByGene.ordered_by_CD45.pdf"),15,10); draw(p1+p2); dev.off()

    #p3 <- Heatmap(SCT_data_CD45pos,name="SCT",column_title=paste0("CD45+ ",ncol(SCT_data_CD45pos)," cells"),top_annotation = col_ha_CD45pos,cluster_rows = F, cluster_columns = SCT_col_clust_CD45pos,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    #p4 <- Heatmap(SCT_data_CD45neg,name="SCT",column_title=paste0("CD45- ",ncol(SCT_data_CD45neg)," cells"),top_annotation = col_ha_CD45neg,cluster_rows = F, cluster_columns = SCT_col_clust_CD45neg,row_split=gene_clusters$gene_cluster,show_row_dend=T,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    #pdf(paste0(outDir,"/","heatmap_SCTdata_kmeansByGene.clustered_by_cells.pdf"),15,10); draw(p3+p4); dev.off()


    ########################################################################################
    # make cell-cell correlation heatmap
    ########################################################################################

    
    # make heatmap annotations
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

    row_ha <- make_heatmap_row_annotation_for_cells(cell_info)

    p <- Heatmap(mat_cor,name="Correlation Coefficient",column_title=paste0(ncol(mat_cor)," cells"),left_annotation = row_ha, top_annotation = col_ha,cluster_rows = F, cluster_columns = F,show_row_dend=F,show_row_names = F,show_column_names = F,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))
    pdf(paste0(outDir,"/","heatmap_RNAdata_correlation.pdf"),14,8); draw(p); dev.off()


}


