# annotate genes to pathways

library(clusterProfiler)
library(org.Hs.eg.db)
library(magrittr)
library(ReactomePA)
library(stringr)
library(enrichplot)
library(pathview)
library(ggplot2)
library(dplyr)

pathwayAnnotator.gmt <- function(genes,pval_thres=1,qval_thres=1,gmtfile=NULL,selected_pathways=NULL,clean_desc=FALSE,gene_symbol=FALSE){

    # filter ORA results
    filter_ORA <- function(x){
            x.result <- x@result
            x.incutoff <- x.result[x.result$p.adjust <= x@pvalueCutoff,]
            x.incutoff <- x.incutoff[!is.na(x.incutoff$ID),]
            x@result <- x.incutoff
            return(x)
    }
  
    #genemap <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
    #geneID      <- genemap$ENTREZID
    geneID = genes

    gmt <- read.gmt(gmtfile)
    if( !is.null(selected_pathways) ){
        gmt = gmt[gmt$ont %in% selected_pathways,]
    }
    
    enrichRes <- enricher(gene=geneID,
                    pvalueCutoff=pval_thres,
                    qvalueCutoff = qval_thres,
                    minGSSize = 1,
                    pAdjustMethod = "BH",
                    TERM2GENE=gmt)
    print("...")
    
    if(!is.null(enrichRes)){
        if(!gene_symbol){
            enrichRes <- setReadable(enrichRes,org.Hs.eg.db,keyType = "ENTREZID")
        }

        if( isTRUE(clean_desc) ){
            enrichRes@result$Description <- gsub("^GO_","",
                                        gsub("^KEGG_","",
                                        gsub("^REACTOME_","",
                                        enrichRes@result$Description)))
        }
        #enrichRes@result$Description <- wrap.labels(enrichRes@result$Description,35)
        enrichRes <- filter_ORA(enrichRes)

        return(enrichRes@result)
    }
  
}

annotate_posDEG_to_pathway <- function(markers_pos,outDir,topn=5){

        library(data.table)

        # load hgnc gene id map
        hgnc = fread("seurat_workflow/data/hgnc_20220106.txt")
        hgnc = hgnc[,c('Approved symbol','NCBI Gene ID')]
        colnames(hgnc) = c('gene_symbol','ENTREZID')

        # gene sets
        gmtfile        <- "seurat_workflow/data/msigdb_bp_kegg_reactome_h.symbol.gmt"

        # 
        y = table(markers_pos$cluster); print(y)
        clusters = names(y[y>=5]); print(clusters)

        annotated_clusters <- lapply(sort(clusters),function(cluster){

                cat(cluster,': ',y[cluster],"\n")
                genes = markers_pos[markers_pos$cluster==cluster,]$gene
                genes = genes[!grepl('^RP[SL]',genes)]
                genes = genes[!grepl('^MT-',genes)]
                #genemap  <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F) %>% dplyr::filter(!is.na(ENTREZID))
                #genemap  <- hgnc[hgnc$gene_symbol %in% genes,]
                res <- pathwayAnnotator.gmt(genes,gmtfile=gmtfile,gene_symbol = TRUE)
                res <- res %>% 
                        arrange(p.adjust) %>% 
                        mutate(gene_cluster=cluster)
        })

        out_annotated_clusters = do.call(rbind,annotated_clusters)

        write.table(out_annotated_clusters,paste0(outDir,"/DE_analysis/","DEG_annotation.tsv"),quote=F,sep="\t",row.names=F)

        # ------------------------------------------------------------------------
        # make pathway FDR heatmap by clusters

        out_annotated_clusters <- read.table(paste0(outDir,"/DE_analysis/","DEG_annotation.tsv"),header=T)

        # take top 5 pathway for each cluster
        top5_pathways <- out_annotated_clusters %>% group_by(gene_cluster) %>% arrange(p.adjust) %>% filter(p.adjust<0.05) %>% slice_head(n=topn)

        selected_pathways <- out_annotated_clusters %>% filter(Description %in% top5_pathways$Description)

        library(reshape2)

        pathway_adjp <- acast(selected_pathways,Description~gene_cluster,value.var="p.adjust")  
        pathway_adjp[is.na(pathway_adjp)] <- 1
        pathway_adjp_symbol <- pathway_adjp
        pathway_adjp_symbol[pathway_adjp<=0.0001] <- '****'
        pathway_adjp_symbol[pathway_adjp<=0.001 & pathway_adjp>0.0001] <- '***'
        pathway_adjp_symbol[pathway_adjp<=0.01 & pathway_adjp>0.001] <- '**'
        pathway_adjp_symbol[pathway_adjp<=0.05 & pathway_adjp>0.01] <- '*'
        pathway_adjp_symbol[pathway_adjp>0.05] <- 'ns'

        col_fun <- c('ns'='#d2dce2','*'='#107C81','**'='#1C864E','***'='#f3964f','****'='#f1643c')

        library(ComplexHeatmap)
        legend_params <- list(title = "FDR Sig.", col_fun= col_fun, border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), legend_height = grid::unit(4, "mm"), legend_width = grid::unit(6, "mm"), grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

        heatmap_width=length(unique(out_annotated_clusters$gene_cluster)) * 0.3
        heatmap_heigth=nrow(pathway_adjp) * 0.3

        pathway_clustering <- hclust(dist(pathway_adjp)) # clustering pathways using FDR

        p <- Heatmap(pathway_adjp_symbol,name="Top sig. pathways per cluster",col=col_fun,cluster_rows=pathway_clustering,cluster_columns=F,row_dend_side="right", row_names_side="left",column_names_side="top",heatmap_legend_param = legend_params,width=unit(heatmap_width,"inch"),height=unit(heatmap_heigth,"inch"))
        
        pdf(paste0(outDir,"/DE_analysis/","DEG_annotation.top_pathways_heatmap.top",topn,".pdf"),20,15);print(p);dev.off()

}

# annotate to bp,kegg,reactome,hallmark separately, then combine
annotate_posDEG_to_pathway_v1 <- function(markers_pos,outDir){

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

        # ------------------------------------------------------------------------
        # make pathway FDR heatmap by clusters

        out_annotated_clusters <- read.table(paste0(outDir,"/DE_analysis/","DEG_annotation.tsv"),header=T)

        # take top 5 pathway for each cluster
        top5_pathways <- out_annotated_clusters %>% group_by(gene_cluster) %>% arrange(p.adjust) %>% filter(p.adjust<0.05) %>% slice_head(n=5)

        selected_pathways <- out_annotated_clusters %>% filter(Description %in% top5_pathways$Description)

        library(reshape2)

        pathway_adjp <- acast(selected_pathways,Description~gene_cluster,value.var="p.adjust")  
        pathway_adjp[is.na(pathway_adjp)] <- 1
        pathway_adjp_symbol <- pathway_adjp
        pathway_adjp_symbol[pathway_adjp<=0.0001] <- '****'
        pathway_adjp_symbol[pathway_adjp<=0.001 & pathway_adjp>0.0001] <- '***'
        pathway_adjp_symbol[pathway_adjp<=0.01 & pathway_adjp>0.001] <- '**'
        pathway_adjp_symbol[pathway_adjp<=0.05 & pathway_adjp>0.01] <- '*'
        pathway_adjp_symbol[pathway_adjp>0.05] <- 'ns'

        col_fun <- c('ns'='#d2dce2','*'='#107C81','**'='#1C864E','***'='#f3964f','****'='#f1643c')

        library(ComplexHeatmap)
        legend_params <- list(title = "FDR Sig.", col_fun= col_fun, border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), legend_height = grid::unit(4, "mm"), legend_width = grid::unit(6, "mm"), grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

        heatmap_width=length(unique(out_annotated_clusters$gene_cluster)) * 0.3
        heatmap_heigth=nrow(pathway_adjp) * 0.3

        pathway_clustering <- hclust(dist(pathway_adjp)) # clustering pathways using FDR

        p <- Heatmap(pathway_adjp_symbol,name="Top5 sig. pathways per cluster",col=col_fun,cluster_rows=pathway_clustering,cluster_columns=F,row_dend_side="right", row_names_side="left",column_names_side="top",heatmap_legend_param = legend_params,width=unit(heatmap_width,"inch"),height=unit(heatmap_heigth,"inch"))
        
        pdf(paste0(outDir,"/DE_analysis/","DEG_annotation.top5_pathways_heatmap.pdf"),20,15);print(p);dev.off()

}