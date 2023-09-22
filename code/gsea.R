# using singleseqgset package to calculate logFC and do GSEA comparing one cluster to all other clusters
# https://arc85.github.io/singleseqgset/articles/singleseqgset.html 
library(Seurat)
library(singleseqgset)
library(GSEABase)
library(clusterProfiler)
library(dplyr)

col_fun <- c('ns'='#d2dce2','*'='#107C81','**'='#1C864E','***'='#f3964f','****'='#f1643c')

# make barplot for GSEA results
# Columns required in gsea_results: Description, NES, p.adjust,cluster
make_gsea_barplot <- function(gsea_results){

        source("seurat_workflow/code/lib_clusterProfiler.R")

        p <- ggplot(data=gsea_results,aes(x=Description,y=NES,fill=sigSymbol(p.adjust),width=0.6)) + geom_bar(stat="identity") +
                                xlab("") + ylab("Normalized enrichment score") +
                                # ylim(ymin,ymax) +
                                #labs(fill="p.adjust") +
                                #ggtitle(title) +
                                coord_flip() +
                                theme_bw() +
                                #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
                                scale_fill_manual("FDR",values=col_fun) +
                                theme(
                                        axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=10), 
                                        axis.text.y=element_text(color="black",size=6), 
                                        text = element_text(size=12),
                                        axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                        axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
                                        axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                        axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
                                        axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                        axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                        axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                        axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
                                        #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
                                        legend.position="right",
                                        #legend.title=element_blank(),
                                        plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
                                        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
                                ) + 
                                facet_wrap(~cluster,nrow=1)
}

make_gsea_heatmap <- function(gsea_results){
        
        library(reshape2)
        library(ComplexHeatmap)
        
        pathway_adjp <- acast(gsea_results,Description~cluster,value.var="p.adjust") 
        pathway_adjp[is.na(pathway_adjp)] <- 1
        pathway_adjp_symbol <- pathway_adjp
        pathway_adjp_symbol[pathway_adjp<=0.0001] <- '****'
        pathway_adjp_symbol[pathway_adjp<=0.001 & pathway_adjp>0.0001] <- '***'
        pathway_adjp_symbol[pathway_adjp<=0.01 & pathway_adjp>0.001] <- '**'
        pathway_adjp_symbol[pathway_adjp<=0.05 & pathway_adjp>0.01] <- '*'
        pathway_adjp_symbol[pathway_adjp>0.05] <- 'ns'

        legend_params <- list(title = "FDR Sig.", col_fun= col_fun, border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), legend_height = grid::unit(4, "mm"), legend_width = grid::unit(6, "mm"), grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

        heatmap_width=ncol(pathway_adjp) * 0.3
        heatmap_heigth=nrow(pathway_adjp) * 0.3

        pathway_clustering <- hclust(dist(pathway_adjp)) # clustering pathways using FDR

        p1 <- Heatmap(pathway_adjp_symbol,name="Top5 sig. pathways per cluster",col=col_fun,#cluster_rows=pathway_clustering,
                        cluster_columns=F,#row_dend_side="right", 
                        show_row_dend = FALSE,row_names_side="left",column_names_side="top",heatmap_legend_param = legend_params,width=unit(heatmap_width,"inch"),height=unit(heatmap_heigth,"inch"))


        pathway_nes <- acast(gsea_results,Description~cluster,value.var="NES") 
        pathway_nes[is.na(pathway_nes)] <- 0

        legend_params <- list(title = "Enrichment Score", border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), legend_height = grid::unit(20, "mm"), legend_width = grid::unit(6, "mm"), grid_height = grid::unit(20, "mm"), grid_width = grid::unit(8, "mm"))

        p2 <- Heatmap(pathway_nes,cluster_rows=pathway_clustering,cluster_columns=F,row_dend_side="right", show_row_dend = TRUE,row_names_side="left",column_names_side="top",heatmap_legend_param = legend_params,width=unit(heatmap_width,"inch"),height=unit(heatmap_heigth,"inch"),
                        cell_fun = function(j, i, x, y, w, h, fill) {
                                        if(pathway_adjp[i, j] < 0.05) {
                                                grid.text("*", x, y)
                                        }
                        }
               )

        p2


}        

make_gsea_plots <- function(combined_gsea_res,outDir,gene.sets.name){

        # // make facets barplot for all clusters together

        # first, let's create gene sets what we want to plot
        # 1) all gene sets, if gene sets < 60
        # 2) top 5 gene sets for each cluster, whatever up- or down-regulated, or significance
        # 3) top 5 up-regulated gene sets, whatever significance
        # 4) top 5 up-regulated gene sets, only significant

        all_pathways = combined_gsea_res

        sig_pathways <- combined_gsea_res %>% group_by(cluster) %>% arrange(p.adjust) %>% filter(p.adjust<0.05)
        sig_pathways <- combined_gsea_res %>% filter(Description %in% sig_pathways$Description)
        
        top5_pathways <- combined_gsea_res %>% group_by(cluster) %>% arrange(p.adjust) %>% slice_head(n=5)
        top5_pathways <- combined_gsea_res %>% filter(Description %in% top5_pathways$Description)

        top5_pathways_up <- combined_gsea_res %>% group_by(cluster) %>% arrange(p.adjust) %>% filter(NES>0) %>% slice_head(n=5)
        top5_pathways_up <- combined_gsea_res %>% filter(Description %in% top5_pathways_up$Description)

        top5_pathways_up_sig <- combined_gsea_res %>% group_by(cluster) %>% arrange(p.adjust) %>% filter(NES>0 & p.adjust<0.05) %>% slice_head(n=5)
        top5_pathways_up_sig <- combined_gsea_res %>% filter(Description %in% top5_pathways_up_sig$Description)

        # ------------------------------------------------------------------------------------------------------------------
        # if gene sets less than 60, plot all gene sets

        clusters <- unique(combined_gsea_res$cluster)

        if(length(unique(combined_gsea_res$Description))<60){

                p <- make_gsea_barplot(all_pathways)

                np <- length(unique(all_pathways$Description))
                width=5+length(clusters)*1
                height=3+np*0.2
                pdf(paste0(outDir,"/",gene.sets.name,"/gsea_all_genesets.barplot.pdf"),width=width,height = height);print(p);dev.off()

                p <- make_gsea_heatmap(all_pathways)

                pdf(paste0(outDir,"/",gene.sets.name,"/gsea_all_genesets.heatmap.pdf"),20,15);print(p);dev.off()                

        }

        if(length(unique(sig_pathways$Description))<60){

                p <- make_gsea_barplot(sig_pathways)

                np <- length(unique(sig_pathways$Description))
                width=5+length(clusters)*1
                height=3+np*0.2
                pdf(paste0(outDir,"/",gene.sets.name,"/gsea_sig_genesets.barplot.pdf"),width=width,height = height);print(p);dev.off()

                p <- make_gsea_heatmap(sig_pathways)

                pdf(paste0(outDir,"/",gene.sets.name,"/gsea_sig_genesets.heatmap.pdf"),20,15);print(p);dev.off()                

        }

        # ------------------------------------------------------------------------------------------------------------------
        # top 5 pathways for plotting

        p <- make_gsea_barplot(top5_pathways)

        np <- length(unique(top5_pathways$Description))
        width=5+length(clusters)*1
        height=3+np*0.2
        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_genesets.barplot.pdf"),width=width,height = height);print(p);dev.off()

        p <- make_gsea_heatmap(top5_pathways)

        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_genesets.heatmap.pdf"),20,15);print(p);dev.off()                


        # ------------------------------------------------------------------------------------------------------------------
        # top 5 up-regulated pathways

        p <- make_gsea_barplot(top5_pathways_up)

        np <- length(unique(top5_pathways_up$Description))
        width=5+length(clusters)*1
        height=3+np/length(clusters)*0.3
        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_up-regulated_genesets.barplot.pdf"),width=width,height = height);print(p);dev.off()

        p <- make_gsea_heatmap(top5_pathways_up)

        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_up-regulated_genesets.heatmap.pdf"),20,15);print(p);dev.off()

        # ------------------------------------------------------------------------------------------------------------------
        # top 5 up-regulated sig. pathways

        p <- make_gsea_barplot(top5_pathways_up_sig)

        np <- length(unique(top5_pathways_up_sig$Description))
        width=5+length(clusters)*1
        height=3+np*0.2
        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_up-regulated_sig_genesets.barplot.pdf"),width=width,height = height);print(p);dev.off()

        p <- make_gsea_heatmap(top5_pathways_up_sig)
        pdf(paste0(outDir,"/",gene.sets.name,"/gsea_top5_up-regulated_sig_genesets.heatmap.pdf"),20,15);print(p);dev.off()

}

# gene.sets can be loaded from gmt file by read.gmt 
do_gsea <- function(obj,outDir, gene.sets=NULL,gene.sets.name=NULL,plot_nterms=NULL){

        print(paste0(outDir,'/',gene.sets.name))

        source("seurat_workflow/code/lib_clusterProfiler.R")

        dir.create(paste0(outDir,'/',gene.sets.name),recursive = T)

        # create a vector for clusters named by cells 
        cell_clusters <- obj@meta.data$seurat_clusters
        names(cell_clusters) <- rownames(obj@meta.data)

        obj <- NormalizeData(obj,assay='RNA') # normalize data

        expr.mat <- as.matrix(obj$RNA@data); print(dim(expr.mat))

        # exclude ribosomal and mitochondrial genes
        genes = rownames(expr.mat)
        genes_RP = !grepl('^RP[SL]',genes)
        genes_MT = !grepl('^MT-',genes)

        expr.mat <- expr.mat[genes_RP & genes_MT,]; print(dim(expr.mat))

        # // First, we will calculate our metric for the gene set enrichment test, which in this case is log fold change between clusters for all genes. 
        # We choose to use the log-normalized data that have been corrected for library size. This is stored in the “@data” slot of Seurat.

        logfc.data <- logFC(cluster.ids=cell_clusters,expr.mat=expr.mat)

        # The logFC function returns a list of length 2, containing the cells in each cluster and the log fold changes in genes between the cluster of interest and all other clusters. 
        # This data is required for the next step, where we calcuate enrichment scores and p values.


        # hallmark
        #gmt.file       <- "seurat_workflow/data/h.all.v7.2.symbols.gmt"
        #gene.sets      <- getGmt(gmt.file)

        # gse.res <- wmw_gsea(expr.mat=expr.mat,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)


        # // Do GSEA using clusterProfilter
    
        #gene.sets <- read.gmt(gmt.file)

        clusters <- names(logfc.data$log.fc.cluster)

        combined_gsea_res <- lapply(clusters,function(cluster){

                print(cluster)

                geneList <- sort(logfc.data$log.fc.cluster[[cluster]],decreasing = T)
                geneList <- geneList[!is.infinite(geneList)]
                
                ggmt <- GSEA(geneList,TERM2GENE=gene.sets,minGSSize = 10, maxGSSize = 1000,pvalueCutoff = 1, eps=0,verbose=FALSE) # previous minGSSize = 1, maxGSSize = 1000
                ggmt@result$cluster <- cluster
                
                plot_nterms <- ifelse(is.null(plot_nterms),20,plot_nterms)  # already ordered by pvalue 
                p <- suppressMessages(my.barplot.gsea.enrichment_score(ggmt,title=paste0('Cluster ',cluster),n=plot_nterms))

                height=3+plot_nterms*0.1

                pdf(paste0(outDir,"/",gene.sets.name,"/gsea_cluster_",cluster,".barplot.pdf"),width=10,height = height);print(p);dev.off()

                ggmt@result

        })

        combined_gsea_res         <- do.call(rbind,combined_gsea_res)

        # Order clusters for visualization
        # remotes::install_github("jmw86069/jamba")
        library(jamba)
        clusters <- clusters[mixedOrder(clusters)]
        combined_gsea_res$cluster <- factor(combined_gsea_res$cluster,levels=clusters)

        write.table(combined_gsea_res,paste0(outDir,"/",gene.sets.name,"/gsea_results.tsv"),quote=F,sep="\t",row.names=F)

        save(logfc.data,combined_gsea_res,file=paste0(outDir,"/",gene.sets.name,"/gsea_results.Rdata"))


        make_gsea_plots(combined_gsea_res,outDir,gene.sets.name)


}
