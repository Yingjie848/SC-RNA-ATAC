# do DE analysis

DEanalysis <- function(obj,markers_pos,markers_pos_top50,outDir){

        #########################################################################################################
        # // annotate up-regulated genes for each cluster

        cat("Annotate up-regulated genes for each cluster\n")

        source("seurat_workflow/code/pathway.R")

        annotate_posDEG_to_pathway(markers_pos,outDir)



        #########################################################################################################
        # make gene heatmap for differentially expressed genes

        print("Make heatmap for top 50 up-regulated DEGs")

        source("seurat_workflow/code/make_heatmap.R")

        gene_info <- markers_pos_top50 %>% group_by(gene) %>% arrange(p_val_adj) %>% summarise(cluster=cluster[1])

        make_heatmap(obj=obj,
                 outDir=paste0(outDir,"/heatmap_top50_posDEGs_logcounts"),
                 genes=unique(markers_pos_top50$gene),
                 show_row_names = F,
                 normalize_to_logcounts = TRUE,
                 gene_info=gene_info)

       
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

}

