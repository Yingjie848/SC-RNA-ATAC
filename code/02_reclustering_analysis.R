
# perform subclusering for cells

source("seurat_workflow/code/clustering_analysis.R")
source("seurat_workflow/code/rename_clusters.R")

# main
main <- function(){


    #################################################################
    # // MBC6_highCov_atac_20230522
    # exclude low quality cells, npc30, res0p8
    obj            <- readRDS("initial_clustering/seurat_rna.rds"); dim(obj)
    obj            <- subset(obj,subset=nCount_RNA>500 & MALAT1 > 1 & percent.mt < 3 & (scDblFinder.doublets != 'doublet' | scrublet_doublets!=1)); dim(obj)
    outdir         <- "initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/" 
    clustering_analysis(obj,outDir=outdir,npc=30,resolution=0.8,clustering_assay=c('RNA'),secondary_analysis=FALSE)
    obj            <- readRDS(paste0(outdir,"/seurat_rna.rds"))
    secondary_analysis(obj,outDir=outdir,stop_after_marker_genes = FALSE, stop_after_DEA=TRUE,clustering_assay=c('RNA'), evaluate_stemness=FALSE)

    # exclude C10, it is cell cycle related, and not found in ATAC data
    obj            <- readRDS("initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/seurat_rna.rds"); dim(obj)
    obj            <- subset(obj,subset= seurat_clusters != 10); dim(obj)
    outdir         <- "initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/reclustered_exclude_C10" 
    clustering_analysis(obj,outDir=outdir,npc=30,resolution=0.8,clustering_assay=c('RNA'),secondary_analysis=FALSE)
    obj            <- readRDS(paste0(outdir,"/seurat_rna.rds"))
    secondary_analysis(obj,outDir=outdir,stop_after_marker_genes = FALSE, stop_after_DEA=TRUE,clustering_assay=c('RNA'), evaluate_stemness=FALSE)

    # annotate clusters
    named_clusters <- c(
        "0"="Epithelial cells",
        "1"="Epithelial cells",
        "2"="Epithelial cells",
        "3"="Epithelial cells",
        "4"="Epithelial cells",
        "5"="Macrophages",
        "6"="Epithelial cells",
        "7"="T cells",
        "8"="Fibroblast cells",
        "9"="Endothelial cells",
        "10"="Pericytes"
    )

    outDir         <- "initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/reclustered_exclude_C10/annotated_clusters" 

    obj <- add_cell_identity(obj,named_clusters=named_clusters,outDir)

}
