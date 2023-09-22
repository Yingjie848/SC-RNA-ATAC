

# load libraries
lapply(c("dplyr","Seurat","HGNChelper","ggpubr"), library, character.only = T)


# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4",
                  "#66bdf0","#7ae5a7","#7adee5","#e57ae0","#dc4a8e","#ddd259","#7140e1")
options(ggplot2.discrete.colour = cluster_cols)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

################################################################################################################################
# annotate by Immune system

# DB file
#db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
db_ = "seurat_workflow/data/ScTypeDB_full_any.xlsx"
tissue = "any" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#
obj <- readRDS("MBC6_highCov_atac_20230522/seurat_cellRangerRawData/initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/reclustered_exclude_C10/seurat_rna.rds")

# get cell-type by cell matrix
RNA_data <- as.matrix(obj[["RNA"]]@data)
es.max = sctype_score(scRNAseqData = RNA_data, scaled = FALSE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# We can also overlay the identified cell types on UMAP plot:
obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  obj@meta.data$customclassif[obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  + ggtitle("") 

outdir="MBC6_highCov_atac_20230522/seurat_cellRangerRawData/initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/reclustered_exclude_C10/annotated_cells_scTypeDB"; dir.create(outdir)

pdf(paste0(outdir,"/DimPlot_umap_scTypeDB.group_by_annotated_cells.pdf"),9,7); print(p); dev.off()

# save results
write.csv(cL_resutls,paste0(outdir,"/cL_results.csv"),row.names=F)
write.csv(sctype_scores,paste0(outdir,"/sctype_scores.csv"),row.names=F)
saveRDS(obj,file=paste0(outdir,"/seurat.rds"))