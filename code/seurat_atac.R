# atac only analysis
# https://stuartlab.org/signac/articles/pbmc_vignette.html
# https://stuartlab.org/signac/articles/mouse_brain_vignette.html
# https://stuartlab.org/signac/articles/pbmc_multiomic.html

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)


library(Seurat)
library(SeuratData) 
library(sctransform)
library(Signac)
library(data.table)
library(dplyr)
library(ggpubr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)


outdir="initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/atac_analysis"; dir.create(outdir)

# load RNA+ATAC seurat object, including seurat_clusters
obj <- readRDS("initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/seurat_rna.rds")

# check DimPlot
p <- DimPlot(obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("UMAP (RNA + ATAC)")
pdf(paste0(outdir,"/DimPlot.pdf")); print(p);dev.off()


#############################################################################
# examine QC metrics
# set default assay to ATAC
DefaultAssay(obj) <- 'ATAC'

# compute nucleosome signal score per cell
obj <- NucleosomeSignal(obj)

# compute TSS enrichment score per cell
obj <- TSSEnrichment(object = obj, fast = FALSE)

obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')

p1 <- TSSPlot(obj)
p2 <- TSSPlot(obj, group.by = 'high.tss') + NoLegend()

# We can also look at the fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength. You can see that cells that are outliers for the mononucleosomal / nucleosome-free ratio (based on the plots above) have different nucleosomal banding patterns. The remaining cells exhibit a pattern that is typical for a successful ATAC-seq experiment.
obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p3 <- FragmentHistogram(object = obj, group.by = 'nucleosome_group')

# make violon plots
p4 <- VlnPlot(
  object = obj,
  features = c(
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


pdf(paste0(outdir,"/TSSplot_by_clusters.pdf")); print(p1);dev.off()
pdf(paste0(outdir,"/TSSplot_by_tss_enrichment.pdf")); print(p2);dev.off()
pdf(paste0(outdir,"/FragmentHist.pdf")); print(p3);dev.off()
pdf(paste0(outdir,"/qc_vlnPlot.pdf")); print(p4);dev.off()


#############################################################################
# identify DA peaks for clusters
da_peaks <- FindAllMarkers(
  object = obj,
  test.use = 'LR'
)

da_peaks <- cbind(data.frame(peak=rownames(da_peaks)),da_peaks)

# find closest genes for each peak
closest_genes <- unique(ClosestFeature(obj, regions = da_peaks$peak))

closest_genes <- closest_genes %>% filter(distance < 1000)

da_peaks <- merge(da_peaks,closest_genes,by.x='peak',by.y='query_region',all.x=T) %>% arrange(p_val)

da_peaks_pos <- da_peaks %>% filter(avg_log2FC > 1 & p_val_adj < 0.001)
da_peaks_neg <- da_peaks %>% filter(avg_log2FC < -1 & p_val_adj < 0.001)

dir.create(paste0(outdir,"/DA_analysis"))
write.table(da_peaks,paste0(outdir,"/DA_analysis/DA_peaks.txt"),quote=F,sep="\t",row.names=F)
write.table(da_peaks_pos,paste0(outdir,"/DA_analysis/DA_peaks_pos.txt"),quote=F,sep="\t",row.names=F)
write.table(da_peaks_neg,paste0(outdir,"/DA_analysis/DA_peaks_neg.txt"),quote=F,sep="\t",row.names=F)

# make heatmap for positive peaks
source("seurat_workflow/code/make_heatmap_atac.R")

make_heatmap_atac(obj=obj,
                 outDir=paste0(outdir,"/DA_analysis/da_peaks_pos"),
                 peaks=da_peaks_pos[!is.na(da_peaks_pos$gene_name),]$peak,
                 peak_info = da_peaks_pos,
                 show_row_names = F,
                 cluster_rows = T
    )

# make heatmap for top 20 peaks of each cluster
da_peaks_pos_top20 <- da_peaks_pos %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(desc(abs(avg_log2FC))) %>%
          dplyr::filter(p_val_adj < 0.001) %>%
          dplyr::slice_head(n = 20)

make_heatmap_atac(obj=obj,
                 outDir=paste0(outdir,"/DA_analysis/da_peaks_pos_top20"),
                 peaks=da_peaks_pos_top20[!is.na(da_peaks_pos_top20$gene_name),]$peak,
                 peak_info = da_peaks_pos_top20,
                 show_row_names = T,
                 cluster_rows = T
    )

# make heatmap for top 50 peaks of each cluster
da_peaks_pos_top50 <- da_peaks_pos %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(desc(abs(avg_log2FC))) %>%
          dplyr::filter(p_val_adj < 0.001) %>%
          dplyr::slice_head(n = 50)

make_heatmap_atac(obj=obj,
                 outDir=paste0(outdir,"/DA_analysis/da_peaks_pos_top50"),
                 peaks=da_peaks_pos_top50[!is.na(da_peaks_pos_top50$gene_name),]$peak,
                 peak_info = da_peaks_pos_top50,
                 show_row_names = T,
                 cluster_rows = T
    )




#############################################################################
# motif analysis

# Adding motif information to the Seurat object
# To add the DNA sequence motif information required for motif analyses, we can run the AddMotifs() function:
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# get top differentially accessible peaks
top.da.peak <- da_peaks_pos %>% filter(avg_log2FC > 1 & p_val_adj < 0.001)

# test enrichment
clusters <- as.vector(unique(top.da.peak$cluster))
enriched.motifs.allClusters <- lapply(clusters,function(x){
    top.da.peak.cluster <- top.da.peak %>% filter(cluster==x)
    print(nrow(top.da.peak.cluster))
    enriched.motifs <- FindMotifs(
                            object = obj,
                            features = top.da.peak.cluster$peak)
    enriched.motifs$cluster <- x
    enriched.motifs
})
enriched.motifs <- do.call(rbind,enriched.motifs.allClusters) %>% filter(p.adjust<0.05)

write.table(top.da.peak,paste0(outdir,"/motif_analysis/top_da_peak.txt"),quote=F,sep="\t",row.names=F)
write.table(enriched.motifs,paste0(outdir,"/motif_analysis/enriched.motifs.txt"),quote=F,sep="\t",row.names=F)

# We can also plot the position weight matrices for the motifs, so we can visualize the different motif sequences.
clusters <- as.vector(unique(enriched.motifs$cluster))
cache <- lapply(clusters,function(x){
    
    enriched.motifs.cluster <- enriched.motifs %>% filter(cluster==x)
    
    plot_motif_top <- MotifPlot(
      object = obj,
      motifs = enriched.motifs.cluster$motif[1:12]
    )

    dir.create(paste0(outdir,"/motif_analysis"))
    pdf(paste0(outdir,"/motif_analysis/motifPlot.cluster_",x,".top1-12.pdf")); print(plot_motif_top);dev.off()

    plot_motif_top <- MotifPlot(
      object = obj,
      motifs = enriched.motifs.cluster$motif[13:24]
    )

    pdf(paste0(outdir,"/motif_analysis/motifPlot.cluster_",x,".top13-24.pdf")); print(plot_motif_top);dev.off()
})




#############################################################################
# save

saveRDS(obj,file=paste0(outdir,"/seurat_atac.rds"))

#############################################################################
# // for chondroid DEGs, what about their ATAC signal

DefaultAssay(obj) <- 'ATAC'

# load DE genes
chondroid_genes <- data.frame(gene=c('FN1','SMOC2','CHSY1','CHST9','CSGALNACT1','CHST11','CHSY1','CREB3L2','XYLT1','SDC2','RARB','MEF2C','SOX5','SOX6','SOX9'))

# find closest genes for each peak
peaks <- rownames(obj$ATAC@data)
closest_genes <- unique(ClosestFeature(obj, regions = peaks))
closest_genes <- closest_genes %>% filter(distance < 1000)

genes_and_peaks <- merge(chondroid_genes,closest_genes,by.x='gene',by.y='gene_name',all.x=T)

# visualize gene ATAC signal
dir.create(paste0(outdir,"/coveragePlots_chondroid_genes"));
genes_and_peaks.overlapped <- genes_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(genes_and_peaks.overlapped$gene)){
  print(gene_name)
  
  p <- CoveragePlot(obj,features=gene_name,region=gene_name,extend.upstream = 50000,extend.downstream = 50000)
  pdf(paste0(outdir,"/coveragePlots_chondroid_genes/coveragePlot.",gene_name,".pdf"),7,5); print(p);dev.off()

}


# make rna and atac heatmap
source("seurat_workflow/code/make_heatmap_rna_and_atac.R")

# visualize gene ATAC signal
dir.create(paste0(outdir,"/heatmap_rna_and_atac_chondroid_genes"));
genes_and_peaks.overlapped <- genes_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(genes_and_peaks.overlapped$gene)){
  print(gene_name)
  t <- genes_and_peaks.overlapped %>% filter(gene==gene_name)
  p <- make_heatmap_rna_and_atac(obj,gene_name=gene_name,peaks=t$query_region)
  pdf(paste0(outdir,"/heatmap_rna_and_atac_chondroid_genes/heatmap_rna_and_atac.",gene_name,".pdf"),15,5); print(p);dev.off()

}



#############################################################################
# // for up-regulated genes, what about their ATAC signal

DefaultAssay(obj) <- 'ATAC'

# load DE genes
DEGs <- read.table("initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/DE_analysis/posDE_genes.txt",header=T)

# find closest genes for each peak
peaks <- rownames(obj$ATAC@data)
closest_genes <- unique(ClosestFeature(obj, regions = peaks))
closest_genes <- closest_genes %>% filter(distance < 1000)

DEGs_and_peaks <- merge(DEGs,closest_genes,by.x='gene',by.y='gene_name',all.x=T) %>% arrange(p_val)

# visualize gene ATAC signal
dir.create(paste0(outdir,"/coveragePlots_posDE_genes_50kUP_20kDOWN"));
DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  
  p <- CoveragePlot(obj,features=gene_name,region=gene_name,extend.upstream = 50000,extend.downstream = 20000)
  pdf(paste0(outdir,"/coveragePlots_posDE_genes_50kUP_20kDOWN/coveragePlot.",gene_name,".pdf"),7,5); print(p);dev.off()

}


# make rna and atac heatmap
source("seurat_workflow/code/make_heatmap_rna_and_atac.R")

# visualize gene ATAC signal
dir.create(paste0(outdir,"/heatmap_rna_and_atac_posDE_genes"));
DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  t <- DEGs_and_peaks.overlapped %>% filter(gene==gene_name)
  p <- make_heatmap_rna_and_atac(obj,gene_name=gene_name,peaks=t$query_region)
  pdf(paste0(outdir,"/heatmap_rna_and_atac_posDE_genes/heatmap_rna_and_atac.",gene_name,".pdf"),15,5); print(p);dev.off()

}

#############################################################################
# // for down-regulated genes, what about their ATAC signal

DefaultAssay(obj) <- 'ATAC'

# load DE genes
DEGs <- read.table("initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/DE_analysis/DE_genes.txt",header=T); nrow(DEGs)
DEGs <- DEGs %>% filter(avg_log2FC<0 & p_val_adj < 0.001); nrow(DEGs)


# find closest genes for each peak
peaks <- rownames(obj$ATAC@data)
closest_genes <- unique(ClosestFeature(obj, regions = peaks))
closest_genes <- closest_genes %>% filter(distance < 1000)

DEGs_and_peaks <- merge(DEGs,closest_genes,by.x='gene',by.y='gene_name',all.x=T) %>% arrange(p_val)

# visualize gene ATAC signal
dir.create(paste0(outdir,"/coveragePlots_negDE_genes"));
DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  
  p <- CoveragePlot(obj,features=gene_name,region=gene_name,extend.upstream = 10000,extend.downstream = 10000)
  pdf(paste0(outdir,"/coveragePlots_negDE_genes/coveragePlot.",gene_name,".pdf"),7,5); print(p);dev.off()

}


# make rna and atac heatmap
source("seurat_workflow/code/make_heatmap_rna_and_atac.R")

# visualize gene ATAC signal
dir.create(paste0(outdir,"/heatmap_rna_and_atac_negDE_genes"));
DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  t <- DEGs_and_peaks.overlapped %>% filter(gene==gene_name)
  p <- make_heatmap_rna_and_atac(obj,gene_name=gene_name,peaks=t$query_region)
  pdf(paste0(outdir,"/heatmap_rna_and_atac_negDE_genes/heatmap_rna_and_atac.",gene_name,".pdf"),15,5); print(p);dev.off()

}

#############################################################################
# intersect DE genes and DA peaks, find DE genes having any peak differentially accessible 

DefaultAssay(obj) <- 'ATAC'

DE_genes <- read.table("initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/DE_analysis/posDE_genes.txt",header=T); nrow(DE_genes)

DA_peaks <- read.table("initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30/atac_analysis/DA_analysis/DA_peaks.txt",header=T)

DE_genes <- DE_genes %>% filter(gene %in% DA_peaks$gene_name); nrow(DE_genes)

# find closest genes for each peak
peaks <- rownames(obj$ATAC@data)
closest_genes <- unique(ClosestFeature(obj, regions = peaks))
closest_genes <- closest_genes %>% filter(distance < 1000)

DEGs_and_peaks <- merge(DE_genes,closest_genes,by.x='gene',by.y='gene_name',all.x=T) %>% arrange(p_val)


# visualize gene ATAC signal
dir.create(paste0(outdir,"/coveragePlots_posDE_genes_intersected_with_DA_peaks"));

write.table(DEGs_and_peaks,paste0(outdir,"/coveragePlots_posDE_genes_intersected_with_DA_peaks/posDE_genes_intersected_with_DA_peaks.tsv"),quote=F,sep="\t",row.names=F)

DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  
  p <- CoveragePlot(obj,features=gene_name,region=gene_name,extend.upstream = 10000,extend.downstream = 10000)
  pdf(paste0(outdir,"/coveragePlots_posDE_genes_intersected_with_DA_peaks/coveragePlot.",gene_name,".pdf"),7,5); print(p);dev.off()

}


# make rna and atac heatmap
source("seurat_workflow/code/make_heatmap_rna_and_atac.R")

# visualize gene ATAC signal
dir.create(paste0(outdir,"/heatmap_rna_and_atac_posDE_genes_intersected_with_DA_peaks"));

write.table(DEGs_and_peaks,paste0(outdir,"/heatmap_rna_and_atac_posDE_genes_intersected_with_DA_peaks/posDE_genes_intersected_with_DA_peaks.tsv"),quote=F,sep="\t",row.names=F)

DEGs_and_peaks.overlapped <- DEGs_and_peaks %>% filter(!is.na(query_region))
for(gene_name in unique(DEGs_and_peaks.overlapped$gene)){
  print(gene_name)
  t <- DEGs_and_peaks.overlapped %>% 
          filter(gene==gene_name) %>% 
          mutate(peak_chr=as.numeric(gsub("chr","",sapply(query_region,function(x) strsplit(x,"-")[[1]][1]))),
                 peak_start=as.numeric(sapply(query_region,function(x) strsplit(x,"-")[[1]][2])),
                 peak_end=as.numeric(sapply(query_region,function(x) strsplit(x,"-")[[1]][3]))
          ) %>% arrange(peak_chr,peak_start)
  p <- make_heatmap_rna_and_atac(obj,gene_name=gene_name,peaks=t$query_region)
  pdf(paste0(outdir,"/heatmap_rna_and_atac_posDE_genes_intersected_with_DA_peaks/heatmap_rna_and_atac.",gene_name,".pdf"),15,5); print(p);dev.off()

}



stop()


plot1 <- VlnPlot(
  object = obj,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("0","1","2")
)

plot2 <- FeaturePlot(
  object = obj,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

open_C0 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_others <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_ <- ClosestFeature(obj, regions = da_peaks_pos$peak)
closest_genes_others <- ClosestFeature(obj, regions = open_others)

CoveragePlot(
  object = obj,
  region = 'chr1-159929891-159929975',features='IGSF9',
  extend.upstream = 40000,
  extend.downstream = 20000
)