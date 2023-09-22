# compare CNV and RNA expression to see if differential expression of some genes are caused by CNV rather than cell differentiation/development 

rm(list=ls())

library(circlize)
library(ComplexHeatmap)
library(dplyr)

dataDir <- "MBC6_highCov_atac/seurat_cellRangerRawData/initial_clustering/reclustered_nCount500-30000_MALAT1-1pct_mt3pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_1_3_res0p2_npc30"
dataDir <- "MBC218_highCov_atac/seurat_cellRangerRawData/initial_clustering/reclustered_nCountRNA500-30000_nCountATAC100-30000_MALAT1-1pct_mt10pct_scDblFinder0p5_npc30_res0p3/reclustered_tumor_C0_3_6_9_12_res0p3_npc30/reclustered_tumor_C0_1_3_4_5_6_res0p3_npc30/"

infercnvDir <- paste0(dataDir,'/infercnv/') 
dir.create(paste0(dataDir,"/analysis/associate_cnv_and_rna/"),recursive = T)

################################################################################################################
# loading inputs

# load cnv for ref and obs cells
references <- read.table(paste0(infercnvDir,"/results/infercnv.references.txt"),header=T,check.names = F)
observations <- read.table(paste0(infercnvDir,"/results/infercnv.observations.txt"),header=T,check.names = F)

# load DE results
DEres <- read.table(paste0(dataDir,"/DE_analysis/posDE_genes.txt"),header=T)

# load metaData including seurat_clusters
metaData <- read.table(paste0(dataDir,"/metaData.txt"),header=T)

# load gene annotations
gene_annot <- read.table("seurat_workflow/data/GRCh38_protein_coding_genes.processed.txt") %>% mutate(chr_id=as.numeric(gsub('Y',24,gsub('X',23,gsub('chr','',V2))))) %>% arrange(chr_id,V3)


################################################################################################################
# make cnv and rna heatmap for same genes and cells

# get copy numebr for DEGs
obs_DEGs <- observations[rownames(observations) %in% DEres$gene,]

clusters <- metaData[match(colnames(obs_DEGs),metaData$cell),]$seurat_clusters
names(clusters) <- colnames(obs_DEGs)

# // make cnv heatmap
cnvdata_plot <- t(obs_DEGs)

# get chr, order matrix by chr coordinates
genes <- colnames(cnvdata_plot)
gene_annot.2 <- gene_annot %>% filter(V1 %in% genes)
cnvdata_plot <- cnvdata_plot[,gene_annot.2$V1]
chrs <- gene_annot.2$V2

# set legend
col_fun <- colorRamp2(c(seq(from=min(cnvdata_plot),to=0.8999,length.out=10),seq(from=0.9,to=1.1,length.out=5),seq(from=1.1001,to=max(cnvdata_plot),length.out=10)), colorRampPalette(c("darkblue", "white","darkred"))(25))
legend_params <- list(title = "CNV", at = c(round(min(cnvdata_plot),1), 1,round(max(cnvdata_plot),1)), col_fun= col_fun, 
                        border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), 
                        legend_height = grid::unit(4, "mm"), legend_width = grid::unit(6, "mm"), 
                        grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

# set column annotation

library(RColorBrewer)
chr_cols <- c(brewer.pal(12, "Paired"),brewer.pal(11, "Set3"),'black')
names(chr_cols) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                                          "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                                                          "chr21","chr22","chrX","chrY")
                                
column_ha = HeatmapAnnotation('Chrom'=chrs,col=list('Chrom'=chr_cols),show_annotation_name = FALSE)
h_cnv <- Heatmap(cnvdata_plot,'CNV',show_row_names=F,show_column_names=F,cluster_rows=F,cluster_columns=F,row_split=clusters,row_gap=unit(2, "mm"),col=col_fun,heatmap_legend_param=legend_params,top_annotation = column_ha)

# // make gene expression heatmap
obj <- readRDS(paste0(dataDir,"/seurat_rna.rds"))
rna_data <- as.matrix(obj$RNA@data)
idx_genes <- match(colnames(cnvdata_plot),rownames(rna_data))
idx_cells <- match(rownames(cnvdata_plot),colnames(rna_data))
rnadata_plot <- t(rna_data[idx_genes,idx_cells])

legend_params <- list(title = "RNA",
                        border = "black",title_gp = grid::gpar(fontsize = 8, fontface = "bold"),labels_gp = grid::gpar(fontsize = 8), 
                        legend_height = grid::unit(4, "mm"), legend_width = grid::unit(6, "mm"), 
                        grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

h_rna <- Heatmap(rnadata_plot,name='RNA',show_row_names=F,show_column_names=F,cluster_rows=F,cluster_columns=F,row_split=clusters,row_gap=unit(2, "mm"),heatmap_legend_param=legend_params,top_annotation = column_ha)

dir.create("analysis/associate_cnv_and_rna")
pdf(paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_rna_heatmap.pdf"),20,10)
h_cnv + h_rna
dev.off()

################################################################################################################
# for each gene, compare the correlation between CNV and RNA
genes <- colnames(cnvdata_plot)
corrs <- lapply(1:length(genes),function(i){
    cnv_i <- cnvdata_plot[,i]
    rna_i <- rnadata_plot[,i]
    corr <- cor.test(cnv_i,rna_i)
    data.frame(gene=genes[i],correlation=corr$estimate,p.value=corr$p.value)
})

corrs <- do.call(rbind,corrs)

pdf(paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_rna_correl.hist.pdf"))
hist(corrs$correlation,breaks=50,xlab="Pearson correlation",col="royalblue",main="CNV vs. RNA level correlation")
dev.off()

write.table(corrs,paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_rna_correl.tsv"),quote=F,sep="\t",row.names=F)

################################################################################################################
# check DEGs, do they have higher copy number
# compare one cluster vs. all others for CNV

clusters_distinct <- sort(unique(clusters))

res_cnv_fc <- data.frame()
for(cluster in clusters_distinct){
    test_cells <- clusters==cluster
    cnvdata_plot.test <- cnvdata_plot[test_cells,]
    cnvdata_plot.others <- cnvdata_plot[!test_cells,]

    # for each gene comparing means using wilcox test
    for(i in 1:ncol(cnvdata_plot)){
        gene <- colnames(cnvdata_plot)[i]
        d_test <- cnvdata_plot.test[,i]
        d_others <- cnvdata_plot.others[,i]
        p.value <- wilcox.test(d_test,d_others)$p.value
        mean_1 <- mean(d_test,na.rm=T)
        mean_2 <- mean(d_others,na.rm=T)
        FC <- mean_1/mean_2
        res_cnv_fc <- rbind(res_cnv_fc,data.frame(gene=gene,cell_cluster=cluster,mean_1=mean_1,mean_2=mean_2,FC=FC,logFC=log2(FC),p.value=p.value))
    }
}

res_cnv_fc <- res_cnv_fc %>% group_by(cell_cluster) %>% arrange(desc(FC))
write.table(res_cnv_fc,paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_fc.tsv"),quote=F,sep="\t",row.names=F)

# find genes having high copy number, and exclude them from DEGs, annotate the rest genes
genes_highCNV <- res_cnv_fc %>% filter(FC>=1.1 & p.value < 1e-3)

write.table(genes_highCNV,paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_fc.highCNV.tsv"),quote=F,sep="\t",row.names=F)

DEres_filtered <- DEres %>% filter(!gene %in% genes_highCNV$gene)


# pathway
source("seurat_workflow/code/pathway.R")

dir.create(paste0(dataDir,"/analysis/associate_cnv_and_rna/DEGs_excluded_highCNV/DE_analysis"),recursive = T)

DEres_filtered_top10 <- DEres_filtered %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(desc(abs(avg_log2FC))) %>%
          dplyr::filter(p_val_adj < 0.001) %>%
          dplyr::slice_head(n = 10)

write.table(DEres_filtered,paste0(dataDir,"/analysis/associate_cnv_and_rna/DEGs_excluded_highCNV/DE_analysis/DE_excluded_highCNV.tsv"),quote=F,sep="\t",row.names=F)
write.table(DEres_filtered_top10,paste0(dataDir,"/analysis/associate_cnv_and_rna/DEGs_excluded_highCNV/DE_analysis/DE_excluded_highCNV.top10.tsv"),quote=F,sep="\t",row.names=F)

annotate_posDEG_to_pathway(DEres_filtered,outDir=paste0(dataDir,"/analysis/associate_cnv_and_rna/DEGs_excluded_highCNV"))

# make heatmap for top 10 up-regulated genes
source("seurat_workflow/code/make_heatmap.R")
gene_info <- DEres_filtered_top10 %>% group_by(gene) %>% arrange(p_val_adj) %>% summarise(cluster=cluster[1])
obj <- readRDS(paste0(dataDir,'/seurat_rna.rds'))
make_heatmap(obj=obj,
            outDir=paste0(dataDir,"/analysis/associate_cnv_and_rna/DEGs_excluded_highCNV/DE_analysis/DE_excluded_highCNV.top10.heatmap"),
            genes=unique(DEres_filtered_top10$gene),
            show_row_names = F,
            normalize_to_logcounts = TRUE,
            gene_info=gene_info)


################################################################################################################
# save data
save(cnvdata_plot,rnadata_plot,clusters,DEres,file=paste0(dataDir,"/analysis/associate_cnv_and_rna/cnv_rna.Rdata"))
