# evaluate stemness using stemness index (average expression values log2(counts/libSize*10000+1))
library(GSEABase)
library(GSVA)
library(dplyr)
library(ggpubr)

evaluate_stemness <- function(obj,outDir=NULL,genes=NULL,gene.sets.name=NULL){

    source("seurat_workflow/code/calc_sigscore.R")
    source("seurat_workflow/code/make_heatmap.R")

    dir.create(outDir,recursive = T)

    # using 35 stem cell markers collected by Wu et al. ONCOLOGY REPORTS 2019 Table SI
    #gmt.file   <- "seurat_workflow/data/marker_genes/cell_origin/stemness_index_wu_oncology_reports_2019.gmt" 
    #genesets   <- getGmt(gmt.file)
    #genes <- geneIds(genesets[['stem_cell']])

    # calculate signature score according to average expression values of log2(counts/libSize*10000+1)
    obj@meta.data$stemness_index <- calc_sigscore(obj$RNA@counts,genes)

    # calculate ssGSEA score
    genes.list <- list(genes); names(genes.list) <- gene.sets.name
    gsva_output <- gsva(as.matrix(obj$RNA@data),genes.list,method="ssgsea",ssgsea.norm=TRUE)
    obj@meta.data$stemness_ssgsea <- t(gsva_output)[,1]


    # make violin plot for stemness_index
    cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
    cell_clusters     <- sort(unique(obj@meta.data$seurat_clusters))
    cluster_cols             = cluster_cols[1:length(cell_clusters)]
    names(cluster_cols)      = cell_clusters

    p <- obj@meta.data %>% 
            ggviolin(x='seurat_clusters',y='stemness_index',color='seurat_clusters',title='Stemness index',xlab='Clusters',ylab='Stemness Index',add='boxplot') + 
            scale_color_manual(values=cluster_cols) + 
            theme(legend.position = 'none')
    pdf(paste0(outDir,"/stemness_index.pdf"),7,5);print(p);dev.off()

    p <- obj@meta.data %>% 
            ggviolin(x='seurat_clusters',y='stemness_index',color='seurat_clusters',title='Stemness index',xlab='Clusters',ylab='Stemness Index',add='boxplot') + 
            scale_color_manual(values=cluster_cols) + 
            theme(legend.position = 'none') +
            stat_compare_means(ref.group = '0',label='p.signif') +
            geom_hline(yintercept = mean(obj@meta.data$stemness_index,na.rm=T),color='grey')
    pdf(paste0(outDir,"/stemness_index.add_signif.pdf"),7,5);print(p);dev.off()

    # make violin plot for stemness_ssgsea
    p <- obj@meta.data %>% 
            ggviolin(x='seurat_clusters',y='stemness_ssgsea',color='seurat_clusters',title='Stemness ssGSEA',xlab='Clusters',,ylab='Stemness ssGSEA',add='boxplot') + 
            scale_color_manual(values=cluster_cols) + 
            theme(legend.position = 'none')
    pdf(paste0(outDir,"/stemness_ssgsea.pdf"),7,5);print(p);dev.off()

    p <- obj@meta.data %>% 
            ggviolin(x='seurat_clusters',y='stemness_ssgsea',color='seurat_clusters',title='Stemness ssGSEA',xlab='Clusters',ylab='Stemness ssGSEA',add='boxplot') + 
            scale_color_manual(values=cluster_cols) + 
            theme(legend.position = 'none') +
            stat_compare_means(ref.group = '0',label='p.signif') +
            geom_hline(yintercept = mean(obj@meta.data$stemness_ssgsea,na.rm=T),color='grey')
    pdf(paste0(outDir,"/stemness_ssgsea.add_signif.pdf"),7,5);print(p);dev.off()


    # make gene heatmap
    make_heatmap(obj,outDir=paste0(outDir,'/gene_heatmap'),genes=genes,show_row_names=T,normalize_to_logcounts=TRUE,cluster_rows=T)
    
    # out table
    write.table(obj@meta.data,paste0(outDir,"/metaData_stemness.txt"),quote=F,sep="\t",row.names=F)

    obj

}

