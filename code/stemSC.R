
library(StemSC)
library(data.table)

Run_StemSC <- function(obj,outDir="outdir"){

    dir.create(outDir)

    obj.tmp <- obj

    exp <- as.matrix(obj.tmp$RNA@data)

    genes <- rownames(exp)

    hgnc = fread("seurat_workflow/data/hgnc_20220106.txt")
    hgnc = hgnc[,c('Approved symbol','NCBI Gene ID')]
    colnames(hgnc) = c('gene_symbol','ENTREZID')

    entrezID <- sapply(genes,function(x){
        hgnc[match(x,hgnc$gene_symbol),]$ENTREZID
    })

    notNA <- !is.na(entrezID)

    exp.entrezID <- exp[notNA,]
    rownames(exp.entrezID) <- entrezID[notNA]

    print("Calculae StemSC score")
    out <- StemSC(exp.entrezID)

    cell_info <- obj.tmp@meta.data
    cell_info$StemSC <- out
    obj.tmp@meta.data <- cell_info


    
    # make violin plot
    cell_clusters     <- sort(unique(obj.tmp@meta.data$seurat_clusters))

    cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                              "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                              "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                              "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")

    cluster_cols             = cluster_cols[1:length(cell_clusters)]
    names(cluster_cols)      = cell_clusters
    
    p <- obj.tmp@meta.data %>% 
            ggviolin(x='seurat_clusters',y='StemSC',add='boxplot',color='seurat_clusters') + scale_color_manual(values=cluster_cols) + geom_hline(yintercept = mean(obj.tmp@meta.data$StemSC,na.rm=T),col='grey') + theme(legend.position = 'none') + stat_compare_means(ref.group = '0',label='p.signif')

    # out
    pdf(paste0(outDir,'/StemSC_score.pdf'),7,5);print(p);dev.off()
    write.table(obj.tmp@meta.data,paste0(outDir,'/cell_info_StemSC.txt'),quote=F,sep="\t",row.names=F) 
    
    obj.tmp
}

main <- function(){

    out <- Run_StemSC(obj.3,outDir=paste0(outDir,"/StemSC"))
}