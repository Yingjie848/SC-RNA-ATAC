
# analyze single cell RNA-seq data using seurat
# only use RNA-seq data for clustering, for RNA-seq + ATAC-seq, use seurat_rna_atac.R

# resources:
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters-1
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
# https://satijalab.org/seurat/articles/visualization_vignette.html
# Mouse Cell Atlas: http://bis.zju.edu.cn/MCA/
# Multimodal reference mapping: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# Integrate scATAC-seq and scRNA-seq using seurat: https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
#                                    using seurat WNN: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# Seurat: Quality control: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

library(Seurat)
library(SeuratData) 
library(sctransform)
library(Signac)
library(data.table)
library(dplyr)
library(reshape2)
library(ggpubr)
library(circlize)
library(stringr)
library(magrittr)
#library(Signac)
library(cowplot)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(DropletUtils)
#library(DoubletFinder)
library(scDblFinder)
# library(diem) # time cost can be very high for DIEM when there are 2000+ cells

# for scds
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)


# set parameters
Args <- commandArgs()
path_raw             <- Args[7]
frag_path            <- Args[8]
data_source          <- Args[9]
outdir               <- Args[10]

dir.create(outdir,recursive = T)


#########################################################################################################
# main
#########################################################################################################

# load raw data to seurat
if(data_source=='cellranger_atac'){
  dRaw = Read10X_h5(path_raw)
  rna_rawcounts <- dRaw[['Gene Expression']]
}else if(data_source=='starsolo'){
  dRaw = ReadSTARsolo(path_raw)
  rna_rawcounts <- dRaw
}
objRaw <- CreateSeuratObject(counts = rna_rawcounts)

cat("Raw data: ",ncol(rna_rawcounts)," barcodes loaded\n")

# add percent.mt for mitochonrial RNA
# add percent.ribo for ribosomal RNA
# add percent.hb for hemoglobin RNA, related to red blood cells
# add MALAT1, which is intranuclear localized
objRaw[["percent.mt"]]   <- PercentageFeatureSet(objRaw, pattern = "^MT-")
objRaw[["percent.ribo"]] <- PercentageFeatureSet(objRaw, pattern = "^RP[SL]")
objRaw[["percent.hb"]]   <- PercentageFeatureSet(objRaw, pattern = "^HB[^(P)]")
objRaw[["MALAT1"]]       <- PercentageFeatureSet(objRaw,pattern="MALAT1")


# // create seurat object for ATAC data
atac_counts = dRaw$Peaks

# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts   <- atac_counts[as.vector(grange.use), ]

# ATAC analysis add gene annotation information
annotations                 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC" # solution when it gets error for chrom, update GenomeInfoDb: https://github.com/Bioconductor/GenomeInfoDb/issues/82#issuecomment-1409482932
genome(annotations)         <- "hg38"

chrom_assay   <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"),min.cells = 0, min.features = 0,genome='hg38',annotation=annotations,fragments=frag_path)
objRaw[['ATAC']] <- chrom_assay

#########################################################################################################
# emptyDrops
cat("Running emptyDrops\n")

# // check total counts with rank
bcrank <- barcodeRanks(rna_rawcounts)
pdf(paste0(outdir,"/emptyDrops_totalCounts_rank.pdf"))
plot(bcrank$rank, bcrank$total, log="xy",xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
dev.off()

# // run emptyDrops
# one has to use raw counts before filtering: 
# ?emptyDrops and 
# http://bioconductor.org/books/3.15/OSCA.advanced/droplet-processing.html#qc-droplets
set.seed(100)
if(file.exists(paste0(outdir,"/emptyDrops_output.tsv"))){
    print("Use existing emptyDrops output")
    emptyDrop_output <- read.table(paste0(outdir,"/emptyDrops_output.tsv"),header=T)
}else{
    print("Running emptyDrops")
    emptyDrop_output <- emptyDrops(rna_rawcounts,lower=100,niters=20000) 
    write.table(emptyDrop_output,paste0(outdir,"/emptyDrops_output.tsv"),quote=F,sep="\t",row.names=T)
}
is.cell = emptyDrop_output$FDR < 0.001
summary(is.cell)


# Check if p-values are lower-bounded by 'niters' (increase 'niters' if any Limited==TRUE and Sig==FALSE)
table(Sig=is.cell, Limited=emptyDrop_output$Limited)
emptyDrops_cells <- rownames(emptyDrop_output)[which(emptyDrop_output$FDR < 0.001)]
cat("Raw data passed emptyDrops ",length(emptyDrops_cells),"\n")

# // raw data

objRaw@meta.data$emptyDrops <- ifelse(rownames(objRaw@meta.data) %in% emptyDrops_cells,'No','Yes')
write.table(objRaw@meta.data,paste0(outdir,"/emptyDrops_rawData_metadata.txt"),quote=F,sep="\t")

# check how many cells in raw barcodes passed emptyDrops
cat("How many cells passed emptyDrops in raw data\n")
table(objRaw@meta.data$emptyDrops)

# compare cells and empty drops
p <- ggboxplot(objRaw@meta.data,x='emptyDrops',y='nCount_RNA',add='jitter',xlab='Empty Droplets') + scale_y_log10()
pdf(paste0(outdir,"/emptyDrops_nCount.rawData.pdf"));print(p);dev.off()


#########################################################################################################
# filter out empty droplets before detecting doublets

obj <- subset(x = objRaw, subset = emptyDrops=="No")
cat("Raw data passed emptyDrops ",nrow(obj@meta.data),"\n")


#########################################################################################################
# detect doublets using scDblFinder

set.seed(123)

if(file.exists(paste0(outdir,"/scDblFinder.rds"))){
  print("Use existing scDblFinder output")
  dbl = readRDS(paste0(outdir,"/scDblFinder.rds"))
}else{
  print("Running scDblFinder")
  dbl <- scDblFinder(obj$RNA@counts)
  saveRDS(dbl,file=paste0(outdir,"/scDblFinder.rds"))
}

cat("scDblFinder results:")
table(dbl$scDblFinder.class)

# add doublets to meta.data
doublets <- data.frame(cell=colnames(dbl),scDblFinder.doublets=dbl$scDblFinder.class,scDblFinder.score=dbl$scDblFinder.score)

obj@meta.data <- cbind(obj@meta.data,
                       doublets[match(rownames(obj@meta.data),doublets$cell),c('scDblFinder.doublets','scDblFinder.score')])
write.table(obj@meta.data,paste0(outdir,"/scDblFinder_metadata.txt"),quote=F,sep="\t")


# check doublets score
p <- ggboxplot(obj@meta.data,x="scDblFinder.doublets",y="scDblFinder.score")
pdf(paste0(outdir,"/scDblFinder_score_boxplot.pdf"));print(p);dev.off()

p <- gghistogram(obj@meta.data,x="scDblFinder.score",fill='steelblue',title="scDblFinder score",ylab="Number of cells")
pdf(paste0(outdir,"/scDblFinder_score.pdf"));print(p);dev.off()

# compare doublet and singlet
p <- ggboxplot(obj@meta.data,x='scDblFinder.doublets',y='nCount_RNA',add='jitter',xlab='') + scale_y_log10()
pdf(paste0(outdir,"/scDblFinder_nCount.pdf"));print(p);dev.off()

p <- ggscatter(obj@meta.data,x='scDblFinder.score',y='nCount_RNA')
pdf(paste0(outdir,"/scDblFinder_nCount.scatterplot.pdf"));print(p);dev.off()


#########################################################################################################
# find doublets using scds
set.seed(30519)

if(file.exists(paste0(outdir,"/scds.rds"))){
  sce = readRDS(paste0(outdir,"/scds.rds"))
}else{
  sce = as.SingleCellExperiment(obj)
  sce = cxds(sce,retRes = TRUE)
  sce = bcds(sce,retRes = TRUE,verb=TRUE)
  sce = cxds_bcds_hybrid(sce)

  saveRDS(sce,file=paste0(outdir,"/scds.rds"))
}

# add scds to meta data
obj@meta.data <- cbind(obj@meta.data,data.frame(cxds_score=sce$cxds_score,bcds_score=sce$bcds_score,hybrid_score=sce$hybrid_score))

write.table(obj@meta.data,paste0(outdir,"/scds_metaData.txt"),quote=F,sep="\t")

# make doublet score vs nCount
p1 <- ggscatter(obj@meta.data,x='hybrid_score',y='nCount_RNA')
p2 <- ggscatter(obj@meta.data,x='cxds_score',y='nCount_RNA')
p3 <- ggscatter(obj@meta.data,x='bcds_score',y='nCount_RNA')

pdf(paste0(outdir,"/scds_nCount.pdf"),12,4);
print(ggarrange(p1,p2,p3,ncol=3))
dev.off()

#########################################################################################################
# find doublets using scrublet

# first write out RNA count matrix
counts <- obj$RNA@counts
write.table(t(counts),paste0(outdir,"/RNA_counts.txt"),quote=F,sep="\t",row.names=F,col.names=F)

# run scrublet
# 1) open another terminal
# 2) load environment: source /home/zhuy1/my_apps/miniconda3/bin/activate /home/zhuy1/my_apps/miniconda3
# 3) start python console: python
# 4) run_scrublet.py


# load 
scrublet_doublets <- read.table("scrublet_doublets.txt")
scrublet_score <- read.table("scrublet_score.txt")

# 
obj@meta.data$scrublet_doublets <- scrublet_doublets$V1
obj$scrublet_score <- scrublet_score$V1

# 
hist(obj$scrublet_score)


#########################################################################################################

saveRDS(obj,file=paste0(outdir,"/seurat_rna.rds"))


#########################################################################################################
# do clustering and find clusters

source("seurat_workflow/code/clustering_analysis.R")

clustering_analysis(obj,outDir=paste0(outdir,'/initial_clustering'),npc=30,resolution = 0.5, secondary_analysis=TRUE,stop_after_marker_genes = TRUE,regress_out_cellcycle=FALSE,clustering_assay=c('RNA') )

#clustering_analysis(obj,outDir=paste0(outdir,'/initial_clustering'),npc=30,resolution = 0.5, stop_after_marker_genes = TRUE,regress_out_cellcycle=FALSE,clustering_assay=c('RNA','ATAC') )


#########################################################################################################
# compare scDblFinder and scrublet
obj <- readRDS(paste0(outdir,'/initial_clustering/seurat_rna.rds'))

out <- table(obj$scDblFinder.doublets,obj$scrublet_doublets)
write.table(out,paste0(outdir,'/initial_clustering/compare_scDblDinder_scrublet.table.txt'))

# make DimPlot by Doublets

p_scDbl    <- DimPlot(obj, reduction = "umap", group.by = "scDblFinder.doublets", label = FALSE, label.size = 3.5, repel = TRUE,order=T) + ggtitle("scDblFinder doublets")
p_scrublet <- DimPlot(obj, reduction = "umap", group.by = "scrublet_doublets", label = FALSE, label.size = 3.5, repel = TRUE,order=T) + ggtitle("scrublet doublets")

pdf(paste0(outdir,'/initial_clustering/compare_scDblDinder_scrublet.DimPlot.pdf'),10,5); print(p_scDbl + p_scrublet); dev.off()

# make venn diagram


#########################################################################################################
# check MALAT1
plot(obj$MALAT1,obj$nCount_RNA)
plot(obj$MALAT1,obj$percent.mt)