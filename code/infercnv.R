# Infer CNV using infercnv based on single cell RNA-seq data

library(infercnv)
library(Seurat)
library(dplyr)

# setwd("/home/zhuy1/my_projects/metaplastic_scRNAseq/run_seurat_MBC1")

Args <- commandArgs()
seurat_obj      <- Args[6]
outdir          <- Args[7]
ref_groups      <- Args[8]
test_groups     <- Args[9]
gene_order_file <- Args[10] # scripts/data/GRCh38_protein_coding_genes.processed.txt, exclude chrM, ordered by chr coordinate, exclude duplicate genes
cell_annot_file <- Args[11]

test_groups <- strsplit(test_groups,',')[[1]]
ref_groups  <- strsplit(ref_groups,',')[[1]]

dir.create(outdir)
dir.create(paste0(outdir,"/results"))

# create RNA counts matrix and cell annotation files
print("prepare RNA counts matrix file")
obj    <- readRDS(seurat_obj)
counts <- obj$RNA@counts

print("prepare cell annotations file")
if(file.exists(cell_annot_file)){
    cell_annot = read.table(cell_annot_file)
    colnames(cell_annot) <- c('cell','cluster')
}else{
    cell_annot=data.frame(cell=rownames(obj@meta.data),
                         cluster=as.vector(obj@meta.data$seurat_clusters) )
}

# extract test cells
if('all' %in% test_groups){
    test_groups=as.vector(unique(cell_annot$cluster))
}

cell_annot <- cell_annot %>% 
                filter(cluster %in% c(ref_groups,test_groups))

counts <- counts[,cell_annot$cell]

gzfile <- gzfile(paste0(outdir,"/counts.txt.gz"),'w')
write.table(as.matrix(counts),gzfile,quote=F,sep="\t")
close(gzfile)

write.table(cell_annot,paste0(outdir,"/cell_annotations.txt"),quote=F,sep="\t",col.names=F,row.names=F)


# Create the InferCNV Object
print("Create the InferCNV Object")
infercnv_obj = CreateInfercnvObject(
                    raw_counts_matrix=paste0(outdir,"/counts.txt.gz"),
                    annotations_file=paste0(outdir,"/cell_annotations.txt"),
                    delim="\t",
                    gene_order_file=gene_order_file, # "~/my_software/infercnv/extdata/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
                    ref_group_names=ref_groups
               )

# Run infercnv
# refere: Slyper et al. 2020, Nature Medicine
print("Run infercnv")

infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=paste0(outdir,"/results"),
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=FALSE,
    no_prelim_plot=FALSE,
    save_rds=FALSE,
    output_format="pdf",
    useRaster=TRUE
)

print("Run infercnv with HMM")
infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=paste0(outdir,"/results"),
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=TRUE,
    no_prelim_plot=FALSE,
    save_rds=FALSE,
    output_format="pdf",
    useRaster=TRUE
)

