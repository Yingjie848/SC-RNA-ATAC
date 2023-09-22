
# secondary analysis for archR results

rm(list=ls())

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

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)

setwd("../run_ArchR/MBC6_excludeC10/")

outdir="secondary_analysis"; dir.create(outdir)


seurat_atac <- readRDS("MBC6_archrData_in_seurat.rds")
seurat_rna <- readRDS("../../run_seurat/MBC6_highCov_atac_20230522/seurat_cellRangerRawData/initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/reclustered_exclude_C10/seurat_rna.rds")


str(seurat_rna[['ATAC']])


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


atac_data <- GetAssayData(seurat_atac, assay = "peaks", slot = "data")

seurat_atac_ <- CreateChromatinAssay(atac_data, fragments = frags.500)


# 
combined <- merge(x=seurat_rna,y=seurat_atac)


# make violon plots
p <- VlnPlot(
  object = seurat_atac,
  features = c(
               'TSSEnrichment', 'NucleosomeRatio','PromoterRatio'),
  pt.size = 0.1,
  ncol = 5
)

gene_name="FN1"
p <- CoveragePlot(combined,assay='ATAC',features=gene_name,region=gene_name,extend.upstream = 50000,extend.downstream = 50000)


