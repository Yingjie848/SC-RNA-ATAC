
# Quick tutorial: https://www.archrproject.com/articles/Articles/tutorial.html
# ArchRtoSignac : an Object Conversion Package for ArchR to Signac: https://github.com/swaruplabUCI/ArchRtoSignac
library(ArchR)
set.seed(1)
addArchRThreads(threads = 8)

outdir="../run_ArchR/MBC6/"
setwd(outdir)

inputFiles="../../run_cellRanger/MBC6/MBC6/outs/atac_fragments.tsv.gz"
names(inputFiles) <- "MBC6"

addArchRGenome("hg38")

# load seurat results 
seur = readRDS("../../run_seurat/MBC6_highCov_atac_20230522/seurat_cellRangerRawData/initial_clustering/reclustered_nCount500_MALAT1-1pct_mt3pct_scd-scr-olp_npc30_res0p8/seurat_rna.rds")

# Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  validBarcodes=rownames(seur@meta.data)
)

# Creating an ArchRProject 
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "MBC6",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

################################################################################################
# Dimensionality Reduction 

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.8), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=TRUE
)

# We can visualize our scATAC-seq data using a 2-dimensional representation such as Uniform Manifold Approximation and Projection (UMAP). To do this, we add a UMAP embedding to our ArchRProject object with the addUMAP() function.
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE)

# find clusters
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

pdf("umap_archr_clusters.pdf",6,5);print(p);dev.off()


# make marker gene UMAP
p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = "PTPRC", embedding = "UMAP",imputeWeights=NULL,quantCut = c(0.01, 0.95),)


################################################################################################
# integrate scATAC-seq and scRNA-seq, map scRNA clusters to scATAC embeddings
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seur,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    sampleCellsATAC = 10000,
    sampleCellsRNA = 10000,
    corCutOff = 0.75,
    nGenes = 2000
)

# 
getCellColData(proj)
getAvailableMatrices(proj)

# plot embeddings
p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")


# Now that we are satisfied with the results of our scATAC-seq and scRNA-seq integration, we can re-run the integration with addToArrow = TRUE to add the linked gene expression data to each of the Arrow files. 

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seur,
    addToArrow = TRUE,
    groupRNA = "seurat_clusters",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
)

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

pdf("umap_rna_clusters.pdf",6,5);print(p);dev.off()

#markerGenes  <- c(
#    #"CD34", #Early Progenitor
#    #"GATA1", #Erythroid
#    #"PAX5", "MS4A1", #B-Cell Trajectory
#    #"CD14", #Monocytes
#    "PTPRC",
#    "CD3D", 
#    "CD8A", "TBX21", "IL7R" #TCells
#  )
#p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = "PTPRC", embedding = "UMAP",imputeWeights=NULL,quantCut = c(0.01, 0.95),)


#Now that we are confident in the alignment of our scATAC-seq and scRNA-seq, we can label our scATAC-seq clusters with the cell types from our scRNA-seq data.

# First, we will create a confusion matrix between our scATAC-seq clusters and the predictedGroup obtained from our integration analysis.
#cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)


################################################################################################
# identify peaks

# first Add Group Coverages to an ArchRProject object
# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates and then merge these replicates into a single insertion coverage file.

proj <- addGroupCoverages(ArchRProj=proj,groupBy="predictedGroup_Un",force=TRUE)

# generate a reproducible peak set in ArchR
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "predictedGroup_Un", 
    pathToMacs2 = pathToMacs2,
    force=TRUE
)

# To retrieve this peak set as a GRanges object
# This peak set contains an annotation for the group from which each peak originated. 
# However, these annotations do not inherently mean that the given peak was only called in that group, rather that the annotated group had the highest normalized significance for that peak call.
peakGR <- getPeakSet(proj)
peakDF <- cbind(data.frame(seqname=as.vector(seqnames(peakGR)), start=as.vector(start(peakGR)), end=as.vector(end(peakGR))),mcols(peakGR))

write.table(peakDF,"peaks.txt",quote=F,sep="\t",row.names=F)


# Add Peak Matrix

# To prepare for downstream analyses, we can create a new ArchRProject called projHeme5 and add a new matrix to it containing insertion counts within our new merged peak set.
proj <- addPeakMatrix(proj)
# 
getAvailableMatrices(proj)

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)


# Visualizing Genome Browser Tracks 
p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "predictedGroup_Un", 
    geneSymbol = 'FN1', 
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$FN1)

################################################################################################
# Identifying Marker Peaks with ArchR
# It requires PeakMatrix
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup_Un",
    #bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# The object returned by the getMarkerFeatures() function is a SummarizedExperiment that contains a few different assays.
markersPeaks

# We can use the getMarkers() function to retrieve particular slices of this SummarizedExperiment that we are interested in. The default behavior of this function is to return a list of DataFrame objects, one for each cell group.
markerList <- getMarkers(markersPeaks,cutOff = "FDR <= 0.01")
markerDf <- do.call(rbind,
                        lapply(1:length(markerList),
                            function(i){
                                df = markerList[[i]]
                                if(nrow(df)==0){
                                    df=data.frame(seqnames=NULL,idx=NULL,start=NULL,end=NULL,Log2FC=NULL,FDR=NULL,MeanDiff=NULL)
                                }else{
                                    df$group=names(markerList)[i]
                                }
                                df
                            }))

markerDfUP <- as.data.frame(markerDf) %>% dplyr::filter(FDR<= 0.01 & Log2FC >= 0.5)
markerDfDOWN <- as.data.frame(markerDf) %>% dplyr::filter(FDR<= 0.01 & Log2FC <= -0.5)

write.table(markerDf,"marker_FDR0.01.txt",quote=F,sep="\t",row.names=F)
write.table(markerDfUP,"markerUP.txt",quote=F,sep="\t",row.names=F)
write.table(markerDfDOWN,"markerDOWN.txt",quote=F,sep="\t",row.names=F)


################################################################################################
saveRDS(proj,file="MBC6_archr.rds")

################################################################################################
# convert ArchR object to signac object
library(ArchRtoSignac)
source("/home/zhuy1/my_projects/metaplastic_scRNAseq/run_seurat/seurat_workflow/code/ArchRtoSignac.r")

# STEP 0 - Check all required dependencies have been installed and load them automatically.
packages <- c("ArchR","Seurat", "Signac","stringr") # required packages
loadinglibrary(packages)

# STEP 1 - Obtain ArchRProject peak matrix for object conversion.
pkm <- getPeakMatrix(proj)

# STEP 2 - Extract appropriate Ensembl gene annotation and convert to UCSC style.
library(EnsDb.Hsapiens.v86)
annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38") # "UCSC" is the default style to change to but can be changed with argument seqStyle

# STEP 3 - Convert ArchRProject to Signac SeuratObject.
fragments_dir <- "fragments/"

seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = fragments_dir,
  pm = pkm, # peak matrix from getPeakMatrix()
  fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
  fragments_file_extension = ".tsv.gz", # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
  annotation = annotations # annotation from getAnnotation()
)

# STEP 4 - Transfer ArchRProject gene score matrix to Signac SeuratObject.
gsm <- getGeneScoreMatrix(ArchRProject = proj, SeuratObject = seurat_atac)

# STEP 5 - Transfer ArchRProject dimension reduction ("IterativeLSI", "IterativeLSI2" or "Harmony") and UMAP to Signac SeuratObject.
seurat_atac <- addDimRed(
  ArchRProject = proj,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims = "IterativeLSI"
) # default is "IterativeLSI"

# add UMAP embeddings to meta.data
seurat_atac@meta.data <- cbind(seurat_atac@meta.data,seurat_atac@reductions$umap@cell.embeddings)

saveRDS(seurat_atac,file="MBC6_archrData_in_seurat.rds")

################################################################################################


chondro_genes <- c('ACAN',
'BMPR1B',
'CHST3',
'CHST9',
'CHSY1',
'COL27A1',
'CREB3L2',
'CSGALNACT1',
'DDR2',
'EFEMP1',
'FGF2',
'FGFR1',
'FGFR3',
'GLI2',
'GPC6',
'IFT80',
'MMP16',
'NFIB',
'RARB',
'SDC2',
'SNAI2',
'SOX5',
'SOX6',
'TRPS1',
'UST',
'XYLT1')

chondro_genes <- c(chondro_genes,'FN1','SMOC2')

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "predictedGroup_Un", 
    geneSymbol = chondro_genes, 
    features = getPeakSet(proj),
    upstream = 100000,
    downstream = 100000
)
grid::grid.newpage()
grid::grid.draw(p$FN1)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Chondro-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)