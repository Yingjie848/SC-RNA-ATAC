
# secondary analysis for archR results

# Quick tutorial: https://www.archrproject.com/articles/Articles/tutorial.html
# ArchRtoSignac : an Object Conversion Package for ArchR to Signac: https://github.com/swaruplabUCI/ArchRtoSignac

rm(list=ls())

library(ArchR)
set.seed(1)
addArchRThreads(threads = 8)

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)

setwd("../run_ArchR/MBC6_excludeC10/")

outdir="secondary_analysis"; dir.create(outdir)


proj <- readRDS("MBC6_archr.rds")



################################################################################################
# plot browser track 

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

#grid::grid.newpage()
#grid::grid.draw(p$FN1)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Chondro-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 7, height = 5)