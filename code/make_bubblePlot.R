# make gene bubble plot
library(Seurat)
library(ggplot2)
library(ggpubr)

# Args <- commandArgs()
# file_seurat <- Args[6]
# outDir      <- Args[7]

# obj <- readRDS(file_seurat)

make_bubblePlot <- function(obj,outDir){

    dir.create(outDir,recursive = T)

    if(class(obj)!='Seurat' & is.character(obj)){
      if(file.exists(obj))
        obj <- readRDS(file_seurat)
    }

    dir.create(outDir,recursive = T)

    # get normalized counts
    rna_data <- as.matrix(obj$RNA@data)

    panck_genes <- rownames(rna_data)
    panck_genes <- panck_genes[grep('^KRT',panck_genes)]
    panck_genes <- panck_genes[!grepl('CAP',panck_genes)]
    obj@meta.data$PanCK_score <- colMeans(rna_data[panck_genes,])

    obj@meta.data$CD45 <- t(rna_data)[,'PTPRC']
    obj@meta.data$CD68 <- t(rna_data)[,'CD68']
    obj@meta.data$CD163 <- t(rna_data)[,'CD163']
    obj@meta.data$CD3D <- t(rna_data)[,'CD3D']
    obj@meta.data$CD3E <- t(rna_data)[,'CD3E']
    obj@meta.data$CD3G <- t(rna_data)[,'CD3G']
    obj@meta.data$CD3 <- rowMeans(cbind(obj@meta.data$CD3D,obj@meta.data$CD3E,obj@meta.data$CD3G))
    obj@meta.data$GZMA <- t(rna_data)[,'GZMA']
    obj@meta.data$PRF1 <- t(rna_data)[,'PRF1']
    obj@meta.data$CD4 <- t(rna_data)[,'CD4']
    obj@meta.data$CD8A <- t(rna_data)[,'CD8A']
    obj@meta.data$CD8B <- t(rna_data)[,'CD8B']
    obj@meta.data$CD8 <- rowMeans(cbind(obj@meta.data$CD8A,obj@meta.data$CD8B))
    obj@meta.data$CD20 <- t(rna_data)[,'MS4A1']
    
    HLA_genes <- rownames(rna_data)
    HLA_genes <- HLA_genes[grep('^HLA',HLA_genes)]
    obj@meta.data$HLA_score <- colMeans(rna_data[HLA_genes,])
    obj@meta.data$`HLA-A` <- t(rna_data)[,'HLA-A']
    obj@meta.data$`HLA-B` <- t(rna_data)[,'HLA-B']
    obj@meta.data$`HLA-C` <- t(rna_data)[,'HLA-C']

    obj@meta.data$PECAM1 <- t(rna_data)[,'PECAM1']
    obj@meta.data$VWF <- t(rna_data)[,'VWF']

    obj@meta.data$EPCAM <- t(rna_data)[,'EPCAM'] # epithelial
    obj@meta.data$CD49f <- t(rna_data)[,'ITGA6'] # epithelial, CD49f, Pal et al. EMBO 2021, epithelial cells can be CD49f-high(basal) or CD49f+(luminal) or CD49f-(mature luminal)
    obj@meta.data$CD24 <- t(rna_data)[,'CD24'] # epithelial

    obj@meta.data$PDGFRA <- t(rna_data)[,'PDGFRA']
    obj@meta.data$COL5A1 <- t(rna_data)[,'COL5A1']
    obj@meta.data$FN1 <- t(rna_data)[,'FN1']
    obj@meta.data$MKI67 <- t(rna_data)[,'MKI67']

    # myoepithelial markers
    obj@meta.data$ACTA2 <- t(rna_data)[,'ACTA2']
    obj@meta.data$MYH11 <- t(rna_data)[,'MYH11']
    obj@meta.data$CALD1 <- t(rna_data)[,'CALD1']
    obj@meta.data$S100A1 <- t(rna_data)[,'S100A1']
    obj@meta.data$TP63 <- t(rna_data)[,'TP63']
    obj@meta.data$SERPINB5 <- t(rna_data)[,'SERPINB5']
    obj@meta.data$KRT5 <- t(rna_data)[,'KRT5']
    obj@meta.data$KRT7 <- t(rna_data)[,'KRT7']
    obj@meta.data$KRT14 <- t(rna_data)[,'KRT14']
    obj@meta.data$KRT17 <- t(rna_data)[,'KRT17']

    # pericyte
    obj@meta.data$PDGFRB <- t(rna_data)[,'PDGFRB']
    obj@meta.data$NOTCH3 <- t(rna_data)[,'NOTCH3']

    # plasma cells
    obj@meta.data$PIM2 <- t(rna_data)[,'PIM2']

    write.table(obj@meta.data,paste0(outDir,"/metaData.txt"),quote=F,sep="\t",row.names=F)

    #########################################################################################################
    # basic marker genes for immune, macrophage, leukocyte, epithelial, endothelial, fibroblast
    print("Basic markers")
    basic_marker_genes = c(
      'CD45', # CD45
      'CD163', # macrophage antigen
      'CCL8', # macrophage
      'CD14', # leukocyte 
      'FCGR3A', # leukocyte, CD16
      'CD20', # B cells,
      'CD8A', # CD8 T cells
      'CD4', # CD4 T cells 
      'CD2', # T cells
      'CD3D', # T cells
      'CD3E', # T cells
      'CD68', # myeloid cells
      'SPI1', # myeloid
      'FCER1G', # myeloid
      'CSF1R', # myeloid

      'JCHAIN', # plasmablasts

      'PECAM1', # Endothelial cell markers, CD31
      'CLDN5', # Endothelial cell markers, Claudin 5
      'VWF', # Endothelial cell markers

      # epithelia cells
      # EPCAM (CD326) and ITGA6 (CD49f) are used to define subpopulations of epithelial cells. Flow cytometry based on CD49f and EpCAM staining separates lineage-negative breast tissue cells into stromal (CD49f–EpCAM–) and epithelial cells, which includes basal (CD49f+EpCAMlo/–), luminal progenitor (LP) (CD49f+EpCAM+), and mature luminal (ML) (CD49f–EpCAM+) cells.
      'EPCAM', # epithelial
      'CD49f', # epithelial, CD49f, Pal et al. EMBO 2021
      'EGFR', # epithelial

      'MKI67' # mesenchymal cells, cycling cells

      #'KRT5','KRT6A','KRT6B','KRT6C','KRT14', # Pan-CK-highMoleculeWeight
      #'KRT7','KRT8','KRT18', # Pan-CK-lowMoleculeWeight
      #'VIM', # Vimentin

    )
    print(basic_marker_genes[!basic_marker_genes %in% rownames(obj$RNA@data)])

    # // make gene bubble plot
    p1 <- DotPlot(obj, assay='RNA', features = basic_marker_genes,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Marker genes") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.basic_marker_genes.pdf"),8,7); print(p1); dev.off()


    #########################################################################################################
    # core marker genes for immune, macrophage, leukocyte, epithelial, endothelial, fibroblast
    print("Core markers")
    core_marker_genes = c(
      'CD45', # CD45
      'CD68', # myeloid cells
      'CD163', # macrophage antigen
      #'CCL8', # macrophage
      #'CD14', # leukocyte 
      #'FCGR3A', # leukocyte, CD16
      'CD20', # B cells,
      'CD8', # CD8A + CD8B
      'CD8A', # CD8 T cells
      'CD8B', # CD8 T cells
      'CD4', # CD4 T cells 
      'CD3', # T cells
      'CD3D', # T cells
      'CD3E', # T cells
      'CD3G', # T cells
      'GZMA', # cytolytic T cells
      'PRF1', # cytolytic T cells
      'HLA_score',
      
      #'JCHAIN', # plasmablasts

      'PECAM1', # Endothelial cell markers, CD31
      #'CLDN5', # Endothelial cell markers, Claudin 5
      'VWF', # Endothelial cell markers

      'PDGFRA', # fibroblast
      'COL5A1', # fibroblast

      'PDGFRB', # pericyte
      'NOTCH3', # pericyte

      'PIM2', # plasma cells: https://www.sciencedirect.com/science/article/abs/pii/S0165247823000767

      # epithelia cells
      # EPCAM (CD326) and ITGA6 (CD49f) are used to define subpopulations of epithelial cells. Flow cytometry based on CD49f and EpCAM staining separates lineage-negative breast tissue cells into stromal (CD49f–EpCAM–) and epithelial cells, which includes basal (CD49f+EpCAMlo/–), luminal progenitor (LP) (CD49f+EpCAM+), and mature luminal (ML) (CD49f–EpCAM+) cells.
      'EPCAM', # epithelial
      'CD24',
      'CD49f', # epithelial, CD49f, Pal et al. EMBO 2021

      #'EGFR', # epithelial

      'MKI67', # mesenchymal cells, cycling cells

      #'KRT5','KRT6A','KRT6B','KRT6C','KRT14', # Pan-CK-highMoleculeWeight
      #'KRT7','KRT8','KRT18', # Pan-CK-lowMoleculeWeight
      #'VIM', # Vimentin
      'PanCK_score'
      

    )
    print(core_marker_genes[!core_marker_genes %in% rownames(obj$RNA@data)])

    # // make gene bubble plot
    p1 <- DotPlot(obj, assay='RNA', features = core_marker_genes,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Marker genes") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.core_marker_genes.pdf"),10,7); print(p1); dev.off()

    # make violin plot
    plots <- lapply(core_marker_genes,function(feature){
       VlnPlot(obj,features=feature)
    })
    p <- ggarrange(plotlist=plots,nrow=2,ncol=1)
    pdf(paste0(outDir,"/gene_VlnPlot.core_marker_genes.pdf"),7,7); print(p); dev.off()

    

    #########################################################################################################
    # fibroblast markers
    print("Fibroblast markers")
    genes_for_fibroblast <- c(
      # The short-list of common fibroblast markers, Muhl et al. Nature Communications, 2020, Results https://www.nature.com/articles/s41467-020-17740-1
      'PDGFRA','FBLN1','FBLN2','LOXL1','LUM',
      #'PDGFRB',

      'FBN1', # HPA, Fibroblasts - ECM organization (mainly), https://www.proteinatlas.org/ENSG00000166147-FBN1

      # markers from Pilling et al. Plos One 2019 Introduction
      'FAP',
      #'HABP4', # HA-BP, hyaluronan
      #'ITGA3', # CD49C
      'FN1', # cellular fibronectin
      #'S100A4' # FSP1
      'THY1',  #CD90

      # fibroblast markers, Pal et al. EMBO, 2021, Figure 2F
      'MMP2',
      'CLMP',
      'DPT',
      'DCN',
      'TWIST2',

      # 'COL5A1', # it expressed in fibroblast cells in breast according to HPA: https://www.proteinatlas.org/ENSG00000130635-COL5A1/single+cell+type
      'COL5A1'
      
    )
    genes_for_fibroblast[!genes_for_fibroblast %in% rownames(obj$RNA@data)]
    table(genes_for_fibroblast)

    # // make gene bubble plot
    p2 <- DotPlot(obj, assay='RNA', features = genes_for_fibroblast,group.by='seurat_clusters',) + RotatedAxis() + ggtitle("Fibroblast") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_fibroblasts.pdf"),6,7); print(p2); dev.off()

    #########################################################################################################
    # Properties used to identify fibrocytes, macrophages and fibroblasts
    # Reilkoff et al. Nature Reviews, immunology https://www.nature.com/articles/nri2990
    # Pilling et al Plos One: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2759556/#:~:text=CD45%20and%20LSP%2D1%20are,27%5D%E2%80%93%5B30%5D.
    print("Properties used to identify fibrocytes, macrophages and fibroblasts")
    genes_for_fibrocytes_macrophages_fibroblasts <- c(
      
      # Adhesion and motility markers
      'LSP1',
      'CD9',
      'ITGAL', # CD11a
      'ITGAM', # CD11b
      'ITGAX', # CD11c
      'SPN', # CD43
      'CD164',
      'LGALS3', # galectin 3
      'CD34',
      # Cell surface enzymes
      'MME', # CD10
      'SIRPA', # CD172a
      'CD45', # CD45
      'FAP',
      # Scavenging receptors and molecules involved in host defence
      'CD14', 'CD68', 'CD163', 
      'MRC1', # CD206
      'CD209', 
      'CR1', # CR1
      'CD36',
      # Fc γ receptors
      'FCGR3A', # CD16
      'FCGR2A', # CD32a
      'FCGR2B', # CD32b
      'FCGR2C', # CD32c
      # Chemokine receptors
      'CCR4', 'CCR5', 'CCR7', 'CXCR1', 'CXCR4', 'CX3CR1',
      # Cell-surface molecules involved in antigen presentation
      'ICAM1', # CD54
      'CD80', 'CD86','HLA-A','HLA-DQA1',
      # Extracellular matrix proteins
      #'COL1A1','COL1A2',# Collagen I
      #'COL3A1', # Collagen III
      #'COL4A1', # Collagen IV
      'VIM', # vimentin
      'TNC', # tenascin

      # Endothelial cell markers
      'PECAM1', # Endothelial cell markers, CD31
      'CLDN5', # Endothelial cell markers, Claudin 5
      'VWF', # Endothelial cell markers

      # markers from Pilling et al.
      'HABP4', # HA-BP, hyaluronan
      'ITGA3', # CD49C
      'FN1', # cellular fibronectin,
      
      # 
      'PDGFRA','FBLN1','FBLN2'
    )
    genes_for_fibrocytes_macrophages_fibroblasts[!genes_for_fibrocytes_macrophages_fibroblasts %in% rownames(obj$RNA@data)]
    table(genes_for_fibrocytes_macrophages_fibroblasts)

    # // make gene bubble plot
    p3 <- DotPlot(obj, assay='RNA', features = genes_for_fibrocytes_macrophages_fibroblasts,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Fibrocytes, Macrophages, Fibroblasts") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_fibrocytes_macrophages_fibroblasts.pdf"),15,7); print(p3); dev.off()

    #########################################################################################################
    # CK
    print("CKs")
    CK_genes = c(
      
      'KRT5','KRT6A','KRT6B','KRT6C','KRT14', # Pan-CK-highMoleculeWeight
      'KRT7','KRT8','KRT18' # Pan-CK-lowMoleculeWeight

    )
    CK_genes[!CK_genes %in% rownames(obj$RNA@data)]

    # // make gene bubble plot
    p4 <- DotPlot(obj, assay='RNA', features = CK_genes,group.by='seurat_clusters',scale=TRUE,scale.by='radius') + RotatedAxis() + ggtitle("CKs") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.CKs.pdf"),6,7); print(p4); dev.off()

    #########################################################################################################
    # cytotoxic T cells gene list
    print("Cytotoxic T cells")
    cytotoxic_genes <- c(
      'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM', 'GNLY', 'PRF1','FASLG',
      # cytotoxic T cell activation markers
      'IFNG', 'TNF','IL2R','IL2'
    )
    cytotoxic_genes[!cytotoxic_genes %in% rownames(obj$RNA@data)]

    # // make gene bubble plot
    p5 <- DotPlot(obj, assay='RNA', features = cytotoxic_genes,group.by='seurat_clusters',scale=TRUE,scale.by='radius') + RotatedAxis() + ggtitle("Cytotoxic T cells") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.cytotoxic_T_cells.pdf"),6,7); print(p5); dev.off()


    #########################################################################################################
    # Tumor associated macrophages gene list
    # Cassetta, Cancer Cell, 2019
    print("TAM")
    TAM_genes <- c(
      'IRF8','CCL2','C1QC','GBP5','HCST','LILRB4','AIF1','PSMB9','GBP4','GBP1','HLA-DOA','C1QA','CCL4','NCF1C','LAP3','TNFAIP3','ITGB2','LAIR1','FOLR2','CD83','SIGLEC1','TCN2','PLTP','C1QB','DOK2',
      'GIMAP6','CD40','CCL3','CCL8','FCN1','CD4','VAV1','TLR7','FGD2','LST1','VSIG4','CLEC7A'
    )
    TAM_genes[!TAM_genes %in% rownames(obj$RNA@data)]

    # // make gene bubble plot
    p6 <- DotPlot(obj, assay='RNA', features = TAM_genes,group.by='seurat_clusters',scale=TRUE,scale.by='radius') + RotatedAxis() + ggtitle("TAM") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.TAM.pdf"),12,7); print(p6); dev.off()


    #########################################################################################################
    # fibroblast markers, Pal et al. EMBO, 2021, Figure 2F
    print("Fibroblast markers, Pal et al. EMBO, 2021, Figure 2F")
    genes_for_fibroblast_Pal <- c(
      'MMP2', 'APOD', 'CLMP','PDGFRA', 'MMP3', 'TWIST2','DPT', 'IGFBP6', 'DCN', 'LUM'
    )
    genes_for_fibroblast_Pal[!genes_for_fibroblast_Pal %in% rownames(obj$RNA@data)]
    table(genes_for_fibroblast_Pal)

    # // make gene bubble plot
    p7 <- DotPlot(obj, assay='RNA', features = genes_for_fibroblast_Pal,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Fibroblast") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_fibroblast_Pal.pdf"),6,7); print(p7); dev.off()


    #########################################################################################################
    # endothelial markers, Pal et al. EMBO, 2021, Figure 2F
    print("Endothelial markers, Pal et al. EMBO, 2021, Figure 2F")
    genes_for_endothelial_Pal <- c(
      'TFF3',
      'CYLD',
      'MMRN1',
      'CLIC2',
      'TFP',
      'SOX4',
      'GJA1',
      'PTPRB',
      'CDH5',
      'ERG',
      'TSPAN12',
      'HEY1',
      'PRRG4',
      # 'COL5A1', # it expressed in fibroblast cells in breast according to HPA: https://www.proteinatlas.org/ENSG00000130635-COL5A1/single+cell+type
      'SNAI1',
      'SPARC',
      'ENG',
      'FLT1',
      'CD93',
      'PLVAP',
      'SOX17'
    )
    genes_for_endothelial_Pal[!genes_for_endothelial_Pal %in% rownames(obj$RNA@data)]
    table(genes_for_endothelial_Pal)

    # // make gene bubble plot
    p8 <- DotPlot(obj, assay='RNA', features = genes_for_endothelial_Pal,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Endothelial") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_Endothelial_Pal.pdf"),7,7); print(p8); dev.off()

    #########################################################################################################
    # myoepithelial markers, Batistatou et al. In vivo, 2003
    print("Myoepithelial markers, Batistatou et al. In vivo, 2003")
    genes_for_myoepithelial <- c(
      'TP63','ACTA2','MYH11','CALD1','S100A1','SERPINB5','KRT5','KRT7','KRT14','KRT17'
    )
    genes_for_myoepithelial[!genes_for_myoepithelial %in% rownames(obj$RNA@data)]
    table(genes_for_myoepithelial)

    # // make gene bubble plot
    p <- DotPlot(obj, assay='RNA', features = genes_for_myoepithelial,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Myoepithelial") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_myoepithelial.pdf"),7,7); print(p); dev.off()


    #########################################################################################################
    # Natural Killer cells
    print("Natural Killer")
    genes_for_NK <- c(
      'NCAM1', # CD56'
      'KLRB1', # CD161
      'FCGR3A', # CD16
      'KLRD1', # CD94
      'B3GAT1' # CD57
    )
    genes_for_NK[!genes_for_NK %in% rownames(obj$RNA@data)]
    table(genes_for_NK)

    # // make gene bubble plot
    p9 <- DotPlot(obj, assay='RNA', features = genes_for_NK,group.by='seurat_clusters') + RotatedAxis() + ggtitle("NK") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_NK.pdf"),4,7); print(p9); dev.off()


    #########################################################################################################
    # B cells
    # https://www.abcam.com/primary-antibodies/b-cells-basic-immunophenotyping#:~:text=For%20most%20mature%20B%20cells,CD30%2C%20a%20regulator%20of%20apoptosis.
    print("B cells")
    genes_for_Bcell <- c(
      # all
      'IGHM', # IgM. Immature B cells express CD19, CD 20, CD34, CD38, and CD45R, but not IgM
      # memory cell
      'CD20', # CD20
      # activated B cell
      'TNFRSF8', # CD30
      'CD19', 
      'IL2RA', # CD25
      'CD79A',
      'CD79B'
    )
    genes_for_Bcell[!genes_for_Bcell %in% rownames(obj$RNA@data)]
    table(genes_for_Bcell)

    # // make gene bubble plot
    p10 <- DotPlot(obj, assay='RNA', features = genes_for_Bcell,group.by='seurat_clusters') + RotatedAxis() + ggtitle("B cells") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_B_cells.pdf"),5,7); print(p10); dev.off()

    #########################################################################################################
    # breast cell origin markers
    # Pal et al TMBO 2021
    print("Breast cell origin markers")
    genes_for_cell_origin <- c(
      # Pal et al. EMBO 2021 figure 1
      # EPCAM (CD326) and ITGA6 (CD49f) are used to define subpopulations of epithelial cells. 
      # Flow cytometry based on CD49f and EpCAM staining separates lineage-negative breast tissue cells into stromal (CD49f–EpCAM–) and epithelial cells, 
      # which includes 
      # basal (CD49f+EpCAMlo/–), luminal progenitor (LP) (CD49f+EpCAM+), and mature luminal (ML) (CD49f–EpCAM+) cells.
      'EPCAM', # epithelial, high in luminal, low in basal
      'CD49f', # CD49f, epithelial, expressed in basal and LP, not expressed in ML 
      'KRT5', # high in basal, low in luminal
      'ACTA2', # high in basal, low in luminal
      'MYLK', # high in basal, low in luminal
      'SNAI2', # high in basal, low in luminal
      'NOTCH4', # high in basal, low in luminal 
      'DKK3', # high in basal, low in luminal 
      'ESR1', # high in ML
      'PGR', # high in ML
      'FOXA1', # high in ML
      'TNFRSF11A', # high in LP
      'KIT', # high in LP
      'SOX10', # high in LP
      # Pal et al. Nat Commun, 2017. Basal cells (Lin–CD29hiCD24+) and luminal cells (Lin–CD29loCD24+)
      'CD24', 
      'ITGB1', # CD29, 
      # Pal et al. Nat Commun, 2017. CD55 marks distinct populations through development
      'CD55', 
      # canonical basal marker
      'KRT14', 
      'VIM', 
      # canonical luminal marker
      'KRT8',  
      'KRT18' 

    )
    print(genes_for_cell_origin[!genes_for_cell_origin %in% rownames(obj$RNA@data)])
    table(genes_for_cell_origin)

    # // make gene bubble plot
    p11 <- DotPlot(obj, assay='RNA', features = genes_for_cell_origin,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Cell origin markers") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_cell_origin.pdf"),8,7); print(p11); dev.off()


    #########################################################################################################
    # chondroid markers
    # Daugaard APMIS 2009. PMID: 19594492
    # osteonectin, bcl-2, cox-2, actin, calponin, D2-40 (podoplanin), mdm-2, CD117 (c-kit) and YKL-40
    print("Chondroid markers")
    genes_for_chondroid <- c(
      'SPARC', # osteonectin
      'BCL2', # Bcl2
      'PTGS2', # cox-2
      'ACTA1', # Actin Alpha1
      'ACTA2', # Actin Alpha2
      'ACTB', # Actin beta
      'ACTC1',
      'ACTG1',
      'ACTG2',
      'ACTR3',
      'CNN1', # Calponin-1
      'CNN2', # Calponin-2
      'CNN3', # Calponin-3
      'PDPN', # D2-40 (podolanin)
      'MDM2', # mdm-2
      'KIT', # CD117
      'CHI3L1' # Chitinase 3 Like 1, YKL-40
    )
    print(genes_for_chondroid[!genes_for_chondroid %in% rownames(obj$RNA@data)])
    table(genes_for_chondroid)

    # // make gene bubble plot
    p12 <- DotPlot(obj, assay='RNA', features = genes_for_chondroid,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Chondroid markers") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_chondroid.pdf"),8,7); print(p12); dev.off()


    #########################################################################################################
    # tumor initiating and epithelial to mesenchymal transition
    # Zhang, Modern Pathology, 2011
    # Gibbons, Dev Dyn, 2019. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5503821/#SD1
    # Creighton et al., 2013. EMT-ness: VIM + CDH2 + FOXC2 + SNAI1 + SNAI2 + TWIST1 + FN1 + ITGB6 + MMP2 + MMP3 + MMP9 + SOX10 + GCS − CDH1 − DSP − OCLN 
    print("Tumor initiating and epithelial to mesenchymal transition")
    genes_for_TI_EMT <- c(
      'ALDH1A1', # ALDH1
      'CD44',
      'CD24',
      # mesenchymal markers
      'VIM', # vimentin
      'CDH2', # N-cadherin
      'FOXC2',
      'SNAI1',
      'SNAI2',
      'TWIST1',
      'FN1', # fibronectin
      'ITGB6', #integrin αvβ6
      'MMP2',
      'MMP3',
      'MMP9',
      'SOX10',
      'GCLC', # GCS
      # epithelial markers
      'CDH1', # E-cadherin
      'DSP', # desmoplakin
      'OCLN', # occludin
      # EMT
      'ZEB1',
      'ZEB2'
    )
    print(genes_for_TI_EMT[!genes_for_TI_EMT %in% rownames(obj$RNA@data)])
    table(genes_for_TI_EMT)

    # // make gene bubble plot
    p13 <- DotPlot(obj, assay='RNA', features = genes_for_TI_EMT,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Tumor initiating and EMT") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_TI_EMT.pdf"),8,7); print(p13); dev.off()


    #########################################################################################################
    # spindle cells
    # Spindle cell lesions: A review on immunohistochemical markers
    # Surbhi, Rashmi Metgud, Smitha Naik, Shrikant Patel, 2017
    # https://www.cancerjournal.net/article.asp?issn=0973-1482;year=2017;volume=13;issue=3;spage=412;epage=418;aulast=Surbhi,#:~:text=Spindle%20cells%20may%20be%20focally,those%20at%20high%20malignant%20risk. 

    print("Spindle cells")
    genes_for_spindle <- c(
      'VIM', # vimin
      'SLC4A1', # AE1
      'SLC4A3', # AE3
      'MUC1', # EMA
      'TP63', # p63
      'PLAG1',
      'KRT7',
      'KRT14',
      'CEACAM1', # CEA
      'CEACAM2', # CEA
      'CEACAM3', # CEA
      'CEACAM4', # CEA
      'CEACAM5', # CEA
      'MYH11', # SMMHC
      'CNN1', # Calponin-1
      'CNN2', # Calponin-2
      'CNN3', # Calponin-3
      'S100A1', # S-100
      'S100B', # S-100
      'WT1', # Wilms tumor 1
      'GFAP'
    )
    print(genes_for_spindle[!genes_for_spindle %in% rownames(obj$RNA@data)])
    table(genes_for_spindle)

    # // make gene bubble plot
    p14 <- DotPlot(obj, assay='RNA', features = genes_for_spindle,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Spindle markers") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_spindle_markers.pdf"),8,7); print(p14); dev.off()
    
    #########################################################################################################
    # breast epithelial subtype markers, Pal et al. marker genes for each epithelial cluster
    # Pal et al TMBO 2021
    print("Breast epithelial subtype markers")
    genes_for_epi_subtype <- c(
      'LAMA3','EPAS1','ID3','ERG', 'TPM2', 'COL17A1','OXTR', 'POSTN', 'WIF1','SNAI2', 'DLL1', 'COL4A2','CRYAB', 'ID4', 'SPHK1','MMP3', 'TAGLN', 'KRT14', # basal
      'UGDH', 'ALCAM', 'IRS1','ENPP1', 'KLHL5', 'HERC4','TOX3', 'UCP2', 'PREX1','TRAFD1', 'PKIB', # ML
      'LTF', 'SERPINB4', 'CMPK1','IL27RA', 'PML', 'SESTD1','CXCL16', 'MFGE8', 'CAVIN3','SFRP1', 'PROM1','SOX10', 'FGFBP1' # LP
    )
    print(genes_for_epi_subtype[!genes_for_epi_subtype %in% rownames(obj$RNA@data)])
    table(genes_for_epi_subtype)

    # // make gene bubble plot
    p15 <- DotPlot(obj, assay='RNA', features = genes_for_epi_subtype,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Epithelial subtype markers (Pal et al.)") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_cell_origin_Pal.pdf"),15,7); print(p15); dev.off()


    #########################################################################################################
    # monocytes
    # https://www.bio-rad-antibodies.com/monocyte-cd-markers-antibodies.html?JSESSIONID_STERLING=B8E76C09819A31702331CECE39FBD52D.ecommerce2&evCntryLang=US-en&cntry=US&thirdPartyCookieEnabled=true
    print("Monocytes")
    genes_for_monocytes <- c(
      'CD2',
      'ITGAM', # CD11b
      'CD14',
      'FCGR3A', # CD16
      'PECAM1', # CD31
      'NCAM1', # CD56
      'SELL', # CD62L
      'CSF1R', # CD115
      'CCR2', # CD192
      'CX3CR1',
      'CXCR3',
      'CXCR4'
    )
    print(genes_for_monocytes[!genes_for_monocytes %in% rownames(obj$RNA@data)])
    table(genes_for_monocytes)

    # // make gene bubble plot
    p16 <- DotPlot(obj, assay='RNA', features = genes_for_monocytes,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Epithelial subtype markers (Pal et al.)") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_monocytes.pdf"),8,7); print(p16); dev.off()

    #########################################################################################################
    # stem cell markers
    # Wu ONCOLOGY REPORTS 2019 Table SI
    print("stem cell markers")
    genes_for_stemness <- c(
        'EPHB2',
        'LATS2',
        'CXCL2',
        'MLLT10',
        'METTL3',
        'EZH2',
        'VEGFA',
        'EGFR',
        'GPSM2',
        'DNMT3A',
        'HES1',
        'MYC',
        'RNF43',
        'AXIN2',
        'CDK6',
        'RGMB',
        'CDCA7',
        'AQP1',
        'DPP4',
        'CFTR',
        'TSPAN6',
        'FERMT1',
        'LRIG1',
        'CA2',
        'STMN1',
        'ASCL2',
        'OLFM4',
        'LGR5',
        'BMI1',
        'SMOC2',
        'ETS2',
        'GPX2',
        'FAM84A',
        'PTPRG',
        'SLC12A2'
    )

    print(genes_for_stemness[!genes_for_stemness %in% rownames(obj$RNA@data)])
    table(genes_for_stemness)

    # // make gene bubble plot
    p16 <- DotPlot(obj, assay='RNA', features = genes_for_stemness,group.by='seurat_clusters') + RotatedAxis() + ggtitle("Stem cell markers") + xlab("") + theme(legend.position = "right")
    pdf(paste0(outDir,"/gene_bubblePlot_rna.genes_for_stemness.pdf"),10,7); print(p16); dev.off()

}

