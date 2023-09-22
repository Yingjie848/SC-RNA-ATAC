# reduction_methods: umap, wnn.umap
library(Seurat)
library(dplyr)
make_featurePlots <- function(obj,outDir,reduction_method='umap',order=FALSE){

    # normalize PanCK and CD45
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
        obj@meta.data$CD3 <- base::rowMeans(cbind(obj@meta.data$CD3D,obj@meta.data$CD3E,obj@meta.data$CD3G))
        obj@meta.data$GZMA <- t(rna_data)[,'GZMA']
        obj@meta.data$PRF1 <- t(rna_data)[,'PRF1']
        obj@meta.data$CD4 <- t(rna_data)[,'CD4']
        obj@meta.data$CD8A <- t(rna_data)[,'CD8A']
        obj@meta.data$CD8B <- t(rna_data)[,'CD8B']
        obj@meta.data$CD8 <- base::rowMeans(cbind(obj@meta.data$CD8A,obj@meta.data$CD8B))
        obj@meta.data$CD20 <- t(rna_data)[,'MS4A1']
        obj@meta.data$PIM2 <- t(rna_data)[,'PIM2']
        
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

        # macrophage
        obj@meta.data$AIM2 <- t(rna_data)[,'AIM2'] # AIM2  The overexpression of AIM2 in macrophages resulted in inflammasome signaling and promoted caspase-1 activation to produce the pro-inflammatory cytokine IL-1β. These cytokines then exert promoting effects on the switch from macrophage M2 polarization phenotype to M1.
        obj@meta.data$IL4R <- t(rna_data)[,'IL4R'] # https://www.cusabio.com/c-20938.html#:~:text=In%20the%20most%20distinct%20genes,1%20and%20CD206%20phenotype%20markers.

        # pericytes markers: https://www.ahajournals.org/doi/full/10.1161/01.res.0000182903.16652.d7#:~:text=Pericyte%20Identification,-The%20morphological%20diversity&text=Several%20markers%20have%20been%20used,the%20promoter%20trap%20transgene%20XlacZ4.
        # However, none of these markers is absolutely specific for pericytes, and none of the markers recognizes all pericytes; their expression is dynamic and varies between organs and developmental stages.
        obj@meta.data$ACTA2 <- t(rna_data)[,'ACTA2'] #  smooth muscle α-actin (SMA)
        obj@meta.data$PDGFRB <- t(rna_data)[,'PDGFRB'] # platelet-derived growth factor receptor (PDGFR)-β
        obj@meta.data$DES <- t(rna_data)[,'DES'] # desmin, DES
        obj@meta.data$CSPG4 <- t(rna_data)[,'CSPG4'] # NG-2, CSPG4
        obj@meta.data$ENPEP <- t(rna_data)[,'ENPEP'] # aminopeptidase A, APA
        obj@meta.data$ANPEP <- t(rna_data)[,'ANPEP'] # aminopeptidase N, ANPEP
        obj@meta.data$RGS5 <- t(rna_data)[,'RGS5'] # RGS5
        # PDGFRβ co-localized only with well-established pericyte markers such as Chondroitin Sulfate Proteoglycan NG2 and the xLacZ4 transgenic reporter. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2936891/
        obj@meta.data$RGS5 <- t(rna_data)[,'CSPG4'] # Chondroitin Sulfate Proteoglycan NG2



        pdf(paste0(outDir,'/FeaturePlot.pdf'),5,5)
        print(FeaturePlot(obj,features='nCount_RNA',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='nFeature_RNA',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='percent.mt',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='MALAT1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='percent.ribo',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='percent.hb',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD45',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PanCK_score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='hybrid_score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='cxds_score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='bcds_score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='scDblFinder.score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='S.Score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='G2M.Score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD68',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD163',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD3',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='GZMA',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PRF1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD4',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD8',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='HLA_score',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD20',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD24',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PECAM1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='VWF',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='EPCAM',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CD49f',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PDGFRA',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='COL5A1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PDGFRB',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='NOTCH3',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PIM2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='FN1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='MKI67',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.ordered.pdf'),5,5)
        print(FeaturePlot(obj,features='CD68',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD163',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD3',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='GZMA',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='PRF1',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD4',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD8',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='HLA_score',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD20',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD24',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='PECAM1',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='VWF',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='EPCAM',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='CD49f',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='PDGFRA',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='COL5A1',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='PDGFRB',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='NOTCH3',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='PIM2',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='FN1',reduction=reduction_method,label=T,order=TRUE))
        print(FeaturePlot(obj,features='MKI67',reduction=reduction_method,label=T,order=TRUE))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Mesenchymal.pdf'),5,5)
        print(FeaturePlot(obj,features='CDH2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='VIM',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='FN1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='ITGB6',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Epithelial.pdf'),5,5)
        print(FeaturePlot(obj,features='CDH1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='DSP',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='OCLN',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Myoepithelial.pdf'),5,5)
        print(FeaturePlot(obj,features='ACTA2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='MYH11',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CALD1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='S100A1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='TP63',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='SERPINB5',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='KRT5',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='KRT7',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='KRT14',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='KRT17',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.EMT.pdf'),5,5)
        print(FeaturePlot(obj,features='SNAI1',reduction=reduction_method,label=T,order=order)) # don't have this gene
        print(FeaturePlot(obj,features='SNAI2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='ZEB1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='ZEB2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='TWIST1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='TWIST2',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.MaSC.pdf'),5,5)
        print(FeaturePlot(obj,features='FN1',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CALD1',reduction=reduction_method,label=T,order=order)) # don't have this gene
        print(FeaturePlot(obj,features='SPTBN1',reduction=reduction_method,label=T,order=order))
        #print(FeaturePlot(obj,features='RICS',reduction=reduction_method,label=T,order=order))
        #print(FeaturePlot(obj,features='MSCP',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='MAML2',reduction=reduction_method,label=T,order=order))
        #print(FeaturePlot(obj,features='KIAA1718',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Macrophages.pdf'),5,5)
        print(FeaturePlot(obj,features='AIM2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='IL4R',reduction=reduction_method,label=T,order=order))
        dev.off()

        pdf(paste0(outDir,'/FeaturePlot.Pericytes.pdf'),5,5)
        print(FeaturePlot(obj,features='ACTA2',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='PDGFRB',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='DES',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CSPG4',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='ENPEP',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='ANPEP',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='RGS5',reduction=reduction_method,label=T,order=order))
        print(FeaturePlot(obj,features='CSPG4',reduction=reduction_method,label=T,order=order))

        pericyte_signature <- 'PECAM1,PDGFRB,CSPG4,ANPEP,ACTA2,DES,RGS5,ABCC9,KCNJ8,CD248,DLK1,TEK,NOTCH3,GLI1,ICAM1,ADM,ANG2,ANGPT1,VEGFA,ZIC1,FOXC1,POSTN,COX4I2,LHFP,HIGD1B,PDZD2,HSD11B1,MCAM,MXRA8,PDE5A,NR1H3,SERPING1,EMID1,ECM1,COLEC11,RARRES2,REM1,ASPN,CYGB,FABP4,VTN,STEAP4,NDUFA4L2,SLC38A11,ATP13A5,AOC3,ANGPT2,INPP4B,GPIHBP1,VIM,PTH1R,IFITM1,TBX18,NT5E,MFGE8,ALPL,COL1A1,MYO1B,COG7,P2RY14,HEYL,GNB4,MSX1,CTGF'
        pericyte_signature <- strsplit(pericyte_signature,',')[[1]]
        pericyte_signature <- pericyte_signature[pericyte_signature %in% rownames(obj$RNA@data)]

        for(perisig in pericyte_signature){
            print(perisig)
            print(FeaturePlot(obj,features=perisig,reduction=reduction_method,label=T,order=order))
        }
        dev.off()

        #make_featurePlots_for_chondroid_markers(obj,outDir,reduction_method='umap')



}

make_featurePlots_for_chondroid_markers <- function(obj,outDir,reduction_method='umap',order=FALSE){

    dir.create(outDir)

    DefaultAssay(obj) <- 'RNA'

    # Daugaard APMIS 2009. PMID: 19594492
    # osteonectin, bcl-2, cox-2, actin, calponin, D2-40 (podoplanin), mdm-2, CD117 (c-kit) and YKL-40
    genesets_Daugaard <- c(
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

    genesets_Daugaard[!genesets_Daugaard %in% rownames(obj$RNA@data)]

    pdf(paste0(outDir,'/FeaturePlot.chondroid_Daugaard.pdf'),5,5)
    for(gene in genesets_Daugaard){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()

    # GO_CHONDROCYTE_DIFFERENTIATION

    genesets_GO_chondrocyte_differentiation <- c(
        'SIX2','ADAMTS7','OSR2','COL2A1','COL6A1','COL6A2','COL6A3','COL7A1','COL11A1', 
        'COL12A1','OSR1','COMP','MAPK14','RFLNA','CCN2','CTNNB1','CHADL','COCH','ACAN','ECM1',
        'EFEMP1','FGF9','FGFR1','FGFR3','CHSY1','SULF1','POC1A','AMELX','GREM1','GLG1','GLI2',
        'GLI3','GPLD1','ANXA2','ANXA2P2','ANXA6','HOXA11','VWA2','IHH','RFLNB','MUSTN1','GDF6',
        'SNX19','LOXL2','LTBP3','SMAD3','MAF','MATN1','MATN2','MATN3','MBL2','MDK','MEF2C','MEF2','DMSX2',
        'NFIB','CCN3','NPPC','CHST11','ZNF219','MEX3C','SCARA3','VIT','CYTL1','POR','IMPAD1','SMPD3','SOX6',
        'SULF2','PTH','PTHLH','PTH1R','IFT80','COL20A1','NKX3-2','RARB','RARG','RELA','SFRP2','SCX','SHOX2',
        'CREB3L2','VWA1','BMP2','BMP4','BMP6','BMPR1A','BMPR1','BBMPR2','SNAI2','SOX5','SOX9','TGFB1','TGFBI',
        'TGFBR1','TGFBR2','TRPS1','COL14A1','WNT7A','WNT10B','WNT2B','WNT9A','ZBTB16','LNPK','HMGA2','WNT5B',
        'ADAMTS12','GDF5','AXIN2','ANXA9','COL27A1','SCIN','RUNX2','RUNX1','RUNX3','SERPINH1','MATN4','FGF18','CCN4','PKDCC','TRIP11','ACVRL1','EIF2AK3'
    )

    genesets_GO_chondrocyte_differentiation <- genesets_GO_chondrocyte_differentiation[genesets_GO_chondrocyte_differentiation %in% rownames(obj$RNA@data)]
    
    pdf(paste0(outDir,'/FeaturePlot.GO_CHONDROCYTE_DIFFERENTIATION.pdf'),5,5)
    for(gene in genesets_GO_chondrocyte_differentiation){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()


    # chondroid-related genes in Daugaard and GO, expressed in MBC5 tumors
    genesets_expressed_in_MBC5_tumors <- c(
        'CNN2','CNN3','CHI3L1', # Daugaard 
        'COL2A1','COL6A1','COL6A2','COL6A3','COL7A1','RFLNA','CCN2','ACAN','FGFR1','FGFR3','GLI2',
        'GLI3','ANXA2','MATN2','NFIB','SOX6','PTH1R','RARB','RARG','CREB3L2','SNAI2','SOX5','SOX9','TRPS1','COL27A1','SCIN'
    )

    genesets_expressed_in_MBC5_tumors <- sort(genesets_expressed_in_MBC5_tumors)

    pdf(paste0(outDir,'/FeaturePlot.chondroid_expressed_in_MBC5.pdf'),5,5)
    for(gene in genesets_expressed_in_MBC5_tumors){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()

    # chondroid-related genes in Daugaard and GO, expressed in MBC6 tumors
    genesets_expressed_in_MBC6_tumors <- c(
        'ACTG1','CNN2','CNN3','KIT','CHI3L1', # Daugaard 
        'COL6A1','ACAN','FGFR3','NFIB','SOX6','RARB','CREB3L2','SNAI2','TRPS1','SCIN',
        'CSGALNACT1','XYLT1','CHSY1','SDC2' # REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM
    )

    genesets_expressed_in_MBC6_tumors <- sort(genesets_expressed_in_MBC6_tumors)

    pdf(paste0(outDir,'/FeaturePlot.chondroid_expressed_in_MBC6.pdf'),5,5)
    for(gene in genesets_expressed_in_MBC6_tumors){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()


    # chondrocytes related genes in PanglaoDB
    chondrocyte_genes <- read.csv("seurat_workflow/data/PanglaoDB_markers_27_Mar_2020.tsv",header=T,sep="\t") %>% dplyr::filter(cell.type=="Chondrocytes") %>% dplyr::arrange(official.gene.symbol) %>% .[['official.gene.symbol']]

    chondrocyte_genes <- chondrocyte_genes[chondrocyte_genes %in% rownames(obj$RNA@data)]

    pdf(paste0(outDir,'/FeaturePlot.PanglaoDB_chondrocytes.pdf'),5,5)
    for(gene in chondrocyte_genes){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()


}



make_featurePlots_for_myoepithelial_markers <- function(obj,outDir,reduction_method='umap',order=FALSE){

    dir.create(outDir)

    DefaultAssay(obj) <- 'RNA'

    # myoepichelial markers from cellMarkers database
    genesets_cellMarkers <- c(
      'ITGB1','CD44','TP63','CNN1','CNN2','CNN3','MME','KRT14','SMN1','LGALS7','CDH3','ALPG'
    )

    genesets_cellMarkers[!genesets_cellMarkers %in% rownames(obj$RNA@data)]

    pdf(paste0(outDir,'/FeaturePlot.myoepithelial_cellMarkers.pdf'),5,5)
    for(gene in genesets_cellMarkers){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()

    # myoepithelial markers, Batistatou et al. In vivo, 2003
    genes_for_myoepithelial <- c(
      'TP63','ACTA2','MYH11','CALD1','S100A1','SERPINB5','KRT5','KRT7','KRT14','KRT17'
    )

    genes_for_myoepithelial[!genes_for_myoepithelial %in% rownames(obj$RNA@data)]

    pdf(paste0(outDir,'/FeaturePlot.myoepithelial_Batistatou_2003.pdf'),5,5)
    for(gene in genes_for_myoepithelial){
        print(gene)
        print(FeaturePlot(obj,features=gene,reduction=reduction_method,label=T,order=order))
    }
    dev.off()

}