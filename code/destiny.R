# To estimate pseudotime, we use destiny package to build diffusion map and use it to estimate pseudotime
# useful webpage: https://danhdtruong.com/Diffusion-Pseudotime/
#                 https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html#comparison-of-the-different-trajectory-inference-methods
library(destiny)
library(Seurat) 
library(SingleCellExperiment)
library(ggplot2)
library(slingshot)

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
options(ggplot2.discrete.colour = cluster_cols)


##############################################################################################################
# here according to http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/diffusionMaps

destiny_mit <- function(sce,outDir){

    dir.create(outDir)

    #this has the cell classification
    table(sce$ident)

    # Running diffusion map
    #this step may take a long time (days) or not finish. It is recommend to send it to the cluster as a script that reads the Seurat or the single cell object, runs DiffusionMap, and saves the object. 
    dm <- DiffusionMap(sce, verbose = TRUE, n_pcs = 50)

    # Plotting the diffusion map
    cellLabels <- sce$ident
    tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                      DC2 = eigenvectors(dm)[, 2],
                      DC3 = eigenvectors(dm)[, 3],
                      DC4 = eigenvectors(dm)[, 4],
                      Samples = cellLabels)
    pdf("./DC1_DC2.pdf", w=11, h=8.5)
    ggplot(tmp, aes(x = DC1, y = DC2, colour = Samples)) +
      geom_point()  + 
      xlab("Diffusion Component 1") + 
      ylab("Diffusion Component 2") +
      theme_classic()
    dev.off()

    # Plotting cell progression along the diffusion map components
    sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # rank cells by their dpt dm1
    sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # rank cells by their dpt dm2
    sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # rank cells by their dpt dm1 reverse order
    sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # rank cells by their dpt dm2 reverse order

    SortedDM1 <- data.frame(DM1Sort = as.data.frame(colData(sce))$pseud_dm1,
                            Samples = as.data.frame(colData(sce))$ident)
    SortedDM2 <- data.frame(DM2Sort = as.data.frame(colData(sce))$pseud_dm2,
                            Samples = as.data.frame(colData(sce))$ident)
    SortedDM1R <- data.frame(DM1SortR = as.data.frame(colData(sce))$pseud_dm1R,
                            Samples = as.data.frame(colData(sce))$ident)
    SortedDM2R <- data.frame(DM2SortR = as.data.frame(colData(sce))$pseud_dm2R,
                            Samples = as.data.frame(colData(sce))$ident)

    ggplot(SortedDM1, aes(x=SortedDM1[,1], y=Samples,color=Samples)) +
      geom_jitter() + xlab("Diffusion Component 1 (DC1)") + ylab("Samples") +
      ggtitle("Cells ordered by DC1")
    ggplot(SortedDM2, aes(x=SortedDM2[,1], y=Samples,color=Samples)) +
      geom_jitter() + xlab("Diffusion Component 2 (DC2)") + ylab("Samples") +
      ggtitle("Cells ordered by DC2")
      
    ggplot(SortedDM1R, aes(x=SortedDM1R[,1], y=Samples,color=Samples)) +
      geom_jitter() + xlab("Minus Diffusion Component 1 (DC1)") + ylab("Samples") +
      ggtitle("Cells ordered by reversed DC1")

    ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Samples,color=Samples)) +
      geom_jitter() + xlab("Minus Diffusion Component 2 (DC2)") + ylab("Samples") +
      ggtitle("Cells ordered by reversed DC2")

}

##############################################################################################################
# here according to https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html#comparison-of-the-different-trajectory-inference-methods
# use log-transformed counts as input

destiny_logcounts <- function(sce,outDir,n_pcs=50){

    dir.create(outDir)
      
    logCounts <- logcounts(sce)  # access log-transformed counts matrix
    colnames(logCounts) <- colData(sce)$seurat_clusters

    # Make a diffusion map
    dm  <- DiffusionMap(as.matrix(t(logCounts)),n_pcs = n_pcs)

    # Plot eigenvalues of diffusion distance matrix. How many diffusion components would you use?
    # This is analagous to the PC elbow plot (scree plot) that we previously used to assess how 
    # many PCs to use in downstream applications like clustering. 
    pdf(paste0(outDir,"/eigenvalues.pdf"),5,5)
    plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
    dev.off()

    # calculate Diffusion pseudotime
    dpt <- DPT(dm)

    # Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
    df <- data.frame( DC1 = eigenvectors(dm)[, 1],
                      DC2 = eigenvectors(dm)[, 2],
                      Cell_type = colData(sce)$seurat_clusters)

    ggplot(df, aes(x = DC1, y = DC2, colour = Cell_type)) +
        geom_point() + #scale_color_tableau() + 
        xlab("Diffusion Component 1") + 
        ylab("Diffusion Component 2") +
        theme_classic() +
        ggtitle("Diffusion map")

    ggsave(paste0(outDir,"/DC1_DC2.pdf"),width=6,height=5)

    # Next, let us use the first diffusion component (DC1) as a measure of pseudotime.
    sce$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt

    ggplot(as.data.frame(colData(sce)), 
          aes(x = pseudotime_diffusionmap, 
              y = seurat_clusters, colour = seurat_clusters)) +
        geom_jitter() +
        theme_classic() +
        xlab("Diffusion Component 1 (DC1)") + ylab("Cell Type") +
        ggtitle("Cells ordered by DC1")

    ggsave(paste0(outDir,"/diffusion_pseudotime.pdf"),width=6,height=5)

    save(sce,dm,dpt,file=paste0(outDir,"/destiny.Rdata"))

}

##############################################################################################################
# use PCA as input
#  here according to https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html#comparison-of-the-different-trajectory-inference-methods
destiny_pca <- function(obj,outDir,reverse_heatmap=FALSE){

    dir.create(outDir)
      
    #pca <- reducedDim(sce, "PCA")
    pca <- Embeddings(obj,'pca')

    rownames(pca) <- obj$seurat_clusters

    dm <- DiffusionMap(pca)

    dpt <- DPT(dm)

    # Plot DC1 vs DC2 and color the cells by their inferred diffusion pseudotime.
    df <- data.frame(DC1 = eigenvectors(dpt@dm)[, 1], 
                     DC2 = eigenvectors(dpt@dm)[, 2], 
                     PC1 = pca[,1],
                     PC2 = pca[,2],
                     seurat_clusters = obj$seurat_clusters
                     )

    ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = seurat_clusters)) + 
        xlab("Diffusion Component 1") + 
        ylab("Diffusion Component 2") +
        theme_classic() +
        ggtitle("Diffusion map")

    ggsave(paste0(outDir,"/DC1_DC2.pdf"),width=6,height=5)

    # Next, let us use the first diffusion component (DC1) as a measure of pseudotime.
    obj$destiny_DC1 <- eigenvectors(dpt@dm)[,1]
    obj$destiny_DC2 <- eigenvectors(dpt@dm)[,2]
    obj$PC1 <- pca[,1]
    obj$PC2 <- pca[,2]

    obj$pseudotime_diffusionmap <- rank(eigenvectors(dpt@dm)[,1]) 

    ggplot(obj@meta.data, 
          aes(x = pseudotime_diffusionmap, 
              y = seurat_clusters, colour = seurat_clusters)) +
        geom_jitter() +
        theme_classic() +
        xlab("Pseudotime") +
        ylab("Clusters") +
        ggtitle("Cells ordered by pseudotime")
    ggsave(paste0(outDir,"/diffusion_pseudotime.pdf"),width=6,height=5)


    # Next, let us use the first diffusion component (DC1) as a measure of pseudotime, but reverse the rank
    obj$pseudotime_diffusionmap_reverse <- rank(-eigenvectors(dpt@dm)[,1]) 

    ggplot(obj@meta.data, 
          aes(x = pseudotime_diffusionmap_reverse, 
              y = seurat_clusters, colour = seurat_clusters)) +
        geom_jitter() +
        theme_classic() +
        xlab("Pseudotime") +
        ylab("Clusters") +
        ggtitle("Cells ordered by pseudotime")
    ggsave(paste0(outDir,"/diffusion_pseudotime.reverse.pdf"),width=6,height=5)

    # For PCA plot, color cells by pseudotime
    colors <- rainbow(10, alpha = 1)
    ggplot(obj@meta.data,aes(x = PC1,y = PC2, colour = pseudotime_diffusionmap)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2") + scale_color_gradientn(colours = colors)
    ggsave(paste0(outDir,"/PC1_PC2_pseudotime.pdf"),width=7.5,height=5)

    # For PCA plot, color cells by reversed pseudotime
    colors <- rainbow(10, alpha = 1)
    ggplot(obj@meta.data,aes(x = PC1,y = PC2, colour = pseudotime_diffusionmap_reverse)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2") + scale_color_gradientn(colours = colors)
    ggsave(paste0(outDir,"/PC1_PC2_pseudotime_reverse.pdf"),width=7.5,height=5)

    # let's try to add DC1 and DC2
    obj$pseudotime_diffusionmap_DC1_DC2 <- rank(eigenvectors(dpt@dm)[,1] + eigenvectors(dpt@dm)[,2]) 

    ggplot(obj@meta.data, 
          aes(x = pseudotime_diffusionmap_DC1_DC2, 
              y = seurat_clusters, colour = seurat_clusters)) +
        geom_jitter() +
        theme_classic() +
        xlab("DC1 + DC2") +
        ylab("Clusters") +
        ggtitle("Cells ordered by DC1 + DC2")
    ggsave(paste0(outDir,"/diffusion_pseudotime_DC1_DC2.pdf"),width=6,height=5)

    # let's try to add all DCs
    obj$pseudotime_diffusionmap_DCs <- rank(rowSums(eigenvectors(dpt@dm))) 

    ggplot(obj@meta.data, 
          aes(x = pseudotime_diffusionmap_DCs, 
              y = seurat_clusters, colour = seurat_clusters)) +
        geom_jitter() +
        theme_classic() +
        xlab("DCs") +
        ylab("Clusters") +
        ggtitle("Cells ordered by DCs")
    ggsave(paste0(outDir,"/diffusion_pseudotime_DCs.pdf"),width=6,height=5)


    # find temporally expressed genes
    # Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
    # Identify the variable genes by ranking all genes by their variance.
    find_temporally_expressed_genes <- function(obj,DC1,reverse_heatmap=FALSE,outDir){

            #Y <- log2(counts(sce) + 1)
            obj <- NormalizeData(obj,assay='RNA')
            Y <- obj$RNA@data
            nexpressed_cells = apply(Y,1,function(x) sum(x>0))
            Y <- Y[nexpressed_cells>=100,]; dim(Y)

            #var1K <- names(sort(apply(Y, 1, var), decreasing = TRUE))[1:1000]
            #Y <- Y[var1K, ]  # only counts for variable genes

            HVGs      <- VariableFeatures(obj[['SCT']])
            HVGs      <- head(HVGs[HVGs %in% rownames(Y)], 1000)
            Y <- Y[HVGs,]
            
            # Fit GAM for each gene using pseudotime as independent variable.
            library(gam)
            t <- DC1
            gam.pval <- apply(Y, 1, function(z){
                    d <- data.frame(z=z, t=t)
                    tmp <- gam(z ~ lo(t), data=d)
                    p <- summary(tmp)[4][[1]][1,5]
                    p
            })

            write.table(data.frame(gene=rownames(Y),gam_pval=gam.pval),paste0(outDir, "/time_genes.tsv"), quote=F,sep="\t",row.names=F)

            gam.pval <- gam.pval[gam.pval<0.05]; length(gam.pval)

            # Identify genes with the most significant time-dependent model fit.
            topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  

            # Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
            require(clusterExperiment)
            if(reverse_heatmap){
                heatdata <- as.matrix(Y[rownames(Y) %in% topgenes, order(t, na.last = NA,decreasing=TRUE)])
                heatclus <- Idents(obj)[order(t, na.last = NA,decreasing=TRUE)]
            }else{
                heatdata <- as.matrix(Y[rownames(Y) %in% topgenes, order(t, na.last = NA)])
                heatclus <- Idents(obj)[order(t, na.last = NA)]
            }
            

            pdf(paste0(outDir, "/heatmap_time_genes.top100_timeDep_genes.ClusterExperiment.pdf"), width=10, height=10)
            ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
            print(clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'original', cexRow = 2.5, fontsize = 15))
            dev.off()

            # using complexheatmap package
            library(ComplexHeatmap)

            cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                                    "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                                    "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                                    "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4")
            clusters <- sort(unique(heatclus))
            cluster_cols             = cluster_cols[1:length(clusters)]
            names(cluster_cols) <- clusters

            cols <- list('Cluster'=cluster_cols)

            col_ha = HeatmapAnnotation(
                                        'Cluster'=heatclus,
                                        col=cols, 
                                        show_annotation_name=TRUE,annotation_name_side = 'left')

            p <- Heatmap(heatdata,name="Log Counts",column_title=paste0(ncol(heatdata)," cells"),top_annotation = col_ha,cluster_rows = T, cluster_columns = F,show_row_dend=F,show_row_names = T ,show_column_names = F,row_names_side='left',row_names_gp = gpar(fontsize = 5),,use_raster = TRUE,raster_device="CairoPNG",gap = unit(1, "mm"),row_title_gp = gpar(fontsize=10))

            pdf(paste0(outDir, "/heatmap_time_genes.top100_timeDep_genes.complexheatmap.pdf"), width=10, height=10); print(p);dev.off()

    }

    find_temporally_expressed_genes(obj,DC1=obj$destiny_DC1,reverse_heatmap=reverse_heatmap,outDir=outDir)

    # 
    save(obj,dm,dpt,file=paste0(outDir,"/destiny.Rdata"))

}

slingshot_pca <- function(obj,outDir,defaultAssay='SCT'){

    dir.create(outDir)

    DefaultAssay(obj) <- defaultAssay

    pca <- Embeddings(obj,'pca')

    sce <- as.SingleCellExperiment(obj)
    colData(sce)$seurat_clusters <- obj$seurat_clusters
    sce <- slingshot(sce,clusterLabels='seurat_clusters',reducedDim='PCA')
    obj$slingshot_PC1 <- pca[,1]
    obj$slingshot_PC2 <- pca[,2]
    obj$slingPseudotime_1 <- sce$slingPseudotime_1
    obj$slingPseudotime_1_rank <- rank(sce$slingPseudotime_1)

    ggplot(obj@meta.data,aes(x = slingshot_PC1,y=slingshot_PC2, colour = seurat_clusters)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2")
    ggsave(paste0(outDir,"/PC1_PC2_clusters.pdf"),width=6,height=5)

    colors <- rainbow(10, alpha = 1)
    ggplot(obj@meta.data,aes(x = slingshot_PC1,y=slingshot_PC2, colour = slingPseudotime_1)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2") + scale_color_gradientn(colours = colors) 
    ggsave(paste0(outDir,"/PC1_PC2_pseudotime.pdf"),width=6,height=5)

    ggplot(obj@meta.data,aes(x = slingPseudotime_1_rank,y=seurat_clusters, colour = seurat_clusters)) + geom_jitter() + theme_classic() + xlab("Pseudotime") + ylab("Clusters")
    ggsave(paste0(outDir,"/slingshot_pseudotime.pdf"),width=6,height=5)

    save(sce,obj,file=paste0(outDir,"/slingshot.Rdata"))

    # find temporally expressed genes
    # Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
    # Identify the variable genes by ranking all genes by their variance.
    find_temporally_expressed_genes <- function(){
    
        #Y <- log2(counts(sce) + 1)
        obj <- NormalizeData(obj,assay='RNA')
        Y <- obj$RNA@data
        nexpressed_cells = apply(Y,1,function(x) sum(x>0))
        Y <- Y[nexpressed_cells>=100,]; dim(Y)

        #var1K <- names(sort(apply(Y, 1, var), decreasing = TRUE))[1:1000]
        #Y <- Y[var1K, ]  # only counts for variable genes

        HVGs      <- VariableFeatures(obj)
        HVGs      <- head(HVGs[HVGs %in% rownames(Y)], 1000)
        Y <- Y[HVGs,]
        

        # Fit GAM for each gene using pseudotime as independent variable.
        library(gam)
        t <- obj$slingPseudotime_1
        gam.pval <- apply(Y, 1, function(z){
                d <- data.frame(z=z, t=t)
                tmp <- gam(z ~ lo(t), data=d)
                p <- summary(tmp)[4][[1]][1,5]
                p
        })

        gam.pval <- gam.pval[gam.pval<0.05]; length(gam.pval)

        # Identify genes with the most significant time-dependent model fit.
        topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  
        # Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
        require(clusterExperiment)
        heatdata <- as.matrix(Y[rownames(Y) %in% topgenes, order(t, na.last = NA)])
        heatclus <- Idents(obj)[order(t, na.last = NA)]

        pdf(paste0(outDir, "/heatmap_time_genes.top100_timeDep_genes.pdf"), width=20, height=10)
        ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
        clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'original', cexRow = 1.5, fontsize = 15)
        #heatmap(log1p(heatdata), Colv = NA, ColSideColors = brewer.pal(9,"Set1")[heatclus])
        dev.off()


    }

}



slingshot_dm <- function(obj,outDir){

    dir.create(outDir)

    pca <- Embeddings(obj,'pca')

    sce <- as.SingleCellExperiment(obj)
    colData(sce)$seurat_clusters <- obj$seurat_clusters
    sce <- slingshot(sce,clusterLabels='seurat_clusters',reducedDim='DM')
    obj$slingshot_PC1 <- dm@eigenvectors[,1]
    obj$slingshot_PC2 <- dm@eigenvectors[,2]
    obj$slingPseudotime_1 <- sce$slingPseudotime_1
    obj$slingPseudotime_1_rank <- rank(sce$slingPseudotime_1)

    ggplot(obj@meta.data,aes(x = slingshot_PC1,y=slingshot_PC2, colour = seurat_clusters)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2")
    ggsave(paste0(outDir,"/PC1_PC2_clusters.pdf"),width=6,height=5)

    colors <- rainbow(10, alpha = 1)
    ggplot(obj@meta.data,aes(x = slingshot_PC1,y=slingshot_PC2, colour = slingPseudotime_1)) + geom_jitter() + theme_classic() + xlab("PC1") + ylab("PC2") + scale_color_gradientn(colours = colors) 
    ggsave(paste0(outDir,"/PC1_PC2_pseudotime.pdf"),width=6,height=5)

    ggplot(obj@meta.data,aes(x = slingPseudotime_1_rank,y=seurat_clusters, colour = seurat_clusters)) + geom_jitter() + theme_classic() + xlab("Pseudotime") + ylab("Clusters")
    ggsave(paste0(outDir,"/slingshot_pseudotime.pdf"),width=6,height=5)

}
