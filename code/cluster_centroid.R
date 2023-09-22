
calculate_cluster_centroid <- function(obj){

    metaData <- obj@meta.data
    clusters = sort(unique(metaData$seurat_clusters))
    cluster_centroid = list()

    for(cluster in clusters){
      
        id=cluster
        print(id)
        cells <- rownames(metaData[metaData$seurat_clusters==cluster,,drop=F])
        cluster_mean <- rowMeans(obj$RNA@data[,cells,drop=F])
        cluster_centroid[[id]] <- cluster_mean
      
    }
    cluster_centroid <- do.call(cbind,cluster_centroid)

}

make_centroid_correlation_heatmap <- function(obj,outDir){

    library(ComplexHeatmap)

    dir.create(outDir)

    cluster_centroid <- calculate_cluster_centroid(obj)

    # exclude genes all centroids are 0
    excluded <- apply(cluster_centroid,1,function(x) sum(x)==0)
    cluster_centroid_nonzero <- cluster_centroid[!excluded,]

    centroid_cor <- cor(cluster_centroid_nonzero)

    h <- Heatmap(centroid_cor,name="Corr.")

    pdf(paste0(outDir,"/centroid_correlation_heatmap.pdf"));print(h);dev.off()

    save(cluster_centroid,centroid_cor,file=paste0(outDir,"/centroid.Rdata"))

}