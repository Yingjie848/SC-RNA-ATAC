
library(Seurat)
library(corrplot)


make_cluster_cluster_correlationPlot <- function(ref_obj,test_obj,genes,outDir="outdir",ref_clusters=NULL,test_clusters=NULL,ref_name='R',test_name='T'){

        dir.create(outDir,recursive = T)

        # extract gene matrix
        ref_mat  <- t(t(as.matrix(ref_obj$RNA@data))[,genes])
        test_mat <- t(t(as.matrix(test_obj$RNA@data))[,genes])

        # make cluster matrix according to seurat_clusters columns in meta.data
        ref_metaData  <- ref_obj@meta.data
        test_metaData <- test_obj@meta.data
 
        if(is.null(ref_clusters)){
            ref_clusters <- sort(unique(as.numeric(as.vector(ref_metaData$seurat_clusters))))
        }
        ref_clust_mat <- lapply(ref_clusters,function(cluster){
            clust_cells <- rownames(ref_metaData[ref_metaData$seurat_clusters==cluster,])
            mat <- ref_mat[,clust_cells]
            rowMeans(mat)
        })
        ref_clust_mat <- do.call(cbind,ref_clust_mat)
        colnames(ref_clust_mat) <- paste0(ref_name,'-C',ref_clusters)

        if(is.null(test_clusters)){
            test_clusters <- sort(unique(as.numeric(as.vector(test_metaData$seurat_clusters))))
        }
        test_clust_mat <- lapply(test_clusters,function(cluster){
            clust_cells <- rownames(test_metaData[test_metaData$seurat_clusters==cluster,])
            mat <- test_mat[,clust_cells]
            rowMeans(mat)
        })
        test_clust_mat <- do.call(cbind,test_clust_mat)
        colnames(test_clust_mat) <- paste0(test_name,'-C',test_clusters)

        # make ref-ref plot

        res <- cor(ref_clust_mat,ref_clust_mat)

        pdf(paste0(outDir,"/ref_ref.pdf"))
        corrplot(res, type = "upper", order = "original", tl.col = "black", tl.srt = 45)
        dev.off()

        # make test-test plot

        res <- cor(test_clust_mat,test_clust_mat)

        pdf(paste0(outDir,"/test_test.pdf"))
        corrplot(res, type = "upper", order = "original", tl.col = "black", tl.srt = 45)
        dev.off()

        # make test-ref plot

        res <- cor(test_clust_mat,ref_clust_mat)
        pdf(paste0(outDir,"/test_ref.pdf"))
        corrplot(res, type = "full", order = "original", tl.col = "black", tl.srt = 45)
        dev.off()

        save(ref_metaData,test_metaData,ref_mat,test_mat,ref_clust_mat,test_clust_mat,file=paste0(outDir,"/data.Rdata"))

}

