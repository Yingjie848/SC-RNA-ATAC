


        #########################################################################################################
        # // determine immune markers if the cluster expressed more than 50% CD45

        if(is.null(immune_markers)){

            cell_info   = obj@meta.data %>% mutate(cell=rownames(obj@meta.data),
                                                   CD45_status=ifelse(CD45>=1,'CD45+','CD45-')) %>% 
                              arrange(as.numeric(as.vector(seurat_clusters)))

            cell_cluster_CD45_freq <- cell_info %>% 
                                          group_by(seurat_clusters) %>% 
                                          summarise(ncells=length(cell),CD45_frequency=sum(CD45>=1)) %>% 
                                          mutate(CD45_cell_pct=CD45_frequency/ncells*100) %>% 
                                          arrange(desc(CD45_cell_pct)) # order cell clusters by CD45+ frequency

            write.table(cell_cluster_CD45_freq,paste0(outDir,"/immune_clusters.cell_cluster_CD45_freq.txt"),quote=F,sep="\t",row.names=F,col.names=F)

            immune_clusters = as.vector(cell_cluster_CD45_freq[cell_cluster_CD45_freq$CD45_cell_pct>50,]$seurat_clusters)

            write.table(immune_clusters,paste0(outDir,"/immune_clusters.txt"),quote=F,row.names=F,col.names=F)


            immune_markers <- unique(markers_pos[markers_pos$cluster %in% immune_clusters,]$gene)

            write.table(immune_markers,paste0(outDir,"/DE_analysis/immune_markers.txt"),quote=F,row.names=F,col.names=F)

        }


        #########################################################################################################
        # caluclate immune GSEA score for each cell using immune markers
        if(length(immune_markers)>=5){

            print("Calculating immune_gsea score")

            library(GSVA)
            library(GSEABase)

            gmt = c('Immune_markers','Immune_markers',immune_markers)
            gmt = toString(paste(gmt,sep=","))
            gmt = gsub(" ","",gmt)
            gmt = gsub(",","\t",gmt)
            gmt = paste0(gmt,"\n")
            gmt.file = paste0(outDir,"/DE_analysis/immune_markers.gmt")
            cat(gmt,file=gmt.file)

            genesets    <- getGmt(gmt.file)
            gsva_output <- gsva(as.matrix(obj$SCT@data),genesets,method="ssgsea")[1,]

            obj@meta.data$immune_gsea <- gsva_output[match(rownames(obj@meta.data),names(gsva_output))]
            write.table(obj@meta.data,paste0(outDir,"/metaData.immune_gsea.txt"),quote=F,sep="\t",row.names=T)

        }else{

            obj@meta.data$immune_gsea <- 0

        }
        