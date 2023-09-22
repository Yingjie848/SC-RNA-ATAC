
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)

# set default colors for ggplot, which are used for cell clusters
cluster_cols <- c("#023FA5","#A30059","#008941","#7D87B9",  "#D6BCC0", "#EF9708","#E07B91","#8DD593", "#11C638","#D33F6A","#8595E1","gold",  "#F79CD4", "#663399","#706563",
                  "#BB7784",  "#4A6FE3", "#E6AFB9", "#0FCFC0", "#8B0000","#008080",
                  "#B5BBE3",  "#C6DEC7", "#EAD3C6", "#F0B98D",  
                  "#9CDED6",  "#F6C4E1", "#866097", "#F79CD4", "#BEC1D4",
                  "#66bdf0","#7ae5a7","#7adee5","#e57ae0","#dc4a8e","#ddd259","#7140e1","black","cyan")
options(ggplot2.discrete.colour = cluster_cols)

rename_clusters <- function(obj,named_clusters=NULL){

  # MBC6
  #named_clusters <- c(
  #  "0"="Epithelial cells",
  #  "1"="Epithelial cells",
  #  "2"="Epithelial cells",
  #  "3"="Epithelial cells",
  #  "4"="Epithelial cells",
  #  "5"="Macrophages",
  #  "6"="Epithelial cells",
  #  "7"="T cells",
  #  "8"="Fibroblast cells",
  #  "9"="Endothelial cells",
  #  "10"="Pericytes"
  #)
  
  obj@meta.data$cell_identity <- named_clusters[as.character(obj@meta.data$seurat_clusters)]

  obj

}

add_cell_identity <- function(obj,named_clusters=NULL,outDir="annotated_cells"){

    dir.create(outDir)
    
    obj <- rename_clusters(obj,named_clusters)

    p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_identity", label = TRUE,repel = TRUE) 
    
    pdf(paste0(outDir,"/DimPlot_grouped_by_cell_identity.pdf"),6.5,5);print(p1);dev.off()

    saveRDS(obj,file=paste0(outDir,"/seurat_rna.annotated_cells.rds"))

    # make cell components plot
    metaData <- obj@meta.data
    cells_per_population = metaData %>% group_by(cell_identity) %>% summarise(cells=length(cell_identity)) %>% mutate(cell_pct=cells/nrow(metaData)*100) %>% arrange(desc(cell_pct))
    cells_per_population$cell_identity <- factor(cells_per_population$cell_identity,levels=cells_per_population$cell_identity)

    p <- cells_per_population %>% ggplot(aes(x=cell_identity,y=cell_pct)) + 
                                geom_bar(stat='identity',position=position_dodge2(),fill='royalblue') + theme_bw() +
                                xlab("") + ylab("Cell Components") + ggtitle("Cell Components (%)") + 
                                geom_text(aes(label=round(cell_pct,1)), vjust= -0.4,color='black') + 
                                theme(axis.text.x=element_text(angle=45,size=8,vjust=1,hjust=1,color='black'))

    pdf(paste0(outDir,"/cells_per_cluster_with_annotated_cells.pdf"),5,5); print(p); dev.off()

    obj

}

