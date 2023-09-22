# annotate cells
library(dplyr)
library(Seurat)

annotate_cells <- function(obj){
    
    if(class(obj)!='Seurat' & is.character(obj)){
      if(file.exists(obj))
        obj <- readRDS(obj)
    }
    
    cat("Number of cells: ",nrow(obj@meta.data),"\n")
    cat("Number of genes: ",nrow(obj),"\n")

    cell_info <- obj@meta.data %>% mutate(cell=rownames(obj@meta.data))

    # remove previous 'PanCK' and 'CD45' columns
    cell_info <- cell_info[,! colnames(cell_info) %in% c('PanCK','CD45')]

    RNA_data  <- as.matrix(obj$RNA@counts)

    # Add Pan-Cytokeratin (Pan-CK) counts 
    print("Add Pan-Cytokeratin (Pan-CK) counts")
    panck_genes <- rownames(RNA_data)[grepl('^KRT',rownames(RNA_data))]
    panck_genes <- panck_genes[!grepl('CAP',panck_genes)]
    #print(panck_genes)

    RNA_data_panck <- RNA_data[panck_genes,]

    df_panck  <- data.frame(cell=colnames(RNA_data_panck),PanCK=colSums(RNA_data_panck))
    cell_info <- merge(cell_info,df_panck,by='cell',all.x=T,sort=F)

    # Add CD45 counts
    print("Add CD45 counts")
    CD45_expr   <- t(RNA_data)[,'PTPRC']
    CD45_counts <- data.frame(cell=names(CD45_expr),CD45=CD45_expr)

    cell_info   <- merge(cell_info,CD45_counts,by='cell',all.x=T,sort=F) 


    # order cells by original matrix
    obj@meta.data <- cell_info[match(colnames(RNA_data),cell_info$cell),]
    rownames(obj@meta.data) <- cell_info$cell

    obj

}