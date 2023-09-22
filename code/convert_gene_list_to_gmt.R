
convert_gene_list_to_gmt <- function(genelist,geneset.name,gmt.file){
        
        genes <- read.table(genelist)$V1
        gmt = c(geneset.name,geneset.name,genes)
        gmt = toString(paste(gmt,sep=","))
        gmt = gsub(" ","",gmt)
        gmt = gsub(",","\t",gmt)
        gmt = paste0(gmt,"\n")
        cat(gmt,file=gmt.file)

}