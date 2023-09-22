library(clusterProfiler)
library(org.Hs.eg.db)
library(magrittr)
library(ReactomePA)
library(stringr)
library(enrichplot)
library(pathview)
library(ggplot2)

# filter ORA results
filter_ORA <- function(x){
  x.result <- x@result
  x.incutoff <- x.result[x.result$p.adjust <= x@pvalueCutoff,]
  x.incutoff <- x.incutoff[!is.na(x.incutoff$ID),]
  x@result <- x.incutoff
  return(x)
}

# filter GSEA results
filter_GSEA <- function(x){
  x.result <- x@result
  x.incutoff <- x.result[x.result$p.adjust <= x@params$pvalueCutoff,]
  x.incutoff <- x.incutoff[!is.na(x.incutoff$ID),]
  x@result <- x.incutoff
  return(x)
}

# Core wrapping function
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

sigSymbol <- function(pvals){
  pvals <- as.numeric(as.vector(pvals))
  pvals.symbol <- ifelse(pvals<=0.0001,"****",
                  ifelse(pvals<=0.001,"***",
                  ifelse(pvals<=0.01,"**",
                  ifelse(pvals<=0.05,"*","ns"))))
  return(pvals.symbol)
}

my.barplot <- function(x,title,n=10,font.size=12,clean_desc=TRUE){
  
  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]
  
  genes <- unique(unlist(strsplit(as.vector(result$geneID),"/")))
  print("Genes:")
  print(length(genes))
  
  n=10
  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]
  #print(result)
  
  result <- result[order(result[,"p.adjust"],decreasing = T),]
  
  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                               gsub("^KEGG ","",
                                    gsub("^REACTOME ","",
                                         result$Description)))
  }
  
  result$Description <- wrap.labels(result$Description,35)
  
  result$Description <- factor(result$Description,levels=unique(result$Description))
  
  #print(result)
  
  # plot
  options(scipen = 1)
  options(digits = 2)
  
  ymax = max(result$Count) * 1.2
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')
  
  p <- ggplot(data=result,aes(x=Description,y=Count,fill=sigSymbol(p.adjust),width=0.6)) + geom_bar(stat="identity") +
    xlab("") + 
    ylab("Gene count") +
    ylim(0,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("Adjusted P-value",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=45, size=font.size), 
      axis.text.y=element_text(color="black",size=font.size), 
      text = element_text(size=font.size),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    ) 
  return(p)
}

my.barplot.gene_ratio <- function(x,title,n=10,font.size=12,clean_desc=TRUE){
  
  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]
  
  genes <- unique(unlist(strsplit(as.vector(result$geneID),"/")))
  print("Genes:")
  print(length(genes))
  
  n=10
  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]
  #print(result)
  
  result <- result[order(result[,"p.adjust"],decreasing = T),]
  
  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                               gsub("^KEGG ","",
                                    gsub("^REACTOME ","",
                                         result$Description)))
  }
  
  result$Description <- wrap.labels(result$Description,35)
  
  result$Description <- factor(result$Description,levels=unique(result$Description))
  
  result$Gene1 <- unlist(sapply(as.vector(result$GeneRatio),function(x) as.numeric(strsplit(x,"/")[[1]][1])))
  result$Gene2 <- unlist(sapply(as.vector(result$BgRatio),function(x) as.numeric(strsplit(x,"/")[[1]][1])))
  result$Gene_ratio <- result$Gene1/result$Gene2
  
  #print(result)
  
  # plot
  options(scipen = 1)
  options(digits = 2)
  
  ymax = max(result$Gene_ratio) * 1.2
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')
  
  p <- ggplot(data=result,aes(x=Description,y=Gene_ratio,fill=sigSymbol(p.adjust),label=Count,width=0.6)) + geom_bar(stat="identity") + geom_text(hjust=-0.1) +
    xlab("") + 
    ylab("Gene Ratio") +
    ylim(0,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("Adjusted P-value",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=font.size), 
      axis.text.y=element_text(color="black",size=font.size), 
      text = element_text(size=font.size),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    ) 
  return(p)
}


my.barplot.gene_ratio.label_pval <- function(x,title,n=10,font.size=12,clean_desc=TRUE){

  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]

  genes <- unique(unlist(strsplit(as.vector(result$geneID),"/")))
  print("Genes:")
  print(length(genes))

  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]

  result <- result[order(result[,"p.adjust"],decreasing = T),]

  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                          gsub("^KEGG ","",
                          gsub("^REACTOME ","",
                          gsub("_"," ",
                          result$Description))))
  }

  result$Description <- wrap.labels(result$Description,30)

  result$Description <- factor(result$Description,levels=unique(result$Description))

  result$Gene1 <- unlist(sapply(as.vector(result$GeneRatio),function(x) as.numeric(strsplit(x,"/")[[1]][1])))
  result$Gene2 <- unlist(sapply(as.vector(result$BgRatio),function(x) as.numeric(strsplit(x,"/")[[1]][1])))
  result$Gene_ratio <- result$Gene1/result$Gene2

  #print(result)

  # plot
  options(scipen = 1)
  options(digits = 2)
  library(scales)

  ymax = max(result$Gene_ratio) * 1.8
  print(paste("ymax:",max(result$Gene_ratio)))
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')

  # bquote(~italic(P)~"="~.(scientific(p.adjust,digits=2)))
  p <- ggplot(data=result,aes(x=Description,y=Gene_ratio,fill=sigSymbol(p.adjust),label=paste0('P = ',scientific(p.adjust,digit=3)," (",Count,")"),width=0.6)) + geom_bar(stat="identity") + geom_text(hjust=-0.1,size=3) +
    xlab("") +
    ylab("Gene Ratio") +
    ylim(0,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    # Add 10% spaces between the p-value labels and the plot border
    #scale_y_continuous(expand = expand_scale(mult = 0.3)) +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("FDR",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=font.size),
      axis.text.y=element_text(color="black",size=font.size),
      text = element_text(size=font.size),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"),
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    )
  return(p)
}


my.barplot.gsea <- function(x,title,n=10,clean_desc=TRUE){
  
  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]
  
  genes <- unique(unlist(strsplit(as.vector(result$core_enrichment),"/")))
  print("Genes:")
  print(length(genes))
  
  n=10
  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]
  
  result <- result[order(result[,"p.adjust"],decreasing = T),]
  
  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                               gsub("^KEGG ","",
                                    gsub("^REACTOME ","",
                                         result$Description)))
  }
  
  result$Description <- wrap.labels(result$Description,35)
  
  result$Description <- factor(result$Description,levels=result$Description)
  
  #print(result)
  
  # plot
  options(scipen = 1)
  options(digits = 2)
  
  ymax = max(result$Gene_ratio) * 1.2
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')
  
  p <- ggplot(data=result,aes(x=Description,y=setSize,fill=sigSymbol(p.adjust),width=0.6)) + geom_bar(stat="identity") +
    xlab("") + 
    ylab("Gene count") +
    ylim(0,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("Adjusted P-value",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=12), 
      axis.text.y=element_text(color="black",size=12), 
      text = element_text(size=12),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    ) 
  return(p)
}

my.barplot.gsea.enrichment_score <- function(x,title,n=10,clean_desc=TRUE){

  cat("Top ",n,"\n")
  
  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]
  
  genes <- unique(unlist(strsplit(as.vector(result$core_enrichment),"/")))
  print("Genes:")
  print(length(genes))
  
  #n=10
  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]
  cat("Top ",nrow(result),"\n")
  
  result <- result[order(result[,"p.adjust"],decreasing = T),]
  
  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                               gsub("^KEGG ","",
                                    gsub("^REACTOME ","",
                                         result$Description))) 
  }
  
  result$Description <- wrap.labels(result$Description,35)
  
  result$Description <- factor(result$Description,levels=result$Description)
  
  #print(result)
  
  # plot
  options(scipen = 1)
  options(digits = 2)
  
  print(class(result$NES))
  print(head(result$NES))
  print(head(as.numeric(as.vector(result$NES))))
  result$NES <- as.numeric(as.vector(result$NES))
  ymin = min(as.numeric(as.vector(result$NES))) * 1.2
  ymax = max(as.numeric(as.vector(result$NES))) * 1.2
  if( ymin > 0 )
    ymin = 0
  if( ymax < 0 )
    ymax = 0
  print(c("my.barplot.gsea.enrichment_score: ",ymin,ymax))
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')
  
  p <- ggplot(data=result,aes(x=Description,y=NES,fill=sigSymbol(p.adjust),width=0.6)) + geom_bar(stat="identity") +
    xlab("") + ylab("Normalized enrichment score") +
    ylim(ymin,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("Adjusted P-value",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=10), 
      axis.text.y=element_text(color="black",size=6), 
      text = element_text(size=12),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    ) 
  return(p)
}


my.barplot.gsea.enrichment_score.facets <- function(x,title,n=10,clean_desc=TRUE){

  cat("Top ",n,"\n")
  
  #result <- x@result
  result <- x
  result <- result[order(result[,"p.adjust"]),]
  
  genes <- unique(unlist(strsplit(as.vector(result$core_enrichment),"/")))
  print("Genes:")
  print(length(genes))
  
  #n=10
  if( nrow(result) < n )
    n = nrow(result)
  result <- result[1:n,]
  cat("Top ",nrow(result),"\n")
  
  result <- result[order(result[,"p.adjust"],decreasing = T),]
  
  if( isTRUE(clean_desc) ){
    result$Description <- gsub("^GO ","",
                               gsub("^KEGG ","",
                                    gsub("^REACTOME ","",
                                         result$Description))) 
  }
  
  result$Description <- wrap.labels(result$Description,35)
  
  result$Description <- factor(result$Description,levels=result$Description)
  
  #print(result)
  
  # plot
  options(scipen = 1)
  options(digits = 2)
  
  print(class(result$NES))
  print(head(result$NES))
  print(head(as.numeric(as.vector(result$NES))))
  result$NES <- as.numeric(as.vector(result$NES))
  ymin = min(as.numeric(as.vector(result$NES))) * 1.2
  ymax = max(as.numeric(as.vector(result$NES))) * 1.2
  if( ymin > 0 )
    ymin = 0
  if( ymax < 0 )
    ymax = 0
  print(c("my.barplot.gsea.enrichment_score: ",ymin,ymax))
  col_fun <- c('ns'='#0473B5','*'='#107C81','**'='#1C864E','***'='#7E863A','****'='#E28627')
  
  p <- ggplot(data=result,aes(x=Description,y=NES,fill=sigSymbol(p.adjust),width=0.6)) + geom_bar(stat="identity") +
    xlab("") + ylab("Normalized enrichment score") +
    ylim(ymin,ymax) +
    #labs(fill="p.adjust") +
    ggtitle(title) +
    coord_flip() +
    theme_bw() +
    #scale_fill_gradient(name="p.adjust",low="red",high="blue",guide = guide_colourbar(reverse=T)) +
    scale_fill_manual("Adjusted P-value",values=col_fun) +
    theme(
      axis.text.x=element_text(color="black", hjust=0.5, vjust=0.5, lineheight=1, angle=0, size=10), 
      axis.text.y=element_text(color="black",size=6), 
      text = element_text(size=12),
      axis.line.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.line.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.line.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      axis.ticks.x = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.x.top = element_line(size = 0.5, linetype = "solid",colour = "black"),
      axis.ticks.y.right = element_line(size = 0.5, linetype = "solid",colour = "black"), 
      #axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
      legend.position="right",
      #legend.title=element_blank(),
      plot.title=element_text(size=12,margin=margin(0,20,20,20)), # adjust the distance between label and axis
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
    ) + 
    facet_wrap(~cluster)
  return(p)
}

plot_ORA <- function(file_DEres,file_ORA,title=""){
  
  geneList <- get_geneList_limma(file_DEres)
  
  ORAres <- read.table(file_ORA,header=T,sep="\t")
  
  #setClass("gseaResult", representation(result = "data.frame", geneList = "numeric",keytype="character",
  #                                      readable="character"))
  #oo <- new("gseaResult",result=ORAres,geneList=geneList,keytype="ENTREZID",readable="TRUE")
  #oo <- setReadable(oo,org.Hs.eg.db,keyType = "ENTREZID")
  
  pdf(paste(file_ORA,"barplot.pdf",sep="."))
  print(my.barplot.gene_ratio.label_pval(ORAres,title=title,n=10,font.size = 14))
  #print(enrichplot::cnetplot(oo, categorySize="pvalue", foldChange=geneList))
  #print(enrichplot::emapplot(oo,showCategory = 10))
  #print(enrichplot::heatplot(oo,foldChange = geneList,showCategory = 10))
  #print(enrichplot::upsetplot(oo))
  dev.off()
  
}


exclude_pathways <- function(x){
  x.result <- x@result
  x.incutoff <- x.result
  x.incutoff <- x.incutoff[!grepl("Leishmaniasis",x.incutoff$Description), ]
  x.incutoff <- x.incutoff[!grepl("Hematopoietic cell lineage",x.incutoff$Description), ]
  x.incutoff <- x.incutoff[!grepl("Asthma",x.incutoff$Description), ]
  x@result <- x.incutoff
  return(x)
}

exclude_pathways_in_MSigDB <- function(x){
  x.result <- x@result
  x.incutoff <- x.result
  x.incutoff <- x.incutoff[!grepl("GSE",x.incutoff$Description), ]
  x.incutoff <- x.incutoff[!grepl("MODULE",x.incutoff$Description), ]
  x.incutoff <- x.incutoff[!grepl("MIR",x.incutoff$Description), ]
  x.incutoff <- x.incutoff[!grepl("chr",x.incutoff$Description), ]
  x@result <- x.incutoff
  return(x)
}

# group GO
GO_group <- function(geneID,input){
  print("Group GO")
  ggo_CC <- groupGO(gene     = geneID,
                    OrgDb    = org.Hs.eg.db,
                    ont      = "CC",
                    level    = 3,
                    readable = TRUE)
  
  ggo_CC <- ggo_CC[ggo_CC$Count>0,]
  print(nrow(ggo_CC))
  write.table(ggo_CC,paste(input,"DEgenes_ggo_CC.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  
  ggo_MF <- groupGO(gene     = geneID,
                    OrgDb    = org.Hs.eg.db,
                    ont      = "MF",
                    level    = 3,
                    readable = TRUE)
  
  ggo_MF <- ggo_MF[ggo_MF$Count>0,]
  print(nrow(ggo_MF))
  write.table(ggo_MF,paste(input,"DEgenes_ggo_MF.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  
  ggo_BP <- groupGO(gene     = geneID,
                    OrgDb    = org.Hs.eg.db,
                    ont      = "BP",
                    level    = 3,
                    readable = TRUE)
  
  ggo_BP <- ggo_BP[ggo_BP$Count>0,]
  print(nrow(ggo_BP))
  write.table(ggo_BP,paste(input,"DEgenes_ggo_BP.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
}

# GO over-representative test
GO_over_representative_test <- function(geneID,geneList,input,useUniverse=FALSE,use_geneList=FALSE){
  print("GO over-representative test")
  if(isTRUE(useUniverse)){
    ego_MF <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }else{
    ego_MF <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       #universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }
  #write.table(ego_MF,paste(input,"DEgenes_ego_MF.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_MF,paste(input,"DEgenes_ego_MF.csv",sep="."),row.names=F)
  ego_MF <- filter_ORA(ego_MF)
  #write.table(ego_MF,paste(input,"DEgenes_ego_MF.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_MF,paste(input,"DEgenes_ego_MF.incutoff.csv",sep="."),row.names=F)

  if(isTRUE(useUniverse)){
    ego_CC <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }else{
    ego_CC <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       #universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }
  #write.table(ego_CC,paste(input,"DEgenes_ego_CC.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_CC,paste(input,"DEgenes_ego_CC.csv",sep="."),row.names=F)
  ego_CC <- filter_ORA(ego_CC)
  #write.table(ego_CC,paste(input,"DEgenes_ego_CC.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_CC,paste(input,"DEgenes_ego_CC.incutoff.csv",sep="."),row.names=F)
  
  if(isTRUE(useUniverse)){
    ego_BP <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }else{
    ego_BP <- enrichGO(gene = geneID,
                       org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       #universe = names(geneList),
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       minGSSize     = 10,
                       readable      = TRUE)
  }
  #write.table(ego_BP,paste(input,"DEgenes_ego_BP.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_BP,paste(input,"DEgenes_ego_BP.csv",sep="."),row.names=F)
  ego_BP <- filter_ORA(ego_BP)
  #write.table(ego_BP,paste(input,"DEgenes_ego_BP.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(ego_BP,paste(input,"DEgenes_ego_BP.incutoff.csv",sep="."),row.names=F)
  
  if( nrow(ego_MF@result) > 0 ){
    pdf(paste(input,"DEgenes_ego_MF.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(ego_MF,title="GO: Molecular Function",n=10,font.size=14))
    #print(barplot(ego_MF,showCategory = 10,color='p.adjust',font.size=14,title="GO: Molecular Function"))
    #if( isTRUE(use_geneList) ){
    #  print(cnetplot(ego_MF, categorySize="pvalue", foldChange=geneList))
    #}else{ print(cnetplot(ego_MF, categorySize="pvalue")) }
    #print(emapplot(ego_MF,showCategory = 10))
    #if( isTRUE(use_geneList) ){
    #  print(heatplot(ego_MF,foldChange = geneList,showCategory = 10))
    #}else{ print(heatplot(ego_MF,showCategory = 10)) }
    #print(upsetplot(ego_MF))
    dev.off()
  }
  if( nrow(ego_CC@result) > 0 ){
    pdf(paste(input,"DEgenes_ego_CC.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(ego_CC,title="GO: Cellular Component",n=10,font.size=14))
    #print(barplot(ego_CC,showCategory = 10,color='p.adjust',font.size=14,title="GO: Cellular Component"))
    #if( isTRUE(use_geneList) ){
    #  print(cnetplot(ego_CC, categorySize="pvalue", foldChange=geneList))
    #}else{ print(cnetplot(ego_CC, categorySize="pvalue")) }
    #print(emapplot(ego_CC,showCategory = 10))
    #if( isTRUE(use_geneList) ){
    #  print(heatplot(ego_CC,foldChange = geneList,showCategory = 10))
    #}else{ print(heatplot(ego_CC,showCategory = 10)) }
    #print(upsetplot(ego_CC))
    dev.off()
  }
  if( nrow(ego_BP@result) > 0 ){
    pdf(paste(input,"DEgenes_ego_BP.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(ego_BP,title="GO: Biological Process",n=10,font.size=14))
    #print(barplot(ego_BP,showCategory = 10,color='p.adjust',font.size=14,title="GO: Biological Process"))
    #if( isTRUE(use_geneList) ){
    #  print(cnetplot(ego_BP, categorySize="pvalue", foldChange=geneList))
    #}else{ print(cnetplot(ego_BP, categorySize="pvalue")) }
    #print(emapplot(ego_BP,showCategory = 10))
    #if( isTRUE(use_geneList) ){
    #  print(heatplot(ego_BP,foldChange = geneList,showCategory = 10))
    #}else{ print(heatplot(ego_BP,showCategory = 10))  }
    #print(upsetplot(ego_BP))
    dev.off()
  }
}

# GO GSEA
GO_GSEA <- function(geneList,input){
  
  print("GO GSEA")
  gseago_MF <- gseGO(geneList   = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "MF",
                     keyType      = "ENTREZID", 
                     nPerm        = 10000,
                     minGSSize    = 10,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose      = TRUE)
  gseago_MF <- setReadable(gseago_MF,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(gseago_MF,paste(input,"DEgenes_gseago_MF.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_MF,paste(input,"DEgenes_gseago_MF.csv",sep="."))
  gseago_MF <- filter_GSEA(gseago_MF)
  #write.table(gseago_MF,paste(input,"DEgenes_gseago_MF.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_MF,paste(input,"DEgenes_gseago_MF.incutoff.csv",sep="."))
  
  gseago_CC <- gseGO(geneList   = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "CC",
                     keyType      = "ENTREZID", 
                     nPerm        = 10000,
                     minGSSize    = 10,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose      = TRUE)
  gseago_CC <- setReadable(gseago_CC,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(gseago_CC,paste(input,"DEgenes_gseago_CC.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_CC,paste(input,"DEgenes_gseago_CC.csv",sep="."),row.names=F)
  gseago_CC <- filter_GSEA(gseago_CC)
  #write.table(gseago_CC,paste(input,"DEgenes_gseago_CC.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_CC,paste(input,"DEgenes_gseago_CC.incutoff.csv",sep="."),row.names=F)
  
  gseago_BP <- gseGO(geneList   = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     keyType      = "ENTREZID", 
                     nPerm        = 10000,
                     minGSSize    = 10,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose      = TRUE)
  gseago_BP <- setReadable(gseago_BP,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(gseago_BP,paste(input,"DEgenes_gseago_BP.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_BP,paste(input,"DEgenes_gseago_BP.csv",sep="."),row.names=F)
  gseago_BP <- filter_GSEA(gseago_BP)
  #write.table(gseago_BP,paste(input,"DEgenes_gseago_BP.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(gseago_BP,paste(input,"DEgenes_gseago_BP.incutoff.csv",sep="."),row.names=F)
  
  #print(gseago_MF@result)
  print(nrow(gseago_MF@result))
  if( nrow(gseago_MF@result) > 0 ){
    pdf(paste(input,"DEgenes_gseago_MF.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(gseago_MF,title="GO: Molecular Function",n=10))
    #print(barplot(gseago_MF,showCategory = 10,color='p.adjust',font.size=14,title="GO: Molecular Function"))
    #print(cnetplot(gseago_MF, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(gseago_MF,showCategory = 10))
    #print(heatplot(gseago_MF,foldChange = geneList,showCategory = 10))
    #print(upsetplot(gseago_MF))
    nmax=5
    if( nrow(gseago_MF@result) < 5 )
      nmax = nrow(gseago_MF@result)
    gseago_MF@result$Description <- wrap.labels(gseago_MF@result$Description,35)
    print(gseaplot2(gseago_MF,geneSetID = 1:nmax))
    dev.off()
  }
  if( nrow(gseago_CC@result) > 0 ){
    pdf(paste(input,"DEgenes_gseago_CC.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(gseago_CC,title="GO: Cellular Component",n=10))
    #print(barplot(gseago_CC,showCategory = 10,color='p.adjust',font.size=14,title="GO: Cellular Component"))
    #print(cnetplot(gseago_CC, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(gseago_CC,showCategory = 10))
    #print(heatplot(gseago_CC,foldChange = geneList,showCategory = 10))
    #print(upsetplot(gseago_CC))
    nmax=5
    if( nrow(gseago_CC@result) < 5 )
      nmax = nrow(gseago_CC@result)
    gseago_CC@result$Description <- wrap.labels(gseago_CC@result$Description,35)
    print(gseaplot2(gseago_CC,geneSetID = 1:nmax))
    dev.off()
  }
  if( nrow(gseago_BP@result) > 0 ){
    pdf(paste(input,"DEgenes_gseago_BP.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(gseago_BP,title="GO: Biological Process",n=10))
    #print(barplot(gseago_BP,showCategory = 10,color='p.adjust',font.size=14,title="GO: Biological Process"))
    #print(cnetplot(gseago_BP, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(gseago_BP,showCategory = 10))
    #print(heatplot(gseago_BP,foldChange = geneList,showCategory = 10))
    #print(upsetplot(gseago_BP))
    nmax=5
    if( nrow(gseago_BP@result) < 5 )
      nmax = nrow(gseago_BP@result)
    gseago_BP@result$Description <- wrap.labels(gseago_BP@result$Description,35)
    print(gseaplot2(gseago_BP,geneSetID = 1:nmax))
    dev.off()
  }
}

# run KEGG over-representative test
KEGG_over_representative_test <- function(geneID,geneList,input,useUniverse=FALSE,use_geneList=TRUE){
  print("KEGG over-representative test")
  if(isTRUE(useUniverse)){
    kk <- enrichKEGG(gene         = geneID,
                     organism     = 'hsa',
                     keyType      = 'ncbi-geneid',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     universe = names(geneList),
                     minGSSize    = 10,
    )
  }else{
    kk <- enrichKEGG(gene         = geneID,
                     organism     = 'hsa',
                     keyType      = 'ncbi-geneid',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     #universe = names(geneList),
                     minGSSize    = 10,
    )
  }
  kk <- setReadable(kk,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(kk,paste(input,"DEgenes_ekegg.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(kk,paste(input,"DEgenes_ekegg.csv",sep="."),row.names=F)
  kk <- filter_ORA(kk)
  kk <- exclude_pathways(kk)
  #write.table(kk,paste(input,"DEgenes_ekegg.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(kk,paste(input,"DEgenes_ekegg.incutoff.csv",sep="."),row.names=F)
  
  if( nrow(kk@result) > 0 ){
    pdf(paste(input,"DEgenes_ekegg.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(kk,title="KEGG ORA",n=10,font.size=14))
    #print(barplot(kk,showCategory = 10,color='p.adjust',font.size=14,title="KEGG over-representative test"))
    #if( isTRUE(use_geneList) ){
    #  print(cnetplot(kk, categorySize="pvalue", foldChange=geneList))
    #}else{ print(cnetplot(kk, categorySize="pvalue")) }
    #print(emapplot(kk,showCategory = 10))
    #if( isTRUE(use_geneList) ){
    #  print(heatplot(kk,foldChange = geneList,showCategory = 10))
    #}else{ print(heatplot(kk,showCategory = 10)) }
    #print(upsetplot(kk))
    dev.off()
    
  }
}

# run KEGG GSEA
KEGG_GSEA <- function(geneList,input){
  print("KEGG GSEA")
  kk2 <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 10000,
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose      = TRUE)
  kk2 <- setReadable(kk2,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(kk2,paste(input,"DEgenes_gsea_kegg.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(kk2,paste(input,"DEgenes_gsea_kegg.csv",sep="."),row.names=F)
  kk2 <- filter_GSEA(kk2)
  #write.table(kk2,paste(input,"DEgenes_gsea_kegg.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(kk2,paste(input,"DEgenes_gsea_kegg.incutoff.csv",sep="."),row.names=F)
  
  if( nrow(kk2@result) > 0 ){
    pdf(paste(input,"DEgenes_gsea_kegg.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(kk2,title="KEGG GSEA",n=10))
    #print(barplot(kk2,showCategory = 10,color='p.adjust',font.size=14,title="KEGG GSEA"))
    #print(cnetplot(kk2, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(kk2,showCategory = 10))
    #print(heatplot(kk2,foldChange = geneList,showCategory = 10))
    #print(upsetplot(kk2,n=10))
    nmax=5
    if( nrow(kk2@result) < 5 )
      nmax = nrow(kk2@result)
    kk2@result$Description <- wrap.labels(kk2@result$Description,35)
    print(gseaplot2(kk2,geneSetID = 1:nmax))
    dev.off()
  }
}  

KEGG_pathview <- function(geneList,input){
  
  hsa04110 <- pathview(gene.data = geneList, pathway.id = "hsa04110", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
  # pathways in cancer
  hsa05200 <- pathview(gene.data = geneList, pathway.id = "hsa05200", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
  # breast cancer
  hsa05224 <- pathview(gene.data = geneList, pathway.id = "hsa05224", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
  # hsa04512:ECM-receptor interaction
  hsa04512 <- pathview(gene.data = geneList, pathway.id = "hsa04512", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
  # hsa04510:Focal adhesion
  hsa04510 <- pathview(gene.data = geneList, pathway.id = "hsa04510", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
  # hsa04514:Cell adhesion molecules (CAMs)
  hsa04514 <- pathview(gene.data = geneList, pathway.id = "hsa04514", species    = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1), gene.idtype = "entrez",
                       out.suffix=paste(input,"kegg",sep="."))
}

# haven't done
WikiPathway_analysis <- function(geneID,geneList,input){
  print("run wikiPathway enricher")
  wpgmtfile <- system.file("/Users/zhuy1/Google_Drive/my_projects_MSKCC/TNBC_NanoString/TNBC_data_analysis/wikipathways-20200110-gpml-Homo_sapiens.zip", package="clusterProfiler")
  wp2gene <- read.gmt(wpgmtfile)
  wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  ewp <- enricher(geneID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  print(head(ewp))
  
  write.table(ewp,paste(input,"DEgenes_wikipath_enricher.txt",sep="."),col.names=T,row.names=F,quote=F)
  
  
  print("run wikiPathway GSEA")
  ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,  minGSSize = 1, pvalueCutoff = 0.05, verbose=FALSE)
  print(head(ewp2))
  
  write.table(ewp2,paste(input,"DEgenes_wikipath_GSEA.txt",sep="."),col.names=T,row.names=F,quote=F)
  
}

# 
MSigDb_ORA <- function(geneID,geneList,input){
  print("run MSigDb ORA analysis")
  gmtfile <- "/home/zhuy1/my_databases/MSigDB/msigdb.v7.1.entrez.gmt"
  gmt <- read.gmt(gmtfile)
  
  egmt <- enricher(gene         = geneID,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize    = 10,
                   pAdjustMethod = "BH",
                   TERM2GENE=gmt)
  #print(head(egmt))
  egmt <- exclude_pathways_in_MSigDB(egmt)
  egmt <- setReadable(egmt,org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(egmt,paste(input,"MSigDb_enricher.csv",sep="."),row.names=F)
  
  if( nrow(egmt@result) > 0 ){
    pdf(paste(input,"MSigDb_enricher.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(egmt,title="MSigDB ORA",n=10,font.size=12))
    #print(barplot(egmt,showCategory = 10,color='p.adjust',font.size=14,title="MSigDB enrichment pathway analysis"))
    #print(cnetplot(egmt, categorySize="pvalue", foldChange=geneList,showCategory = 10))
    #print(emapplot(egmt,showCategory = 10))
    #print(heatplot(egmt,foldChange = geneList,showCategory = 10))
    #print(upsetplot(egmt))
    dev.off()
  }
  
}

MSigDb_GSEA <- function(geneID,geneList,input){
  print("run MSigDb GSEA analysis")
  gmtfile <- "/home/zhuy1/my_databases/MSigDB/msigdb.v7.1.entrez.gmt"
  gmt <- read.gmt(gmtfile)
  
  ggmt <- GSEA(geneList, TERM2GENE=gmt, nPerm = 10000, minGSSize = 10, pvalueCutoff = 0.05, verbose=FALSE)
  ggmt <- setReadable(ggmt,org.Hs.eg.db,keyType = "ENTREZID")
  
  write.csv(ggmt,paste(input,"MSigDb_GSEA.csv",sep="."),col.names=T,row.names=F,quote=F)
  
}

MSigDb_ORA_ARsignatures <- function(geneID,geneList,input){
  print("run MSigDb analysis for AR signatures")
  gmtfile <- "~/my_projects/LAR_project/MSigDB_Androgen/merged_AR_genesets.entrez.gmt"
  gmt <- read.gmt(gmtfile)
  
  egmt <- enricher(gene=geneID,
                   pvalueCutoff=1,
                   qvalueCutoff = 1,
                   minGSSize = 1,
                   pAdjustMethod = "BH",
                   TERM2GENE=gmt)
  #print(head(egmt))
  
  if( is.null(egmt) ){
    return(0)
  }
  
  #egmt <- exclude_pathways_in_MSigDB(egmt)
  egmt <- setReadable(egmt,org.Hs.eg.db,keyType = "ENTREZID")
  
  write.csv(egmt,paste(input,"MSigDb_ARsignatures_enricher.csv",sep="."))

  if( nrow(egmt@result) > 0 ){
    egmt <- filter_ORA(egmt)
  }
  
  if( nrow(egmt@result) > 0 ){
    print(paste("results",nrow(egmt@result)))
    #print(egmt@result)
    pdf(paste(input,"MSigDb_enricher.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(egmt,title="MSigDB ORA",n=10,font.size=14))
    #print(barplot(egmt,showCategory = 10,color='p.adjust',font.size=14,title="MSigDB enrichment pathway analysis"))
    print(cnetplot(egmt, categorySize="pvalue", foldChange=geneList,showCategory = 10))
    print(emapplot(egmt,showCategory = 10))
    print(heatplot(egmt,foldChange = geneList,showCategory = 10))
    print(upsetplot(egmt))
    dev.off()
  }
  
}

Reactome_enrichPathway <- function(geneID,geneList,input,useUniverse=FALSE,use_geneList=TRUE){
  
  print("Reactome enriched pathway analysis")
  if(isTRUE(useUniverse)){
    x <- enrichPathway(gene=geneID,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       minGSSize = 10,
                       organism = "human",
                       universe = names(geneList),
                       pAdjustMethod = "BH",
                       readable=T)
  }else{
    x <- enrichPathway(gene=geneID,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       minGSSize = 10,
                       organism = "human",
                       #universe = names(geneList),
                       pAdjustMethod = "BH",
                       readable=T)
  }
  
  #write.table(x,paste(input,"DEgenes_Reactome_enrich.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(x,paste(input,"DEgenes_Reactome_enrich.csv",sep="."),row.names=F)
  
  x <- filter_ORA(x)
  #write.table(x,paste(input,"DEgenes_Reactome_enrich.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(x,paste(input,"DEgenes_Reactome_enrich.incutoff.csv",sep="."),row.names=F)
  
  if( nrow(x@result) > 0 ){
    pdf(paste(input,"DEgenes_Reactome_enrich.pdf",sep="."))
    print(my.barplot.gene_ratio.label_pval(x,title="Reactome ORA",n=10,font.size=14))
    #print(barplot(x,showCategory = 10,color='p.adjust',font.size=14,title="Reactome enrichment pathway analysis"))
    #if( isTRUE(use_geneList) ){
    #  print(cnetplot(x, categorySize="pvalue", foldChange=geneList,showCategory = 10))
    #}else{ print(cnetplot(x, categorySize="pvalue",showCategory = 10)) }
    #print(emapplot(x,showCategory = 10))
    #if( isTRUE(use_geneList) ){
    #  print(heatplot(x,foldChange = geneList,showCategory = 10))
    #}else{ print(heatplot(x,showCategory = 10)) }
    #print(upsetplot(x))
    dev.off()
  }
}

Reactome_GSEA <- function(geneList,input){
  
  print("Reactome GSEA analysis")
  y <- gsePathway(gene  = geneList, 
                  nPerm = 10000,
                  organism = "human",
                  pvalueCutoff = 1,
                  by = "fgsea",
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  verbose = TRUE)
  y <- setReadable(y,org.Hs.eg.db,keyType = "ENTREZID")
  #write.table(y,paste(input,"DEgenes_Reactome_GSEA.txt",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
  write.csv(y,paste(input,"DEgenes_Reactome_GSEA.csv",sep="."),row.names=F)
  y <- filter_GSEA(y)
  #write.table(y,paste(input,"DEgenes_Reactome_GSEA.incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  write.csv(y,paste(input,"DEgenes_Reactome_GSEA.incutoff.csv",sep="."),row.names=F)
  
  
  if( nrow(y@result)>0 ){
    pdf(paste(input,"DEgenes_Reactome_GSEA.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(y,title="Reactome GSEA",n=10))
    #print(barplot(y,showCategory=10, color='p.adjust',font.size=14))
    #print(cnetplot(y, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(y,showCategory = 10))
    #print(heatplot(y,foldChange = geneList,showCategory = 10))
    #print(upsetplot(y,n=10))
    nmax=6
    if( nrow(y@result) < 5 )
      nmax = nrow(y@result)
    y@result$Description <- wrap.labels(y@result$Description,35)
    print(gseaplot2(y,geneSetID = 1:nmax))
    dev.off()
  }
  
}

# cell cycle related
# R-HSA-69620: Cell Cycle Checkpoints
# R-HSA-69278: Cell Cycle, Mitotic
# R-HSA-69273: Cyclin A/B1/B2 associated events during G2/M transition
# R-HSA-69275: G2/M Transition
# R-HSA-69206: G1/S Transition
# R-HSA-69242: S Phase
# R-HSA-68886: M Phase

# R-HSA-2467813: Separation of Sister Chromatids
# R-HSA-453274: Mitotic G2-G2/M phases
# R-HSA-69618: Mitotic Spindle Checkpoint
# R-HSA-2500257: Resolution of Sister Chromatid Cohesion
# R-HSA-68877: Mitotic Prometaphase
# R-HSA-2555396: Mitotic Metaphase and Anaphase
# R-HSA-68882: Mitotic Anaphase

# R-HSA-1474244: Extracellular matrix organization

# R-HSA-5685942: HDR through Homologous Recombination (HRR)

# R-HSA-5693571: Nonhomologous End-Joining (NHEJ)

Reactome_GSEA_select_pathways <- function(geneList,input){
  
  # cell cycle 1
  pathways1=c("R-HSA-69620","R-HSA-69278","R-HSA-69273","R-HSA-69275","R-HSA-69206","R-HSA-69242","R-HSA-68886")
  # cell cycle 2
  pathways2=c("R-HSA-2467813","R-HSA-453274","R-HSA-69618","R-HSA-68877","R-HSA-2555396","R-HSA-68882")
  # R-HSA-5685942: HDR through Homologous Recombination
  pathways3=c("R-HSA-5685942")
  # R-HSA-5693571: Nonhomologous End-Joining (NHEJ)
  pathways4=c("R-HSA-5693571")
  
  pathway_name = "selected_pathways"
  
  print("Reactome GSEA analysis")
  y <- gsePathway(gene  = geneList, 
                  nPerm = 10000,
                  organism = "human",
                  pvalueCutoff = 1,
                  by = "fgsea",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  verbose = FALSE)
  y <- setReadable(y,org.Hs.eg.db,keyType = "ENTREZID")
  write.table(y,paste(input,"DEgenes_Reactome_GSEA",pathway_name,"txt",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
  #y <- filter_GSEA(y)
  write.table(y,paste(input,"DEgenes_Reactome_GSEA",pathway_name,"incutoff.txt",sep="."),col.names=T,row.names=F,quote=F,sep="\t")
  
  
  if( nrow(y@result)>0 ){
    pdf(paste(input,"DEgenes_Reactome_GSEA",pathway_name,"pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(y,title="Reactome GSEA",n=10))
    #print(barplot(y,showCategory=10, color='p.adjust',font.size=14))
    print(cnetplot(y, categorySize="pvalue", foldChange=geneList))
    print(emapplot(y,showCategory = 10))
    print(heatplot(y,foldChange = geneList,showCategory = 10))
    #print(upsetplot(y,n=10))
    nmax=6
    if( nrow(y@result) < 5 )
      nmax = nrow(y@result)
    
    y@result$Description <- wrap.labels(y@result$Description,35)
    print(gseaplot2(y,geneSetID = 1:nmax))
    
    selected_pathways <- match(pathways1,y@result$ID)
    selected_pathways <- selected_pathways[!is.na(selected_pathways)]
    if( length(selected_pathways) > 0 )
      print(gseaplot2(y,geneSetID = selected_pathways,title="Cell Cycle",pvalue_table = FALSE,subplots=c(1,3),base_size = 15))
    
    selected_pathways <- match(pathways2,y@result$ID)
    selected_pathways <- selected_pathways[!is.na(selected_pathways)]
    if( length(selected_pathways) > 0 )
      print(gseaplot2(y,geneSetID = selected_pathways,title="Cell Cycle",pvalue_table = FALSE,subplots=c(1,3),base_size = 15))
    
    selected_pathways <- match(pathways3,y@result$ID)
    selected_pathways <- selected_pathways[!is.na(selected_pathways)]
    if( length(selected_pathways) > 0 )
      print(gseaplot2(y,geneSetID = selected_pathways,title="HDR through Homologous Recombination",pvalue_table = FALSE,subplots=c(1,2,3),base_size = 15))
    
    selected_pathways <- match(pathways4,y@result$ID)
    selected_pathways <- selected_pathways[!is.na(selected_pathways)]
    if( length(selected_pathways) > 0 )
      print(gseaplot2(y,geneSetID = selected_pathways,title="Nonhomologous End-Joining (NHEJ)",pvalue_table = FALSE,subplots=c(1,2,3),base_size = 15))
    
    dev.off()
  }
  
}

# clean title as output name
clean_output_name <- function(x){
  x <- gsub(" ","_",x)
  x <- gsub("\\/","_",x)
  x <- gsub("ER-","ERneg",x)
  x <- gsub("HER2-","HER2neg",x)
  x <- gsub("\\+","pos",x)
  return(x)
}

#gmtfile <- "~/my_projects/LAR_project/MSigDB_Androgen/merged_AR_genesets.entrez.gmt"
#gmtfile_AR <- "~/my_projects/LAR_project/MSigDB_Androgen/DOANE_NELSON.entrez.gmt"
#gmtfile_AR <- "~/my_projects/LAR_project/MSigDB_Androgen/DOANE.entrez.gmt"
#gmtfile_FA <- "/Users/zhuy1/my_projects/LAR_project/MSigDB_fatty_acid_metabolism_pathway/KEGG_FATTY_ACID_METABOLISM.genelist.entrezid.gmt"

do_GSEA_on_gmt <- function(gmtfile,geneList,title,input){
  
  print(paste("run GSEA for",title))
  
  gmt <- read.gmt(gmtfile)
  
  write.csv(geneList,paste(input,clean_output_name(title),"genelist.csv",sep="."))
  
  ggmt <- GSEA(geneList, TERM2GENE=gmt, nPerm = 10000, 
               minGSSize = 10, pvalueCutoff = 1, verbose=FALSE)
  
  ggmt <- setReadable(ggmt,org.Hs.eg.db,keyType = "ENTREZID")
  
  write.csv(ggmt,paste(input,clean_output_name(title),"GSEA.csv",sep="."),row.names=F)
  
  if( nrow(ggmt@result)>0 ){
    pdf(paste(input,clean_output_name(title),"GSEA.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(ggmt,title="MSigDb ARsignatures GSEA",n=10))
    #print(barplot(ggmt,showCategory=10, color='p.adjust',font.size=14))
    #print(cnetplot(ggmt, categorySize="pvalue", foldChange=geneList))
    #print(emapplot(ggmt,showCategory = 10))
    #print(heatplot(ggmt,foldChange = geneList,showCategory = 10))
    #print(upsetplot(ggmt,n=10))
    nmax=6
    if( nrow(ggmt@result) < nmax )
      nmax = nrow(ggmt@result)
    ggmt@result$Description <- wrap.labels(ggmt@result$Description,30)
    print(gseaplot2(ggmt,geneSetID = 1:nmax,title = title,pvalue_table = TRUE,subplots=c(1,2,3),base_size = 14))
    dev.off()
  }
}

get_geneList <- function(input){
  
  # load output of DE analysis
  data <- read.csv(input)
  
  #data <- data[order(data$pval_t.test),] 
  
  data <- data[data$gene!="HBD" & data$gene!="MEMO1",] # there are two EntrezID for HBD and MEMO1
  data <- data[data$gene!="MMD2" & data$gene!="TEC",] # there are two EntrezID for MMD2 and TEC
  
  genemap <- bitr(data$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  
  data$EntrezID <- genemap$ENTREZID
  
  data <- data[!is.na(data$EntrezID),]
  
  geneList <- data$log2fc
  names(geneList) <- data$EntrezID
  geneList <- sort(geneList,decreasing = TRUE)
  
  return(geneList)
  
}

get_geneList_limma <- function(input){
  
  # load output of DE analysis
  data <- read.csv(input)
  prefix=input
  
  data <- data[order(data$adj.P.Val),] 
  
  data <- data[data$X!="HBD" & data$X!="MEMO1",] # there are two EntrezID for HBD and MEMO1
  data <- data[data$X!="MMD2" & data$X!="TEC",] # there are two EntrezID for MMD2 and TEC
  
  genemap <- bitr(data$X,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  
  data$EntrezID <- genemap$ENTREZID
  
  data <- data[!is.na(data$EntrezID),]
  
  write.table(data$EntrezID,paste(prefix,"ENTREZID.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  
  geneList <- data$logFC
  names(geneList) <- data$EntrezID
  geneList <- sort(geneList,decreasing = TRUE)
  
  return(geneList)
  
}

run_clusterProfiler.limma <- function(input,prefix,fc_cutoff,pval_cutoff,outdir="DEanalysis"){
  
  print("run_clusterProfiler.limma")
  
  dir.create(outdir)
  prefix=paste(outdir,prefix,sep="/")
  
  # load output of DE analysis
  data <- read.csv(input,stringsAsFactors=FALSE)
  data$logFC <- as.numeric(as.vector(data$logFC))
  
  data <- data[order(data$adj.P.Val),] 

  if(!'entrez_id' %in% colnames(data)){
    data <- data[data$X!="HBD" & data$X!="MEMO1",] # there are two EntrezID for HBD and MEMO1
    data <- data[data$X!="MMD2" & data$X!="TEC",] # there are two EntrezID for MMD2 and TEC
    genemap <- bitr(data$X,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
    data$EntrezID <- genemap$ENTREZID
  }else{
    data$EntrezID <- data$entrez_id
  }
  
  data <- data[!is.na(data$EntrezID),]
  write.table(data$EntrezID,paste(prefix,"ENTREZID.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  
  geneList <- data$logFC
  names(geneList) <- data$EntrezID
  geneList <- sort(geneList,decreasing = TRUE)
  write.table(geneList,paste(prefix,"ENTREZID_FC.txt",sep="."),quote=FALSE,row.names = T,col.names = F)
  
  # get DE genes
  data.DE      <- subset(data,abs(logFC) > log2(fc_cutoff) & adj.P.Val < pval_cutoff)
  data.DE.up   <- subset(data,logFC > log2(fc_cutoff) & adj.P.Val < pval_cutoff)
  data.DE.down <- subset(data,logFC < -log2(fc_cutoff) & adj.P.Val < pval_cutoff)
  
  print("Differential expressed genes")
  print(nrow(data.DE))
  print(nrow(data.DE.up))
  print(nrow(data.DE.down))
  
  geneID      <- data.DE$EntrezID
  geneID.up   <- data.DE.up$EntrezID
  geneID.down <- data.DE.down$EntrezID
  
  write.table(geneID,paste(prefix,"DEgenes.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  
  # for all genes
  outprefix = prefix
  
  print(length(geneList))
  if( length(geneList) > 5 ){
    #GO_group(geneID,outprefix)
    GO_over_representative_test(geneID,geneList,outprefix)
    #GO_GSEA(geneList,outprefix)
    
    KEGG_over_representative_test(geneID,geneList,outprefix)
    KEGG_GSEA(geneList,outprefix)
    #KEGG_pathview(geneList,outprefix)
    
    #WikiPathway_analysis(geneID,geneList,outprefix)
    
    MSigDb_ORA(geneID,geneList,outprefix)
    #MSigDb_ORA_ARsignatures(geneID,geneList,outprefix)
    #MSigDb_GSEA(geneID,geneList,outprefix)
    Reactome_enrichPathway(geneID,geneList,outprefix)
    Reactome_GSEA(geneList,outprefix)
  }
  
  # for up-regulated genes
  outprefix = paste(prefix,"up",sep=".")
  write.table(geneID.up,paste(outprefix,"DEgenes.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  if( length(geneID.up) > 5 ){
    GO_over_representative_test(geneID.up,geneList,outprefix)
    KEGG_over_representative_test(geneID.up,geneList,outprefix)
    MSigDb_ORA(geneID.up,geneList,outprefix)
    #MSigDb_ORA_ARsignatures(geneID.up,geneList,outprefix)
    Reactome_enrichPathway(geneID.up,geneList,outprefix)
  }
  
  # for down-regulated genes
  outprefix = paste(prefix,"down",sep=".")
  write.table(geneID.down,paste(outprefix,"DEgenes.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  if( length(geneID.down) > 5 ){
    GO_over_representative_test(geneID.down,geneList,outprefix)
    KEGG_over_representative_test(geneID.down,geneList,outprefix)
    MSigDb_ORA(geneID.down,geneList,outprefix)
    #MSigDb_ORA_ARsignatures(geneID.down,geneList,outprefix)
    Reactome_enrichPathway(geneID.down,geneList,outprefix)
  }
}

pathwayAnnotator <- function(genes,outfile,showCategory=6,pval_thres=1,qval_thres=1,title=""){
  
  print(outfile)
  
  genemap <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  
  geneID      <- genemap$ENTREZID
  
  enrichRes <- enrichPathway(gene=geneID,
                pvalueCutoff=pval_thres,
                qvalueCutoff = qval_thres,
                minGSSize = 1,
                organism = "human",
                pAdjustMethod = "BH",
                readable=T)
  
  if(!is.null(enrichRes)){
    enrichRes@result$Description <- wrap.labels(enrichRes@result$Description,35)
    enrichRes <- filter_ORA(enrichRes)
  
  if( nrow(enrichRes@result)>0 ){
    
  # get first hit of genes
  hits <- sapply(unique(enrichRes@result$geneID),function(x) grep(x,enrichRes@result$geneID)[1])
    
  enrichRes@result <- unique(enrichRes@result[hits,])
    
  grDevices::pdf(outfile,7,5)
  print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=14))
  print(cnetplot(enrichRes,showCategory = showCategory))
  print(emapplot(enrichRes,showCategory = showCategory))
  print(heatplot(enrichRes,showCategory = showCategory))
  print(upsetplot(enrichRes,n=showCategory))
  grDevices::dev.off()
  
  grDevices::pdf(paste(outfile,"pdf",sep="."),6,4)
  print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=14))
  print(cnetplot(enrichRes,showCategory = showCategory))
  print(emapplot(enrichRes,showCategory = showCategory))
  print(heatplot(enrichRes,showCategory = showCategory))
  print(upsetplot(enrichRes,n=showCategory))
  grDevices::dev.off()
  
  grDevices::pdf(paste(outfile,"pdf","pdf",sep="."),10,10)
  print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=14))
  print(cnetplot(enrichRes,showCategory = showCategory))
  print(emapplot(enrichRes,showCategory = showCategory))
  print(heatplot(enrichRes,showCategory = showCategory))
  print(upsetplot(enrichRes,n=showCategory))
  grDevices::dev.off()
  
  }
    
    return(enrichRes@result)
  }

}

pathwayAnnotator.gmt <- function(genes,outfile,showCategory=6,pval_thres=1,qval_thres=1,gmtfile=NULL,title="",clean_desc=TRUE){
  
  print(outfile)
  
  #genes <- c('SRSF11','ACADSB','ADAMTS1','PER1','SLC16A2','SMC5','WDR6','ZNF546','ZNF729','CNTNAP1','DYX1C1','EIF4A2','EP400','IMPA2','MAGEA11','PIWIL1','BDP1','C3orf77','CDH1','CPS1','MXRA5','MYT1','USP9X','ANK2','CMYA5','DLC1','DMD','SAGE1','PIK3CA')
  
  genemap <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  geneID      <- genemap$ENTREZID
  
  if( is.null(gmtfile) ){
    #gmtfile <- "/Users/zhuy1/my_projects/LAR_project/MSigDB/c2.all.v7.1.entrez.gmt"
    gmtfile <- "/home/zhuy1/my_databases/Oncogenic_pathway/cell_2018_173_Sanchez-Vega.entrezid.gmt"
  }
  gmt <- read.gmt(gmtfile)
  
  enrichRes <- enricher(gene=geneID,
                   pvalueCutoff=pval_thres,
                   qvalueCutoff = qval_thres,
                   minGSSize = 1,
                   pAdjustMethod = "BH",
                   TERM2GENE=gmt)
  
  if(!is.null(enrichRes)){
    enrichRes <- setReadable(enrichRes,org.Hs.eg.db,keyType = "ENTREZID")
    
    if( isTRUE(clean_desc) ){
      enrichRes@result$Description <- gsub("_","",
                                      gsub("^GO_","",
                                      gsub("^KEGG_","",
                                      gsub("^REACTOME_","",
                                      enrichRes@result$Description))))
    }
    enrichRes@result$Description <- gsub("_"," ",enrichRes@result$Description)

    enrichRes@result$Description <- wrap.labels(enrichRes@result$Description,35)
    enrichRes <- filter_ORA(enrichRes)
    write.csv(enrichRes,paste0(outfile,".csv"),row.names=F)
    
    if( nrow(enrichRes@result)>0 ){
      
      # get first hit of genes
      hits <- sapply(unique(enrichRes@result$geneID),function(x) grep(x,enrichRes@result$geneID)[1])
      
      enrichRes@result <- unique(enrichRes@result[hits,])
      
      grDevices::pdf(paste0(outfile,".1.pdf"),7,5)
      print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=10))
      print(cnetplot(enrichRes,showCategory = showCategory))
      print(emapplot(enrichRes,showCategory = showCategory))
      print(heatplot(enrichRes,showCategory = showCategory))
      print(upsetplot(enrichRes,n=showCategory))
      grDevices::dev.off()
      
      grDevices::pdf(paste(outfile,"2.pdf",sep="."),6,4)
      print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=10))
      print(cnetplot(enrichRes,showCategory = showCategory))
      print(emapplot(enrichRes,showCategory = showCategory))
      print(heatplot(enrichRes,showCategory = showCategory))
      print(upsetplot(enrichRes,n=showCategory))
      grDevices::dev.off()
      
      grDevices::pdf(paste(outfile,"3.pdf",sep="."),10,10)
      print(my.barplot.gene_ratio.label_pval(enrichRes@result,title=title,n=showCategory,font.size=14))
      print(cnetplot(enrichRes,showCategory = showCategory))
      print(emapplot(enrichRes,showCategory = showCategory))
      print(heatplot(enrichRes,showCategory = showCategory))
      print(upsetplot(enrichRes,n=showCategory))
      grDevices::dev.off()
      
    }
    
    #return(enrichRes@result)
  }
  
}

get_geneList <- function(input){

  # load output of DE analysis
  data <- read.csv(input)
  
  #data <- data[order(data$pval_t.test),] 
  
  data <- data[data$gene!="HBD" & data$gene!="MEMO1",] # there are two EntrezID for HBD and MEMO1
  data <- data[data$gene!="MMD2" & data$gene!="TEC",] # there are two EntrezID for MMD2 and TEC
  
  genemap <- bitr(data$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  
  data$EntrezID <- genemap$ENTREZID
  
  data <- data[!is.na(data$EntrezID),]
  
  geneList <- data$log2fc
  names(geneList) <- data$EntrezID
  geneList <- sort(geneList,decreasing = TRUE)
  
  return(geneList)

}

get_geneList_limma <- function(input){
  
  # load output of DE analysis
  data <- read.csv(input)
  
  data <- data[order(data$adj.P.Val),] 
  
  data <- data[data$X!="HBD" & data$X!="MEMO1",] # there are two EntrezID for HBD and MEMO1
  data <- data[data$X!="MMD2" & data$X!="TEC",] # there are two EntrezID for MMD2 and TEC
  
  genemap <- bitr(data$X,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  
  data$EntrezID <- genemap$ENTREZID
  
  data <- data[!is.na(data$EntrezID),]
  
  write.table(data$EntrezID,paste(input,"ENTREZID.txt",sep="."),quote=FALSE,row.names = F,col.names = F)
  
  geneList <- data$logFC
  names(geneList) <- data$EntrezID
  geneList <- sort(geneList,decreasing = TRUE)
  
  return(geneList)
  
}

pathwayAnnotator.gsea.gmt <- function(limma_output,outprefix,showCategory=6,pval_thres=1,qval_thres=1,gmtfile=NULL,selected_pathways=NULL){
  
  geneList <- get_geneList_limma(limma_output)
  
  if( is.null(gmtfile) ){
    gmtfile <- "~/my_databases/Oncogenic_pathway/cell_2018_173_Sanchez-Vega.entrezid.gmt"
    gmtfile <- "~/my_databases/MSigDB/c2.all.v7.1.entrez.gmt"
  }
  gmt <- read.gmt(gmtfile)
  
  if( !is.null(selected_pathways) ){
    gmt = gmt[gmt$ont %in% selected_pathways,]
  }
  #print(head(gmt))
  write.csv(gmt,paste(outprefix,"gmt",sep="."),row.names=F)
  
  ggmt <- GSEA(geneList, TERM2GENE=gmt, nPerm = 10000, minGSSize = 1, maxGSSize = 1000,pvalueCutoff = 1, verbose=FALSE)
  ggmt <- setReadable(ggmt,org.Hs.eg.db,keyType = "ENTREZID")
  
  write.csv(ggmt,paste(outprefix,"GSEA.csv",sep="."),row.names=F)

  
  if( nrow(ggmt@result)>0 ){
    pdf(paste(outprefix,"GSEA.pdf",sep="."))
    print(my.barplot.gsea.enrichment_score(ggmt,title="GSEA",n=showCategory))
    library(scales)
    print(cnetplot(ggmt, categorySize="enrichmentScore", foldChange=squish(geneList,c(-2,2)),circular=FALSE,showCategory = showCategory))
    print(cnetplot(ggmt, categorySize="enrichmentScore", foldChange=geneList,circular=FALSE,showCategory = showCategory))
    print(emapplot(ggmt,showCategory = showCategory))
    print(heatplot(ggmt,foldChange = geneList,showCategory = showCategory))

    for( i in seq(1,nrow(ggmt@result),1) ){
      j = i + 0
      if( j > nrow(ggmt@result) )
        j = nrow(ggmt@result)
      print(c(i,j))
      ggmt@result$Description <- wrap.labels(ggmt@result$Description,20)
      print(gseaplot2(ggmt,geneSetID = i:j,title=ggmt@result[i,1],pvalue_table = TRUE,subplots=c(1,2,3),base_size = 15,color="blue"))
    }
    
    dev.off()
  }
}

do_something <- function(){
  
  pathwayAnnotator.gsea.gmt("AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_TNBC/DEanalysis_LAR_TNBC_vs_Non-LAR_TNBC.fit.csv","AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_TNBC/DEanalysis_LAR_TNBC_vs_Non-LAR_TNBC.fit.oncopathways",showCategory = 10)
  pathwayAnnotator.gsea.gmt("AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumA/DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumA.fit.csv","AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumA/DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumA.fit.oncopathways",showCategory = 10)
  pathwayAnnotator.gsea.gmt("AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumB/DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumB.fit.csv","AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumB/DEanalysis_LAR_TNBC_vs_Non-LAR_ERpos_LumB.fit.oncopathways",showCategory = 10)
  pathwayAnnotator.gsea.gmt("AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_Non-TNBC_HER2/DEanalysis_LAR_TNBC_vs_Non-LAR_Non-TNBC_HER2.fit.csv","AR_pathway.AR.DEanalysis_LAR_TNBC_vs_Non-LAR_Non-TNBC_HER2/DEanalysis_LAR_TNBC_vs_Non-LAR_Non-TNBC_HER2.fit.oncopathways",showCategory = 10)
  
  
}


pathwayAnnotator.ORA <- function(genes,outfile,showCategory=6,pval_thres=1,qval_thres=1,gmtfile=NULL,title="",clean_desc=TRUE){
  
  print(outfile)
  
  #genes <- c('SRSF11','ACADSB','ADAMTS1','PER1','SLC16A2','SMC5','WDR6','ZNF546','ZNF729','CNTNAP1','DYX1C1','EIF4A2','EP400','IMPA2','MAGEA11','PIWIL1','BDP1','C3orf77','CDH1','CPS1','MXRA5','MYT1','USP9X','ANK2','CMYA5','DLC1','DMD','SAGE1','PIK3CA')
  
  genemap <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
  geneID      <- genemap$ENTREZID
  
  GO_over_representative_test(geneID,geneID,outfile,use_geneList = FALSE)
  
  KEGG_over_representative_test(geneID,geneID,outfile,use_geneList = FALSE)
  
  Reactome_enrichPathway(geneID,geneID,outfile,use_geneList = FALSE)
  
  # oncogenic pathway
  #gmtfile <- "/home/zhuy1/my_databases/Oncogenic_pathway/cell_2018_173_Sanchez-Vega.entrezid.gmt"
  #pathwayAnnotator.gmt(genes,paste0(outfile,".oncopath"),showCategory=10,pval_thres=1,qval_thres=1,gmtfile=NULL,title="Oncogenic pathway",clean_desc=FALSE)
  
  # AR gene sets
  #gmtfile <- "~/my_projects/LAR_project/MSigDB_Androgen/merged_AR_genesets.entrez.gmt"
  #pathwayAnnotator.gmt(genes,paste0(outfile,".AR_genesets"),showCategory=10,pval_thres=1,qval_thres=1,gmtfile=gmtfile,title="AR gene sets",clean_desc=FALSE)
  
  # AR responsive genes
  #gmtfile <- "~/my_projects/LAR_project/MSigDB_Androgen/DOANE_NELSON.entrez.gmt"
  #pathwayAnnotator.gmt(genes,paste0(outfile,".AR_response"),showCategory=10,pval_thres=1,qval_thres=1,gmtfile=gmtfile,title="AR responsive genes",clean_desc=FALSE)
}
