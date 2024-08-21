# returns diffTable containing top two most up-/down-regulated genes based on log2FC
getDiffTop <- function (complete, alpha = 0.05) 
{
  diff.table <- list();
  for (name in names(complete)) {
    complete.name <- complete[[name]]
    sample.ref <- gsub("[^0-9a-zA-Z]","",strsplit(name,'vs')[[1]][2])
    sample.treat <- gsub("[^0-9a-zA-Z]","",strsplit(name,'vs')[[1]][1])
    idx.sample.ref = grep(sample.ref,gsub("[^0-9a-zA-Z]","",names(complete.name)))
    idx.sample.ref.norm = idx.sample.ref[grep('norm',names(complete.name[idx.sample.ref]))]
    idx.sample.ref = idx.sample.ref[!grepl('norm',names(complete.name[idx.sample.ref]))]
    idx.sample.treat = grep(sample.treat,gsub("[^0-9a-zA-Z]","",names(complete.name)))
    idx.sample.treat.norm = idx.sample.treat[grepl('norm',names(complete.name[idx.sample.treat]))]
    idx.sample.treat = idx.sample.treat[!grepl('norm',names(complete.name[idx.sample.treat]))]
    idx.complete = c(1,2,idx.sample.ref,idx.sample.treat,idx.sample.ref.norm,
                     idx.sample.treat.norm,
                     (ncol(complete.name)-9):ncol(complete.name))
    idx.complete=sort(idx.complete)
    # renaming columns to normalized mean of each condition
    colnames(complete.name)[gsub("[^0-9a-zA-Z]","",colnames(complete.name))==sample.ref]=
      paste0("mean.norm.",gsub(" ",".",
                               colnames(complete.name)[gsub("[^0-9a-zA-Z]","",colnames(complete.name))==sample.ref]))
    colnames(complete.name)[gsub("[^0-9a-zA-Z]","",colnames(complete.name))==sample.treat]=
      paste0("mean.norm.",
             gsub(" ",".",colnames(complete.name)[gsub("[^0-9a-zA-Z]","",colnames(complete.name))==sample.treat]))
    if (!is.null(complete.name$betaConv)) {
      up = complete.name[which(complete.name$padj <= 
                                 alpha & complete.name$betaConv & complete.name$log2FoldChange>0),]
      diff.table[[name]]=rbind(diff.table[[name]],head(up[order(up$log2FoldChange,decreasing=TRUE),idx.complete],n=2))
      diff.table[[name]]=rbind(diff.table[[name]],tail(up[order(up$log2FoldChange,decreasing=TRUE),idx.complete],n=2))
      dwn = complete.name[which(complete.name$padj <= 
                                 alpha & complete.name$betaConv & complete.name$log2FoldChange<0),]
      diff.table[[name]]=rbind(diff.table[[name]],head(dwn[order(dwn$log2FoldChange,decreasing=FALSE),idx.complete],n=2))
      diff.table[[name]]=rbind(diff.table[[name]],tail(dwn[order(dwn$log2FoldChange,decreasing=FALSE),idx.complete],n=2))
    }
    else {
      up = complete.name[which(complete.name$padj <= 
                                 alpha & complete.name$log2FoldChange>0),]
      diff.table[[name]]=rbind(diff.table[[name]],head(up[order(up$log2FoldChange,decreasing=TRUE),idx.complete],n=2))
      diff.table[[name]]=rbind(diff.table[[name]],tail(up[order(up$log2FoldChange,decreasing=TRUE),idx.complete],n=2))
      dwn = complete.name[which(complete.name$padj <= 
                                  alpha &  complete.name$log2FoldChange<0),]
      diff.table[[name]]=rbind(diff.table[[name]],head(dwn[order(dwn$log2FoldChange,decreasing=FALSE),idx.complete],n=2))
      diff.table[[name]]=rbind(diff.table[[name]],tail(dwn[order(dwn$log2FoldChange,decreasing=FALSE),idx.complete],n=2))
    }
  }
  
  return(diff.table)
}