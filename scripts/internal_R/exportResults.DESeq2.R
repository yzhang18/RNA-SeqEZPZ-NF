### modified CT 2018/06/20 adding gene symbol changing fold-change
#### modified CT 2018/08/13 changing fold-change when it is less than 0.
### modified CT 2021/02/02 fixed fold-change when it gets rounded to zero
exportResults.DESeq2 <- function (out.DESeq2, group, alpha = 0.05, export = TRUE,fc.cutoff=1) 
{
 
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  counts <- data.frame(Id = rownames(counts(dds)), counts(dds), 
                       round(counts(dds, normalized = TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", 
                                                            colnames(counts(dds))))
  bm <- data.frame(Id = rownames(results[[1]]), baseMean = round(results[[1]][, 
                                                                              "baseMean"], 2))
  base <- merge(counts, bm, by = "Id", all = TRUE)
  tmp <- base[, paste("norm", colnames(counts(dds)), sep = ".")]
  for (cond in levels(group)) {
    base[, cond] <- round(apply(as.data.frame(tmp[, group == 
                                                    cond]), 1, mean), 0)
  }
  ### modified CT 2018/06/20 adding gene symbol changing fold-change
  ### 2021/02/02 fixed fold-change when it's close to zero
  #library(biomaRt)
  library(tibble)
  library(EnsDb.Hsapiens.v86)
  
  #ensembl.hs <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ###
  complete <- list()
  for (name in names(results)) {
    complete.name <- base
    
    res.name <- data.frame(Id = rownames(results[[name]]), 
                           FoldChange = 2^(results[[name]][, "log2FoldChange"]),
                           log2FoldChange = round(results[[name]][,"log2FoldChange"], 3), 
                           pvalue = results[[name]][,"pvalue"], 
                           padj = results[[name]][, "padj"])
    #### modified CT 2018/08/13 changing fold-change when it is less than 1.
    res.name$FoldChange=ifelse(res.name$FoldChange<1,-1/res.name$FoldChange,res.name$FoldChange)
    # round fold-change after converting the sign
    res.name$FoldChange=round(res.name$FoldChange,3)
    complete.name <- merge(complete.name, res.name, by = "Id", 
                           all = TRUE)
    mcols.add <- data.frame(Id = rownames(counts(dds)), 
                            dispGeneEst = round(mcols(dds)$dispGeneEst, 4), 
                            dispFit = round(mcols(dds)$dispFit, 4), 
                            dispMAP = round(mcols(dds)$dispMAP, 4), 
                            dispersion = round(mcols(dds)$dispersion, 4), 
                            betaConv = mcols(dds)$betaConv, 
                            maxCooks = round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by = "Id", 
                           all = TRUE)
    ##-- modified by CT 2019/05/15
    # browser()
    gene.id = complete.name$Id
    gene.id = as.character(gene.id)
    # removed "." in gene.id (CT 2022/02/21)
    gene.id = sapply(strsplit(gene.id,"[.]"), function(x) x[1])
    ## query by the ensembl gene id
    ## filters letting you know what type of input
    ## attributes what should be displayed (listAttributes(ensembl.hs) to get other attributes)
    ## values is the input vector
    #tmp.gene <- getBM(filters= "ensembl_gene_id",
    #                   attributes= c("ensembl_gene_id", "chromosome_name","start_position","end_position","hgnc_symbol"),
    #                    values=gene.id,mart= ensembl.hs)
    #geneSymbol = tmp.gene$hgnc_symbol[match(complete.name$Id,tmp.gene$ensembl_gene_id)]
    
    hsens=EnsDb.Hsapiens.v86
    
    tmp.gene=select(hsens,  
                    keys = gene.id, 
                    columns = c("ENTREZID", "GENENAME", "GENEID"), 
                    keytype = "GENEID")
    complete.name.id.nodot=sapply(strsplit(complete.name$Id,"[.]"),function(x) x[1])
    geneSymbol=tmp.gene$GENENAME[match(complete.name.id.nodot,tmp.gene$GENEID)]
    
    complete.name=add_column(complete.name,geneSymbol,.after=1)
    #--
    complete[[name]] <- complete.name
    if (export) {
      up.name <- complete.name[which(complete.name$padj <= 
                                       alpha & complete.name$betaConv & 
                                      complete.name$log2FoldChange >
                                       (log2(fc.cutoff))), ]
      up.name <- up.name[order(up.name$padj), ]
      down.name <- complete.name[which(complete.name$padj <= 
                                         alpha & complete.name$betaConv & 
                                        complete.name$log2FoldChange < 
                                         (-log2(fc.cutoff))), ]
      down.name <- down.name[order(down.name$padj), ]
      name <- gsub("_", "", name)
      write.table(complete.name, file = paste0("tables/",name, ".complete.txt"), 
                  sep = "\t", row.names = FALSE, 
                  dec = ".", quote = FALSE)
      write.table(up.name, file = paste0("tables/",name, 
                                         ".up.txt"), row.names = FALSE, sep = "\t", dec = ".", 
                  quote = FALSE)
      write.table(down.name, file = paste0("tables/",name, 
                                           ".down.txt"), row.names = FALSE, sep = "\t", 
                  dec = ".", quote = FALSE)
    }
  }
  return(complete)
}
