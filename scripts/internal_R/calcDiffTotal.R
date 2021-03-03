calcDiffTotal <- function (complete, alpha = 0.05) 
{
  nDiffTotal<- matrix(NA, ncol = 5, nrow = length(complete), 
                       dimnames = list(names(complete), c("Test vs Ref","min baseMean", "# down", 
                                                          "# up", "# total")))
  for (name in names(complete)) {
    complete.name <- complete[[name]]
    if (!is.null(complete.name$betaConv)) {
      nDiffTotal[name, 2:4] = c(
        min(complete.name[which(complete.name$padj <= 
                                  alpha & complete.name$betaConv), "baseMean"]),
        nrow(complete.name[which(complete.name$padj <= 
           alpha & complete.name$betaConv & complete.name$log2FoldChange <= 0), ]), 
        nrow(complete.name[which(complete.name$padj <= 
               alpha & complete.name$betaConv & complete.name$log2FoldChange >= 0), ]))
    }
    else {
      nDiffTotal[name, 2:4] = c(
        min(complete.name[which(complete.name$padj <= 
             alpha),"baseMean" ]),
        nrow(complete.name[which(complete.name$padj <= 
             alpha & complete.name$log2FoldChange <= 0), ]), 
        nrow(complete.name[which(complete.name$padj <= 
             alpha & complete.name$log2FoldChange >= 0), ]))
    }
  }
  nDiffTotal[, 5] <- nDiffTotal[, 3] + nDiffTotal[, 4]
  nDiffTotal[, 1] <- gsub("_", " ", rownames(nDiffTotal))
  rownames(nDiffTotal) <- NULL
  return(nDiffTotal)
}