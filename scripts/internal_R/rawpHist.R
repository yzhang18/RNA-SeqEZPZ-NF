rawpHist <- function (complete, outfile = TRUE,fdrtool.group=NULL,out.DESeq2=NULL) 
{

  ncol <- ifelse(length(complete) <= 4, ceiling(sqrt(length(complete))), 
                 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile & is.null(fdrtool.group)) 
    png(filename = "figures/rawpHist.png", width = cairoSizeWrapper(1800 * 
                                                                      ncol), height = cairoSizeWrapper(1800 * nrow), res = 300)
  if (outfile & !is.null(fdrtool.group)) 
    png(filename = "figures/rawpHist_corrected.png", width = cairoSizeWrapper(1800 * 
                                                                      ncol), height = cairoSizeWrapper(1800 * nrow), res = 300)
  
  par(mfrow = c(nrow, ncol))
  for (name in names(complete)) {
    if(!(name %in% fdrtool.group)){
    hist(complete[[name]][complete[[name]]$baseMean>5, "pvalue"], breaks=0:20/20, xlab = "Raw p-value", 
         col = "skyblue", las = 1, main = paste0("Distribution of raw p-values - ", 
                                                 gsub("_", " ", name)))
    }else if(name %in% fdrtool.group) {
      # remove genes filtered out by independent filtering and the dispersion outliers
      # save index where padj will be adjusted
      idx.adj <- which(!is.na(out.DESeq2$results[[name]]$padj) & !is.na(out.DESeq2$results[[name]]$pvalue))
      DESeq2Res <- out.DESeq2$results[[name]][idx.adj,]
      # remove original adjusted p-values
      DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
      FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = F)
      # estimated null model variance. Theoretical is 1
      # FDR.DESeq2Res$param[1, "sd"]
      # add the new BH-adjusted p-values and add results to results data frame
      DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
      out.DESeq2$results[[name]][idx.adj,"padj"] = DESeq2Res[,"padj"]
      out.DESeq2$results[[name]][idx.adj,"pval"] = FDR.DESeq2Res$pval
      # plot histogram of "correct" p-values
      hist(out.DESeq2$results[[name]][out.DESeq2$results[[name]]$baseMean>5,"pval"],breaks=0:20/20,
           col = "royalblue4", 
           main = paste0("Distribution of corrected null model - ", 
                         gsub("_", " ", name)), xlab = "CORRECTED p-values")
    }
  }
  if (outfile) 
    dev.off()
  return(out.DESeq2)
}