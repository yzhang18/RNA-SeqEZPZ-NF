summarizeResults.DESeq2 <- function (out.DESeq2, group, independentFiltering = TRUE, cooksCutoff = TRUE, 
                                     alpha = 0.05, col = c("lightblue", "orange", "MediumVioletRed", 
                                                           "SpringGreen"),fdrtool.group=NULL) 
{
  if (!I("figures" %in% dir())) 
    dir.create("figures", showWarnings = FALSE)
  if (!I("tables" %in% dir())) 
    dir.create("tables", showWarnings = FALSE)
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  diagSizeFactorsPlots(dds = dds, group = group, col = col)
  countsBoxplots(dds, group = group, col = col)
  dispersionsPlot(dds = dds)
  if (independentFiltering) {
    tabIndepFiltering <- tabIndepFiltering(results)
    cat("Number of features discarded by the independent filtering:\n")
    print(tabIndepFiltering, quote = FALSE)
  }
  else {
    tabIndepFiltering <- NULL
  }
  complete <- exportResults.DESeq2(out.DESeq2, group = group, 
                                   alpha = alpha)
  nDiffTotal <- calcDiffTotal(complete = complete, alpha = alpha)
  cat("\nNumber of features down/up and total:\n")
  print(nDiffTotal, quote = FALSE)
  # modified by CT 5/14/2019
  diff.table <- getDiffTop(complete=complete, alpha = alpha) 
  # if fdrtool.group not NULL re-estimate adjusted pval
  if(is.null(fdrtool.group)){
    rawpHist(complete = complete)
  }else{
    out.DESeq2=rawpHist(complete = complete,fdrtool.group=fdrtool.group,out.DESeq2=out.DESeq2)
    complete <- exportResults.DESeq2(out.DESeq2, group = group, 
                                     alpha = alpha)
    nDiffTotal <- calcDiffTotal(complete = complete, alpha = alpha)
    cat("\nNumber of features down/up and total:\n")
    print(nDiffTotal, quote = FALSE)
    diff.table <- getDiffTop(complete=complete, alpha = alpha)
  }
    diffPlot(out.DESeq2 = out.DESeq2,alpha=alpha)
  MAPlot(complete = complete, alpha = alpha)
  volcanoPlot(complete = complete, alpha = alpha)
  return(list(complete = complete, tabIndepFiltering = tabIndepFiltering, 
              nDiffTotal = nDiffTotal,diffTable=diff.table,out.DESeq2=out.DESeq2))
}