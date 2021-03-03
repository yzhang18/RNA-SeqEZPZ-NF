volcanoPlot <- function (complete, alpha = 0.05, outfile = TRUE,fc.cutoff=1) 
{
  ncol <- ifelse(length(complete) <= 4, ceiling(sqrt(length(complete))), 
                 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) 
    png(filename = "figures/volcanoPlot.png", width = cairoSizeWrapper(1800 * 
      ncol), height = cairoSizeWrapper(1800 * nrow), res = 300)
  par(mfrow = c(nrow, ncol))
  for (name in names(complete)) {
    complete.name <- complete[[name]]
    complete.name$padj[which(complete.name$padj == 0)] <- .Machine$double.xmin
    log10pval <- -log10(complete.name$padj)
    ylim <- c(0, 1) * quantile(log10pval, probs = 0.99, na.rm = TRUE)
    plot(complete.name$log2FoldChange, pmin(ylim[2], log10pval), 
         ylim = ylim, las = 1, cex = 0.45, xlab = expression(log[2] ~ 
         fold ~ change), ylab = expression(-log[10] ~ 
         adjusted ~ P ~ value), col = ifelse(complete.name$padj <= alpha &
                                               abs(complete.name$log2FoldChange) > (log2(fc.cutoff)), "red", "black"), 
         pch = ifelse(log10pval >=ylim[2], 2, 20), main = paste0("Volcano plot - ",                                                                                                                                                                                                                gsub("_", " ", name)))
    abline(h = -log10(alpha), lty = 2, col = "lightgray")
  }
  if (outfile) 
    dev.off()
}