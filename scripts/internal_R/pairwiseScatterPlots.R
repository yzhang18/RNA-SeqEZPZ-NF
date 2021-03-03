# modified 2021/02/02 by CT changing size of text inside cells
pairwiseScatterPlots <- function (counts, group, outfile = TRUE) 
{
  ncol <- ncol(counts)
  if (ncol <= 30) {
    if (outfile) 
      png(filename = "figures/pairwiseScatter.png", width = cairoSizeWrapper(700 * 
                                                                               ncol), height = cairoSizeWrapper(700 * ncol), 
          res = 300)
    panel <- function(x, y, ...) {
      points(x, y, pch = ".")
      abline(a = 0, b = 1, lty = 2)
    }
    lower.panel <- function(x, y, ...) {
      horizontal <- (par("usr")[1] + par("usr")[2])/2
      vertical <- (par("usr")[3] + par("usr")[4])/2
      text(horizontal, vertical, round(SERE(2^cbind(x, 
                                                    y) - 1), digits = 2), cex = ncol/3)
    }
    pairs(log2(counts + 1), panel = panel, lower.panel = lower.panel, 
          las = 1, labels = paste(group,rep, 
                                  sep = "\n"), main = "Pairwise scatter plot", 
          cex.labels = ncol/4, cex.main = ncol/4)
    if (outfile) 
      dev.off()
  }
  else {
    warning("No pairwise scatter-plot produced because of a too high number of samples (>30).")
    if (outfile) 
      png(filename = "figures/pairwiseScatter.png", width = 1900, 
          height = 1900, res = 300)
    par(mar = c(1.5, 1.5, 1.5, 1.5))
    plot(0, 0, bty = "o", pch = ".", col = "white", xaxt = "n", 
         yaxt = "n", xlab = "", ylab = "")
    text(0, 0.3, "No pairwise scatter-plot produced because of\na too high number of samples (>30).", 
         pos = 1, cex = 1.5)
    if (outfile) 
      dev.off()
  }
}