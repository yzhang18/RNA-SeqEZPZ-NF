run.DESeq2 <- function (counts, target, varInt, batch = NULL, locfunc = "median", 
          fitType = "parametric", pAdjustMethod = "BH", cooksCutoff = TRUE, 
          independentFiltering = TRUE, alpha = 0.05, ...) 
{
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = target, 
                                design = formula(paste("~", ifelse(!is.null(batch), paste(batch, 
                                                                                          "+"), ""), varInt)))
  cat("Number of rows before pre-filter",nrow(dds),"\n")
  dds <- dds[rowSums(counts(dds))>1,]
  cat("Number of rows after pre-filter", nrow(dds),"\n")
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)), collapse = " "), "\n")
  dds <- estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  dds <- estimateDispersions(dds, fitType = fitType)
  dds <- nbinomWaldTest(dds, ...)
  results <- list()
  for (comp in combn(nlevels(colData(dds)[, varInt]), 2, simplify = FALSE)) {
    levelRef <- levels(colData(dds)[, varInt])[comp[1]]
    levelTest <- levels(colData(dds)[, varInt])[comp[2]]
    results[[paste0(levelTest, "_vs_", levelRef)]] <- results(dds, 
	contrast = c(varInt, levelTest, levelRef), pAdjustMethod = pAdjustMethod, 
	cooksCutoff = cooksCutoff, independentFiltering = independentFiltering, alpha = alpha)
    cat(paste("Comparison", levelTest, "vs", levelRef, "done\n"))
  }
  return(list(dds = dds, results = results, sf = sizeFactors(dds)))
}
