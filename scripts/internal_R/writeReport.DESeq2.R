writeReport.DESeq2 <-
function (target, counts, out.DESeq2, summaryResults, majSequences,
    workDir, projectName, author, targetFile, rawDir, featuresToRemove,
    varInt, condRef, batch, fitType, cooksCutoff, independentFiltering,
    alpha, pAdjustMethod, typeTrans, locfunc, colors)
{
    rmarkdown::render(input = "/scripts/internal_R/report_DESeq2.rmd",
        output_file = paste0(projectName,
        "_report.html"), output_dir = workDir, intermediates_dir = workDir,
        knit_root_dir = workDir, run_pandoc = TRUE, quiet = TRUE,
        clean = TRUE)
    cat("HTML report created\n")
}