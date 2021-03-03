descriptionPlots <- function (counts, group, col , ggplot_theme = theme_gray())
{
    if (!I("figures" %in% dir()))
        dir.create("figures", showWarnings = FALSE)
    barplotTotal(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
    barplotNull(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
    densityPlot(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
    #majSequences <- majSequences(counts = counts, group = group,
    #    col = col, ggplot_theme = ggplot_theme)
    majSequences <- majSequences(counts=counts,col=col,group=group,outfile=TRUE)
	cat("Matrix of SERE statistics:\n")
    print(tabSERE(counts))
    #pairwiseScatterPlots(counts = counts)
	pairwiseScatterPlots(counts=counts, group=group, outfile = TRUE)
    return(majSequences)
}