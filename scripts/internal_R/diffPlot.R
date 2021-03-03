# for diffPlot to work, group names must be part of sample names
# i.e. if group names 'hairpin3', sample names must contains hairpin3

diffPlot <- function (out.DESeq2, alpha = 0.05, outfile = TRUE,fc.cutoff=1) 
{
  dds <- out.DESeq2$dds
  nrow <- 2 # for up and down genes
  if (outfile) 
    pdf(file= "figures/DiffPlot.pdf", height=11,width=8)
  # get the log transformed counts
  rld <- rlog(dds, blind=F)
  plot.lst=list()
  for (comp in combn(nlevels(colData(dds)[, varInt]), 2, simplify = FALSE)) {
    levelRef <- levels(colData(dds)[, varInt])[comp[1]]
    levelTest <- levels(colData(dds)[, varInt])[comp[2]]
    results <- subset(results(dds, contrast=c("group",levelTest,levelRef)), padj < alpha)
    # get appropriate columns
    # correspond to the samples (to pull from rld)
    aCols <- grep(gsub("[^0-9a-zA-Z]","",levelTest),gsub("[^0-9a-zA-Z]","",colData(dds)$group))
    bCols <- grep(gsub("[^0-9a-zA-Z]","",levelRef),gsub("[^0-9a-zA-Z]","",colData(dds)$group))
    colNums=c(aCols,bCols); 
    # make the lists of up genes
    tmp.upgenes <- results[results$log2FoldChange>(log2(fc.cutoff)),]
    upgenes <- rownames(tmp.upgenes[order(tmp.upgenes$log2FoldChange,decreasing=TRUE),])
    
    # this gives us the rows we want
    rows <- match(upgenes, row.names(rld))
    mat <- assay(rld)[rows,colNums]
    # df <- data.frame(factor(colData(rld)["group"][colNums,]),
		# factor(colData(rld)["rep"][colNums,]))
	# names(df)=c("group","rep")
    df <- as.data.frame(colData(rld)[c("group","rep")][colNums,])

	# # number of replicates
	# nbatch=dim(unique(df["rep"]))[1]
	# # number of groups
	# ngroups=dim(unique(df["group"]))[1]
	# # color for annotation columns
	# plot.colors = brewer.pal(nbatch+ngroups+2,"Paired")
	
	# annot.colors = list(
    # plot.colors[3:(2+ngroups)],
    # plot.colors[(2+ngroups+1):(2+ngroups+nbatch)]
    # )
	# names(annot.colors)[1]=names(df["group"])
	# names(annot.colors)[2]=names(df["rep"])
	# names(annot.colors[[1]])=levels(df[,1])
	# names(annot.colors[[2]])=levels(df[,1])
	if(is.null(nrow(mat))){
		plot.lst[[1]]=textGrob("No significant up-regulated genes.");
	}else if(nrow(mat)<2){
		plot.lst[[1]]=textGrob("Significant up-regulated genes < 2. No heatmap was generated");
	}else{
		legend_labels=as.character(seq(round(min(mat)),round(max(mat)),by=2))
		legend_labels[length(legend_labels)]=paste0("rlog normalized counts\n",round(max(mat)))
		plot.lst[[1]]=pheatmap(mat, fontsize=12,annotation_col=df, 
			annotation_names_col = FALSE,
			legend_breaks=seq(round(min(mat)),round(max(mat)),by=2),
			legend_labels=legend_labels,
             main=paste("Up genes. Ref=",levelRef),
             show_rownames = FALSE,show_colnames = FALSE,
             color = colorRampPalette(brewer.pal(n=9,'YlGnBu'))(100),
			 #color = colorRampPalette(c("white",plot.colors[2]))(100),
			 #annotation_colors = annot.colors,
			 silent=TRUE)[[4]]
    }
    # make the lists of down genes
    tmp.dwngenes <- results[-results$log2FoldChange>(log2(fc.cutoff)),]
    dwngenes <- rownames(tmp.dwngenes[order(-tmp.dwngenes$log2FoldChange,decreasing=TRUE),])
    
    # this gives us the rows we want
    rows <- match(dwngenes, row.names(rld))
    mat <- assay(rld)[rows,colNums]
    
    df <- as.data.frame(colData(rld)[c("group","rep")][colNums,])
    if(is.null(nrow(mat))){
		plot.lst[[2]]=textGrob("No significant down-regulated genes.");
	}else if(nrow(mat)<2){
		plot.lst[[2]]=textGrob("Significant down-regulated genes < 2. No heatmap was generated");
	}else{
		legend_labels=as.character(seq(round(min(mat)),round(max(mat)),by=2))
		legend_labels[length(legend_labels)]=paste0("rlog normalized counts\n",round(max(mat)))
		plot.lst[[2]]=pheatmap(mat,annotation_col = df, fontsize=12, 
			annotation_names_col = FALSE,
			legend_breaks=seq(round(min(mat)),round(max(mat)),by=2),
			legend_labels=legend_labels,
            main=paste("Down genes. Ref=",levelRef),
            show_rownames = FALSE,show_colnames = FALSE,
            color = colorRampPalette(brewer.pal(n=9,'YlGnBu'))(100),
			#color = colorRampPalette(rev(plot.colors[1:2]))(100),
			#annotation_colors=annot.colors,
			silent=TRUE)[[4]]
	 }
  grid.arrange(plot.lst[[1]],plot.lst[[2]],nrow=2)
    }
  if (outfile){
	tmp=dev.off()
	while(tmp>1) tmp=dev.off()
  }
  
}
