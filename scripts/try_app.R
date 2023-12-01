if(length(.libPaths())>1) .libPaths(.libPaths()[-1])
library(shiny)
library(GeneOverlap)
library(gridExtra)
library(ggplot2)
library(gridGraphics)
library(eulerr)
library(reshape)
library(ggthemes)
library(ggpointdensity)
library(shiny)
library(RColorBrewer)
library(grid)
library(gridBase)
library(venn)
library(UpSetR)
library(DESeq2)
library(clusterProfiler)
library(msigdbr)
library(stringr)

# filename was changed. Make sure can read both
if(file.exists("/mnt/outputs/diff_analysis_rslt/RNA-seq_differential_analysis.RData")){
 load("/mnt/outputs/diff_analysis_rslt/RNA-seq_differential_analysis.RData")
} else {
 load("/mnt/outputs/diff_analysis_rslt/RNA-seq differential analysis.RData")
}

grp.name=names(out.DESeq2$results)
idx=match(grp.name,names(out.DESeq2$results))
results=out.DESeq2$results[idx]
fdr.co=rep(1,6)
fc.cutoff=rep(1,6)
meanDiff.cutoff=rep(0,6)
normCts <- counts(out.DESeq2$dds,normalized=TRUE)

grp.up=list()
for(i in 1:length(grp.name)){
# filter out by mean difference
grps=strsplit(grp.name[i],"_vs_")
treat.grp=grps[[1]][1]
ref.grp=grps[[1]][2]
treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
diff.mean=treat.mean-ref.mean
grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
!is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i],]
}
grp.dwn=list()
for(i in 1:length(grp.name)){
# filter out by mean difference
grps=strsplit(grp.name[i],"_vs_")
treat.grp=grps[[1]][1]
ref.grp=grps[[1]][2]
treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
diff.mean=ref.mean-treat.mean
grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) &
!is.na(results[[grp.name[i]]]$padj) & diff.mean >= meanDiff.cutoff[i],]
}
compare.df = NULL
grp=grp.up
for(i in 1:length(grp.name)){
compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
group1=rep(grp.name[i],dim(grp[[i]])[1]),
group2=rep("Up-regulated",dim(grp[[i]])[1])))
}

grp=grp.dwn
for(i in 1:length(grp.name)){
compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
group1=rep(grp.name[i],dim(grp[[i]])[1]),
group2=rep("Down-regulated",dim(grp[[i]])[1])))
}
enrich.pval.co=0.05
# Using clusterProfiler to perform hypergeometric test on msigdb signatures
msigdb.species <- "Danio rerio"
msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "MF") %>%
 dplyr::select(gs_name, gene_symbol)
msig.name ="MSigDB GO Molecular Function"
msig.gene.set$gs_name = gsub("GOMF_","",msig.gene.set$gs_name)
msig.gene.set$gs_name = gsub("_"," ",msig.gene.set$gs_name)
msig.gene.set$gs_name = tolower(msig.gene.set$gs_name)
msig.gene.set$gs_name = str_to_title(msig.gene.set$gs_name)
msig.gene.set$gs_name = str_replace_all(msig.gene.set$gs_name, 
                        c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", "Trna" = "tRNA",
                        "Mirna" = "miRNA", "Rrna" = "rRNA","Atp" = "ATP","Adp" = "ADP",
                        "Lncrna" = "lncRNA","Snrna" = "snRNA","Snorna" = "snoRNA"))
formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                              TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                              pAdjustMethod="BH")
dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name) + facet_grid(~group2) +
 scale_y_discrete(labels=function(x) str_wrap(x, width=40)) +
 scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
 scale_color_distiller(palette = 'Blues')
