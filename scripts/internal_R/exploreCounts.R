exploreCounts <- function (object, group, typeTrans = "VST", gene.selection = "pairwise", 
          col , batch,varInt,batchRem=FALSE) 
{
 if (class(object) == "DESeqDataSet") {
  if (typeTrans == "VST") 
   counts.trans <- assay(varianceStabilizingTransformation(object))
  else counts.trans <- assay(rlogTransformation(object))
  if(batchRem){
   PCAPlot(object, group=group,counts.trans, col=col,varInt,typeTrans,batch=batch,batchRem=TRUE)
   clusterPlot(counts.trans = counts.trans, group = group)
  }else{
   PCAPlot(object, group=group,counts.trans,col=col,varInt,typeTrans,batch=batch,batchRem=FALSE)
   clusterPlot(counts.trans = counts.trans, group = group)
  }
 }
 else if (class(object) == "DGEList") {
  MDSPlot(dge = object, group = group, col = col, gene.selection = gene.selection)
  clusterPlot(counts.trans = cpm(object, prior.count = 2, 
   log = TRUE), group = group)
 }
 else {
  stop("The object is not a DESeqDataSet nor a DGEList")
 }
}
