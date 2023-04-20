
#############function
DEGanalysis_t <- function(exprMat, group, thr.filter=log2(0.001), thr.FC=2, thr.p=0.05){
  
  ### This function is for identify DEG from a given matrix
  ### Parameters:
  ###   Exprmat:  should a matrix in stead of a data frame, 
  ###             with column for samples and rows for signature(gene/protein/metabolics...)
  ###             For transcriptome, a matrix after log transformation is recommended.
  ###   Group:    should be a factor whose length is identical to the number of the columns in exprMat,
  ###             describing the group information of each column in exprMat
  ###   thr.filter: genes with average lower than thr.filter in both groups will not be considered as DEGs
  ###
  ### Examples:
  ###   findDEG(exprMat = log2(exprMat_RNA_fpkm_D5D6+0.01), 
  ###           group = factor(dt.meta[colnames(exprMat_RNA_fpkm_D5D6)]$group) )
  ###
  
  
  if(!require("data.table")) install.packages("data.table")
  library(data.table)
  
  group_A = levels(group)[1]
  group_B = levels(group)[2]
  col_A = which(as.numeric(group)==1)
  col_B = which(as.numeric(group)==2)
  exprMat_filt <- exprMat[apply(exprMat,1,sd)>0,]
  
  
  dt.DEG <- apply(exprMat_filt, 1, function(x){
    
    mean = mean(x[c(col_A,col_B)])
    sd = sd(x[c(col_A,col_B)])
    mean_A = mean(x[col_A])
    mean_B = mean(x[col_B])
    logFC = mean_B - mean_A
    
    if( (mean_A > thr.filter) | (mean_B > thr.filter)){
      obj<-try(t.test(x[col_A], x[col_B],var.equal = TRUE), silent=TRUE)
      if (is(obj, "try-error")){p.value =1}	else {p.value = obj$p.value}
      
      # p.value = t.test(x[col_A], x[col_B], var.equal = T)$p.value
      if( p.value < thr.p ){
        DEGtype='not DEG'
        if( logFC > log2(thr.FC) ) DEGtype='Up-regulated'
        if( logFC < (-log2(thr.FC)) ) DEGtype='Down-regulated'
      } else {
        DEGtype='not DEG'
      }
    } else {
      p.value = NA
      DEGtype= 'Low-expressed'
    }
    
    
    return(list(
      mean = mean, 
      sd = sd,
      mean_A = mean_A,
      mean_B = mean_B,
      logFC = logFC,
      p.value = p.value,
      DEGtype = DEGtype
    ))
    
  }) %>% rbindlist()
  
  dt.DEG$group_A = group_A
  dt.DEG$group_B = group_B
  dt.DEG$gene = rownames(exprMat_filt)
  
  
  dt.DEG.passLowFlit <- dt.DEG[DEGtype!='Low-expressed']
  dt.DEG.passLowFlit$p.value.adj <- p.adjust(dt.DEG.passLowFlit$p.value, method='fdr')
  dt.DEG.passLowFilt.sorted <- dt.DEG.passLowFlit[order(ifelse(DEGtype%in%c('Up-regulated','Down-regulated'),0,1), 
                                                        -abs(logFC))]
  
  dt.DEG.lowExpr <- dt.DEG[DEGtype=='Low-expressed']
  
  dt.DEG <- rbindlist(list(dt.DEG.passLowFilt.sorted, dt.DEG.lowExpr), use.names = T, fill = T)
  
  
  dt.DEG.forOutput <-  dt.DEG[,.(group_A, group_B, gene, mean, mean_A, mean_B, sd, logFC, p.value, p.value.adj, DEGtype)]
  return(dt.DEG.forOutput)
  
}


