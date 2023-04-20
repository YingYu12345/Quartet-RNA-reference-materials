

#############
source("library.R")
source("theme_color.R")
source("color_theme.R")
source("func_snrdb.R")



###########DEG 

DEGanalysis <- function(exprMat, group){
  
  library(edgeR)
  library(limma)
  
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)
  
  keep <- filterByExpr(dge, design,min.count = 0)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  
  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$gene = rownames(result)
  result$groupA =  levels(group)[1]
  result$groupB =  levels(group)[2]
  return(as.data.table(result))
}


#############import 
meta<-read.csv("meta_data.csv")

exprMat_fpkm<-readRDS("exprMat_FPKM.rds")
exprMat_count<-readRDS("exprMat_count.rds")

logexpr<-apply(exprMat_fpkm,2,function(x){log2(x+0.01)})

usample<-c("D5","D6","F7","M8")
sample_combn<-data.frame(
  sampleA=rep("D6",3),
  sampleB= c("D5","F7","M8"))

colnames(sample_combn)<-c("sampleA","sampleB")
sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")

ubatch<-sort(unique(as.character(meta$batch)))

for( i in 1:ncol(sample_combn)){
  sample_combn[,i]<-as.character(sample_combn[,i])
}



DEG.cal.D56<-c()
for ( i in 1:length(ubatch)){
  meta_dat<-meta[meta$batch==ubatch[i],]
  expr<-exprMat_count[,match(meta_dat$library,colnames(exprMat_count))]
  
  j=1
    lst.library.sampleA=meta_dat$library[meta_dat$sample==sample_combn$sampleA[j]]
    lst.library.sampleB=meta_dat$library[meta_dat$sample==sample_combn$sampleB[j]]
    lst.library.forDEG <- c(lst.library.sampleA, lst.library.sampleB)
    
    y.DEG <- DEGanalysis(expr[,lst.library.forDEG],
                         factor(as.character(meta_dat$sample[match(lst.library.forDEG,meta_dat$library)])))
    
    y.DEG$batch<-ubatch[i]
    y.DEG<-data.frame(y.DEG)
    y.DEG$compare<-paste(y.DEG$groupB,"/",y.DEG$groupA,sep="")
    
    y.DEG_rev<-data.frame(
      logFC= -y.DEG$logFC,
      AveExpr=y.DEG$AveExpr,
      t=y.DEG$t,
      P.Value=y.DEG$P.Value,
      adj.P.Val=y.DEG$adj.P.Val,
      B=y.DEG$B,
      gene=y.DEG$gene,
      groupA=y.DEG$groupB,
      groupB=y.DEG$groupA,
      batch=y.DEG$batch,
      compare=rep(sample_combn$compare[j],nrow(y.DEG))
    )
    
    DEG.cal.D56<-rbind(DEG.cal.D56,y.DEG_rev)
  
  print(i)
}

DEG.cal.D56<-data.frame(DEG.cal.D56)
DEG.cal.D56$Type<-"non-DEG"
DEG.cal.D56$Type[intersect(which(DEG.cal.D56$P.Value<0.05),which(DEG.cal.D56$logFC>=1))]<-"up-regulate"
DEG.cal.D56$Type[intersect(which(DEG.cal.D56$P.Value<0.05),which(DEG.cal.D56$logFC<=(-1)))]<-"down-regulate"


############F7/D6, M8/D6


DEG.cal.other<-c()
for ( i in 1:length(ubatch)){
  meta_dat<-meta[meta$batch==ubatch[i],]
  expr<-exprMat_count[,match(meta_dat$library,colnames(exprMat_count))]
  
  for (j in 2:nrow(sample_combn)){
    lst.library.sampleA=meta_dat$library[meta_dat$sample==sample_combn$sampleA[j]]
    lst.library.sampleB=meta_dat$library[meta_dat$sample==sample_combn$sampleB[j]]
    lst.library.forDEG <- c(lst.library.sampleA, lst.library.sampleB)
     
    y.DEG <- DEGanalysis(expr[,lst.library.forDEG],
                         factor(as.character(meta_dat$sample[match(lst.library.forDEG,meta_dat$library)])))
    
    y.DEG$batch<-ubatch[i]
    y.DEG$compare<-sample_combn$compare[j]
    
    DEG.cal.other<-rbind(DEG.cal.other,y.DEG)
  }
  print(i)
}

DEG.cal.other<-data.frame(DEG.cal.other)
DEG.cal.other$Type<-"non-DEG"
DEG.cal.other$Type[intersect(which(DEG.cal.other$P.Value<0.05),which(DEG.cal.other$logFC>=1))]<-"up-regulate"
DEG.cal.other$Type[intersect(which(DEG.cal.other$P.Value<0.05),which(DEG.cal.other$logFC<=(-1)))]<-"down-regulate"

DEG.cal3<-rbind(DEG.cal.D56,DEG.cal.other)
saveRDS(DEG.cal3,"DEG_limma.rds")