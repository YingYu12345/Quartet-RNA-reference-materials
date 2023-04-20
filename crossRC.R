

source("library.R")

logexpr<-readRDS("logFPKM.rds")
ratioexpr<-readRDS("ratioD6.rds")

meta<-readRDS("metaQ_r252.rds")

ubatch<-sort(unique(as.character(meta$batch)))


combn_s<-data.frame(
  sampleA=rep("D6",3),
  sampleB=c("D5","F7","M8")
)

colnames(combn_s)<-c("sampleA","sampleB")
combn_s$compare<-paste(combn_s$sampleB,combn_s$sampleA,sep="/")


datnam<-c("logexpr","ratioexpr")
datatype<-c("before","after")

cross_RC<-c()
for (a in 1:2){
  print(datnam[a])
  exprMat<-get(datnam[a])

  for ( i in 1:nrow(batchpairs)){
print(i)
  for(j in 1:nrow(combn_s)){
    
   lib.A=meta$library[intersect(which(meta$batch==batchpairs$batchA[i]),which(meta$sample==combn_s$sampleA[j]))]
   lib.B=meta$library[intersect(which(meta$batch==batchpairs$batchB[i]),which(meta$sample==combn_s$sampleB[j]))]
    
   exprMat_A<-rowMeans(exprMat[,colnames(exprMat) %in% lib.A])
   exprMat_B<-rowMeans(exprMat[,colnames(exprMat) %in% lib.B])
   
   FC_pairs<-data.frame(
     testFC=exprMat_B-exprMat_A,
     gene=rownames(exprMat),
     compare=rep(combn_s$compare[j],nrow(exprMat)))
   FC_pairs$gene_compare<-paste(FC_pairs$gene,FC_pairs$compare)
  
   int<-intersect(FC_pairs$gene_compare,refdata3$gene_compare)
   FC_pairs_f<-FC_pairs[FC_pairs$gene_compare %in% int,]
   refFC_pairs_f<-refdata3[refdata3$gene_compare %in% int,]
   
   merge<-merge(FC_pairs_f,refFC_pairs_f,by="gene_compare")
   
   cor.FC.pearson<-cor(merge$testFC,merge$meanlogFC,method="pearson")
   
   cross_RC<-rbind(cross_RC,c(datatype[a],combn_s$compare[j], as.character(batchpairs[i,c(3,4)]),nrow(merge),cor.FC.pearson))
    
  }
  }
  
}

cross_RC<-data.frame(cross_RC)
colnames(cross_RC)<-c("datatype","compare","batch2","batch_type","N","RC")

cross_RC$RC<-as.numeric(as.character(cross_RC$RC))

saveRDS(cross_RC,"expr_mat/RC_cross.rds")

