############################################

source("library.R")
source("color_theme.R")


####################import data

logexpr<-readRDS("exprMat_log2FPKM.rds")
meta<-read.csv("meta_data.csv")
samplepair<-readRDS("samplepairs_r526.rds")


ubatch<-unique(as.character(meta$batch))
usample<-as.character(unique(meta$sample))

samplepair1<-samplepair[which(samplepair$n_group>1),]
samplepair1<-samplepair1[grep("D6",samplepair1$detail),]

samplepair1$detail<-as.character(samplepair1$detail)

samplepair2<-samplepair1

intra_N_RC<-c()

for ( i in 1:length(ubatch)){
  print(i)
  meta_batch<-meta[meta$batch==ubatch[i],]
  
  for ( j in 1:nrow(samplepair2)){
    meta_batch_sam<-meta_batch[meta_batch$sample_rep %in% unlist(strsplit(samplepair2[j,"detail"],";")),]
    
    #D6lib
    D6.lib<-meta_batch_sam$library[meta_batch_sam$sample=="D6"]
    if(length(D6.lib)==1){
      mat.D6<-logexpr[,colnames(logexpr)%in% D6.lib]
    }else{
      mat.D6<-rowMeans(logexpr[,colnames(logexpr)%in% D6.lib])
    }
    
    
    #not D6
    other.lib<-meta_batch_sam$library[!meta_batch_sam$sample=="D6"]
    uother.sample<-unique(meta_batch_sam$sample[match(other.lib,meta_batch_sam$library)])
    
    mat.other.mean.full<-c()
    for ( k in 1:length(uother.sample)){
      
      lib<-other.lib[grep(uother.sample[k],other.lib)]
      
      if(length(lib)==1){
        mat.other.mean<-data.frame(
          mean=logexpr[,colnames(logexpr)%in% lib],
          sample=rep(uother.sample[k],nrow(logexpr)),
          gene=rownames(logexpr))
        
      }else{
        mat.other.mean<-data.frame(
          mean=rowMeans(logexpr[,colnames(logexpr)%in% lib]),
          sample=rep(uother.sample[k],nrow(logexpr)),
          gene=rownames(logexpr))
      }
      mat.other.mean.full<-rbind(mat.other.mean.full,mat.other.mean)
    }
    
    mat.fc.full<-mat.other.mean.full
    mat.fc.full$D6<-mat.D6[match(mat.fc.full$gene,names(mat.D6))]   
    mat.fc.full$FC<-mat.fc.full$mean-mat.fc.full$D6
    mat.fc.full$compare<-paste(mat.fc.full$sample,"/D6",sep="")
    mat.fc.full$gene_compare<-paste(mat.fc.full$gene,mat.fc.full$compare)
    
    int<-intersect(mat.fc.full$gene_compare,refexpr$gene_compare)
    
    testlab_logFC_f<-mat.fc.full[match(int,mat.fc.full$gene_compare),]
    Ref_RE_f<-refexpr[match(int,refexpr$gene_compare),]
    colnames(Ref_RE_f)<-paste("ref_",colnames(Ref_RE_f),sep="")
    colnames(testlab_logFC_f)<-paste("testlab_",colnames(testlab_logFC_f),sep="")
    
    
    refvstest<-merge(Ref_RE_f[,c("ref_gene_compare","ref_meanlogFC")],testlab_logFC_f[,c("testlab_gene_compare","testlab_FC")],by.x="ref_gene_compare",by.y="testlab_gene_compare")
    
    c1<-cor(refvstest$ref_meanlogFC,refvstest$testlab_FC,use="pairwise.complete.obs")
    intra_N_RC<-rbind(intra_N_RC,c(ubatch[i],samplepair2[j,],c1))
    
    if(j %in% seq(0,nrow(samplepair2),10)) {
      print(j)
    }
    
  }
  
}


colnames(intra_N_RC)<-c("batch","n_sample","detail","n_group","n_rep","groups","RC")

head(intra_N_RC)

intra_N_RC<-data.frame(intra_N_RC)

for ( i in 1:ncol(intra_N_RC)){
  intra_N_RC[,i]<-(as.character(intra_N_RC[,i]))
  
}

for ( i in c(2,4,5,7)){
  intra_N_RC[,i]<-as.numeric(as.character(intra_N_RC[,i]))
  
}

saveRDS(intra_N_RC,"intraN_RC.rds")

