
source("library.R")
source("func_snrdb.R")


logexpr<-readRDS("exprMat_log2FPKM.rds")

meta<-read.csv("meta_data.csv")

ubatch<-unique(as.character(meta$batch))


pcs_all<-c()
snr_table<-c()
for ( i in 1:length(ubatch)){
  name<-as.character(meta$library[meta$batch==ubatch[i]])
  
  logexpr_batch<-logexpr[,colnames(logexpr)%in% name]
  expdgene1_batch<-expdgene1[expdgene1$batch==ubatch[i],]
  gg1<-names(table(expdgene1_batch$gene)==4)
  
  mat<-logexpr_batch[rownames(logexpr_batch)%in% gg1,]
  
  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]
  
  sample<-meta$sample[match(colnames(dat_z),meta$library)]
  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)
  snr_table<-rbind(snr_table,c(ubatch[i],snrvalue,length(gg1)))
  
  pca_prcomp <- prcomp(t(dat_z ), scale=F,center =F)
  pcs <- data.frame(predict(pca_prcomp)) 
  
  pcs<-pcs[,c(1:2)]
  
  pcs$PC1_comp<-round(summary(pca_prcomp)$importance[2,1]*100)
  pcs$PC2_comp<-round(summary(pca_prcomp)$importance[2,2]*100)
  
  pcs$library<-rownames(pcs)
  
  pcs$SNR<-snrvalue
  
  pcs$sample<-meta$sample[match(pcs$library,meta$library)]
  pcs$batch<-ubatch[i]
  
  pcs_all<-rbind(pcs_all,pcs)
  
  print(i)
  
}


saveRDS(pcs_all,"intra_PCA.rds")





