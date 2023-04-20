source("library.R")
source("func_snrdb.R")

logexpr<-readRDS("exprMat_log2FPKM.rds")

meta<-readRDS("data/metaQ_r252_20230406.rds")

ubatch<-unique(as.character(meta$batch))
usample<-as.character(unique(meta$sample))

samplepair<-readRDS("samplepairs_r526.rds")

################calculate SNRs 

samplepair$SNR<-0


for ( s in 1:nrow(samplepair)){
  if(samplepair$n_sample[s]==0){
    ratio_expr<-logexpr
  }else{
    
    ratio_expr<-c()
    for ( i in 1:length(ubatch)){
      
      meta_batch<-meta[meta$batch==ubatch[i],]
      name<-meta_batch$library
      logexpr_batch<-logexpr[,colnames(logexpr)%in% name]
      
      sample_m<-unlist(strsplit(as.character(samplepair$detail[s]),";"))
      libs_m<-meta_batch$library[meta_batch$sample_rep %in% sample_m]
      
      if(length(libs_m)>1){
        m<-rowMeans(logexpr_batch[,colnames(logexpr_batch)%in% libs_m])
      }else{
        m<-logexpr_batch[,colnames(logexpr_batch)%in% libs_m]
      }
      mat<-apply(logexpr_batch,2,function(x){x-m})
      ratio_expr<-cbind(ratio_expr,mat) 
    }
  }
  
  mat<-ratio_expr
  
  mat<-mat[apply(mat,1,function(x){length(which(is.na(x)))==0}),]
  mat<-mat[apply(mat,1,function(x){length(which(is.infinite(x)))==0}),]
  
  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]
  
  sample<-meta$sample[match(colnames(dat_z),meta$library)]
  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)
  
  mat<-c()
  samplepair$SNR[s]<-snrvalue
  
print(s)
  
} 


