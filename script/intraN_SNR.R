
source("func_snrdb.R")


####################import data

logexpr<-readRDS("exprMat_log2FPKM.rds")
meta<-read.csv("meta_data.csv")

ubatch<-unique(as.character(meta$batch))
usample<-as.character(unique(meta$sample))

samplepair<-readRDS("expr_mat/samplepairs_r526_20230410.rds")
samplepair1<-samplepair[intersect(which(samplepair$n_group>1),which(samplepair$n_rep>1)),]

samplepair1$detail<-as.character(samplepair1$detail)

#----------------------------------------------
############ intra-batch N SNR-----------------
#----------------------------------------------

snr_table<-c()
  
for ( i in 1:length(ubatch)){
  print(ubatch[i])  
  
meta_batch<-meta[meta$batch==ubatch[i],]

for ( j in 1:nrow(samplepair1)){
  meta_batch_sam<-meta_batch[meta_batch$sample_rep %in% unlist(strsplit(samplepair1[j,"detail"],";")),]
  
  name<-meta_batch_sam$library
  mat<-logexpr[,colnames(logexpr)%in% name]
  
  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]

  sample<-meta_batch_sam$sample[match(colnames(dat_z),meta_batch_sam$library)]
  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)
  snr_table<-rbind(snr_table,c(ubatch[i],as.character(samplepair1[j,]),snrvalue))

if(j %in% seq(0,nrow(samplepair1),10)) {
  print(j)
}
    }
  }
  
colnames(snr_table)<-c("batch","n_sample","detail","n_group","n_rep","groups","snr")

head(snr_table)

snr_table<-data.frame(snr_table)

for ( i in 1:ncol(snr_table)){
  snr_table[,i]<-(as.character(snr_table[,i]))
  
}

for ( i in c(2,4,5,7)){
  snr_table[,i]<-as.numeric(as.character(snr_table[,i]))
  
}

saveRDS(snr_table,"intraN_SNR.rds")

