
source("library.R")
source("func_snrdb.R")

logexpr<-readRDS("")

meta<-readRDS("")
ubatch<-unique(as.character(meta$batch))

expdgene<-readRDS("")
expdgene1<-expdgene

########

snr_table<-c()
for ( i in 1:length(ubatch)){
  expdgene1_batch<-expdgene1[expdgene1$batch==ubatch[i],]
  gg1<-names(table(expdgene1_batch$gene)==4)
  logexpr_batch<-logexpr[rownames(logexpr)%in% gg1,]
  
  name<-as.character(meta$library[meta$batch==ubatch[i]])
  
  for (j in 1:12){
  
  name1<-setdiff(name,name[j])
  
  mat<-logexpr_batch[,colnames(logexpr_batch)%in% name1]

  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]
  
  sample<-meta$sample[match(colnames(dat_z),meta$library)]
  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)
  snr_table<-rbind(snr_table,c(ubatch[i],snrvalue,length(gg1),name[j]))
  }
  print(i)
  
  }
  

snr_table<-data.frame(snr_table)
colnames(snr_table)<-c("Batch","SNR","N","ExdLibrary")
snr_table$SNR<-as.numeric(as.character(snr_table$SNR))
snr_table$N<-as.numeric(as.character(snr_table$N))

saveRDS(snr_table,"expr_mat/intra_SNR11.rds")


