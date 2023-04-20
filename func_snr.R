################function
#usage: 
##note: expression matrix should be z-scaled before using snrdb_function
#  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
#  sample<-meta_batch$sample[match(colnames(dat_z),meta_batch$library)]
#  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)


#######
snrdb_function<-function(exprMat,group){
  
  library(data.table)
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),center=F,scale=F)
#  pca_prcomp <- prcomp(t(exprMat),retx=T,center=T,scale=F)
  
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  

  signoise_db <- 10*log10(signoise)
  return(signoise_db)
  
}
