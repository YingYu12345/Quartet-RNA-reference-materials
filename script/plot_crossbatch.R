
source("color_theme.R")
source("func_snrdb.R")

####################import data

logexpr252<-readRDS("alldetect_logFPKM.rds")
ratio_expr252<-readRDS("alldetect_ratioD6_f001_r29433c252.rds")


#############import metadata 
meta<-read.csv("meta_data.csv")

ubatch<-sort(unique(as.character(meta$batch)),decreasing=T)



####################################raw pca _ribo

ribo_raw<-logexpr252[,grep("^R_",colnames(logexpr252))]
ribo_ratio<-ratio_expr252[,grep("^R_",colnames(ratio_expr252))]
poly_raw<-logexpr252[,grep("^P_",colnames(logexpr252))]
poly_ratio<-ratio_expr252[,grep("^P_",colnames(ratio_expr252))]
both_raw<-logexpr252
both_ratio<-ratio_expr252

datnam<-c("ribo_raw","poly_raw","both_raw","poly_ratio","ribo_ratio","both_ratio")

pcsall<-c()
for ( i in 1:6){
  mat<-get(datnam[i])
  
  dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]
  
  sample<-meta$sample[match(colnames(dat_z),meta$library)]
  snrvalue<- round(snrdb_function(dat_z,as.factor(sample)),1)
  
  pca_prcomp <- prcomp(t(dat_z ),scale=F,retx = T)
  pcs <- predict(pca_prcomp) %>% data.frame()
  pcs<-pcs[,1:2]
  pcs$sample = meta$sample[match(rownames(pcs),meta$library)]
  pcs$batch = meta$batch[match(rownames(pcs),meta$library)]
  
  pcs$datatype<-datnam[i]
  pcs$SNR<-snrvalue
  
  PC1_comp<-round(summary(pca_prcomp)$importance[2,1]*100)
  PC2_comp<-round(summary(pca_prcomp)$importance[2,2]*100)
  
  pcs$PC1_comp<-PC1_comp
  pcs$PC2_comp<-PC2_comp
  
  
  pcsall<-rbind(pcsall,pcs)
  
  print(i)
}

pcsall$PC2[pcsall$datatype=="poly_ratio"]<-(-pcsall$PC2[pcsall$datatype=="poly_ratio"])

nam<-c("poly_raw","ribo_raw","both_raw","poly_ratio","ribo_ratio","both_ratio")

pcsall$datatype<-factor(pcsall$datatype,levels=nam,ordered=T)


for(i in 1:length(nam)){
  
  pcs<-pcsall[pcsall$datatype==nam[i],]
  xlab = paste("PC1 (",pcs$PC1_comp[1],"%)",sep="")
  ylab = paste("PC2 (",pcs$PC2_comp[1],"%)",sep="")
  tit<-namout[i]
  
  p<-ggplot(pcs, aes(x=PC1, y=PC2))+
  geom_point(aes(color=sample,shape=batch),size=2.5)+
  theme_few()+
  mytheme12+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=9),
        legend.position = "bottom",
        plot.title=element_text(hjust=0.5,size=13),
        plot.subtitle = element_text(hjust=0.5,size=12)
        )+
  scale_color_manual(values = colors.sample.Quartet.fill,name="Group") +
  scale_fill_manual(values = colors.sample.Quartet.fill,name="Group") +
  scale_shape_manual(values=shape_batch,name="Batch")+
  labs(x = xlab,y = ylab,title=tit)+
  guides(colour = guide_legend(override.aes = list(size=3),ncol=2))+
  guides(shape=guide_legend(override.aes = list(size=2),ncol=7))

p

p_n<-p+
  theme(legend.position = "none",
        plot.margin = margin(12, 12, 12, 12))

eval(parse(text=paste("p_",i,"<-p_n",sep="")))



}

p3<-plot_grid(p_1,p_2,p_3,p_4,p_5,p_6,ncol=3)

p3