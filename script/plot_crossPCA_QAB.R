
source("library.R")
source("color_theme.R")

#######################


log2f<-readRDS("exprMatQAB_logFPKM_f001_B20.rds")
meta2<-readRDS("metaQAB_r360.rds")
ubatch<-unique(as.character(meta2$batch))


log2f<-log2f[,!colnames(log2f)%in% "R_ILM_L2_B1_B_3"]
meta2<-meta2[!meta2$library %in% "R_ILM_L2_B1_B_3", ]


#########
up<-unique(meta2$protocol)

mins<-min(log2f)

log2f_both<-log2f
log2f_p<-log2f[,colnames(log2f)%in% meta2$library[meta2$protocol=="P"]]
log2f_r<-log2f[,colnames(log2f)%in% meta2$library[meta2$protocol=="R"]]

nam<-c("log2f_both","log2f_p","log2f_r")
namout<-c("Both","P","R")

pcsall<-c()
for ( i in 1:length(nam)){
    
  logexpr_1<-get(nam[i])
  mat<-logexpr_1[apply(logexpr_1,1,function(x){length(which(x>mins))>0.5*ncol(logexpr_1)}),]
  
  dat_z2<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
  dat_z2<-dat_z2[apply(dat_z2,1,function(x){length(which(is.na(x)))==0}),]
  
  sample<-meta2$sample[match(colnames(dat_z2),meta2$library)]
  
  ######pca
  pca_prcomp2 <- prcomp(t(dat_z2 ), scale=F,center = F)
  pcs <- data.frame(predict(pca_prcomp2))
  
  pcs<-pcs[,c("PC1","PC2")]
  pcs$library<-rownames(pcs)
  
  pcs2<-merge(pcs,meta2,by.x="library",by.y="library")
  
  pcs2$cal<-namout[i]
  
  PC1_comp<-round(summary(pca_prcomp2)$importance[2,1]*100)
  PC2_comp<-round(summary(pca_prcomp2)$importance[2,2]*100)
  
  pcs2$PC1_comp<-PC1_comp
  pcs2$PC2_comp<-PC2_comp
  
  pcsall<-rbind(pcsall,pcs2)
  
  print(i)
  
}


up<-unique(pcsall$protocol)
ucal<-unique(pcsall$cal)
ucal<-c("P","R","Both")
ucal_out<-c("Poly(A)","RiboZero","Both")

for(i in 1:length(ucal)){

  pcs<-pcsall[pcsall$cal==ucal[i],]
  xlab = paste("PC1 (",pcs$PC1_comp[1],"%)",sep="")
  ylab = paste("PC2 (",pcs$PC2_comp[1],"%)",sep="")
  tit<-ucal_out[i]
       
       
  p<-ggplot(pcs, aes(x=PC1, y=PC2))+
  geom_point(aes(color=sample,shape=batch),size=2.5)+
  theme_few()+
  mytheme12+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "right",
        plot.title=element_text(hjust=0.5))+
  scale_color_manual(values = color.sample6,name="Group") +
  scale_fill_manual(values = color.sample6,name="Group") +
  scale_shape_manual(values=shape_batch20,name="Batch")+
  labs(x = xlab,y = ylab,title=tit)+
  guides(colour = guide_legend(override.aes = list(size=3),ncol=1))+
    guides(shape="none")
    

  p
  
  p_n<-p+
  theme(legend.position = "none",
        plot.margin = margin(15, 15, 15, 15))

eval(parse(text=paste("p_",i,"<-p_n",sep="")))


}


ll<-get_legend(p)
p3<-plot_grid(p_1,p_2,ll,ncol=3,rel_widths = c(3,3,1))
