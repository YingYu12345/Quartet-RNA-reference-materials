source("library.R")
source("color_theme.R")


library(scales)
library(ggsci)
library(RColorBrewer) 
my_pal<-pal_jama(alpha =1)(8)
show_col(my_pal)


########SNR 11

SNR11<-readRDS("expr_mat/2b_intra_SNR11_20230406.rds")
SNR12<-readRDS("expr_mat/2b_intra_SNR_20230406.rds")

###################plot

xx<-tapply(SNR11$SNR,SNR11$Batch,function(x){mean(x)+6})
SNR11$cut<-round(xx[match(SNR11$Batch,names(xx))],1)

SNR11$out<-ifelse(SNR11$SNR> SNR11$cut,"outside","inside")
SNR11$ExdLibrary_1<-sapply(strsplit(as.character(SNR11$ExdLibrary),"_"),function(x){paste(x[5],x[6],sep="_")})

SNR11$out<-factor(SNR11$out,levels=c("outside","inside"),ordered=T)
colnames(SNR11)[2]<-"SNR11"

ubatch<-unique(as.character(SNR11$Batch))
SNR11max<-c()
for ( i in 1:length(ubatch)){
  x<-SNR11[SNR11$Batch==ubatch[i],]
  SNR11max<-rbind(SNR11max,x[x$SNR==max(x$SNR),])
}

SNR11max$SNR12<-SNR12$SNR[match(SNR11max$Batch,SNR12$Batch)]

#plot bar
p<-ggplot(data=SNR12,aes(x=reorder(Batch,-SNR),y=SNR))+
   geom_bar(position=position_dodge(),stat="identity",fill= my_pal[3])+
   theme_few()+
   mytheme12+
   labs(y="SNR",x="")+
   geom_hline(yintercept = 12,linetype=2)+
   theme(legend.text = element_text(size=12),
         axis.ticks.x = element_blank())+
   myx90+
   scale_y_continuous(limits =c(0,35),breaks=c(0,10,12,20,30,35))

p

point.colors <- c( "Normal"="#DF8F44FF", "Outlier"="#B24745FF")
point.shape<-c("Normal"=16,"Outlier"=16)
#point

p1<-p+
	geom_point(data=subset(SNR11,out=="inside"),aes(x=Batch,y=SNR,color="Normal",shape="Normal"),size=2)+
	geom_point (data=subset(SNR11max,out=="outside"),aes(x=Batch,y=SNR,color="Outlier",shape="Outlier"),size = 3.5) +#max+
	scale_colour_manual(name="SNR11",values=point.colors)+
	scale_shape_manual(name="SNR11",values=point.shape)+
	geom_text_repel(data=subset(SNR11,out=="outside"),aes(x=Batch,y=SNR,label = ExdLibrary_1 ),size = 3.5) 


p1

