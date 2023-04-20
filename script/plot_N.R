
source("library.R")
source("color_theme.R")
source("func_snrdb.R")



meta<-read.csv("meta_data.csv")

sample12<-unique(meta[,c("sample_rep","sample")])

pasteto<-function(x){
  y<-c()
  for (i in 1:length(x)){
    y<-paste(y,x[i],sep=";")
    
  }
  y<-gsub("^;","",y)
  return(y)
}


#---------------------------------#
###########snr_intra###############
#---------------------------------#



#import
snr_intra<-readRDS("intraN_SNR.rds")
snr_intra$sgr<-""

for ( i in 1:nrow(snr_intra)){
  x<-snr_intra$detail[i]
  
 x1<-sort(unique(sample12$sample[match( unlist(strsplit(x,";")),sample12$sample_rep)]))
  
 snr_intra$sgr[i]<-pasteto(x1)
}



snr_intra$groups<-gsub("G3R2","G3R2*",snr_intra$groups,fix=T)
snr_intra$groups<-gsub("G4R2","G4R2*",snr_intra$groups,fix=T)
snr_intra$groups<-gsub("G3R3","G3R3*",snr_intra$groups,fix=T)
snr_intra$groups<-gsub("G4R3","G4R3*",snr_intra$groups,fix=T)


ll<-names(table(snr_intra$n_sample))
snr_intra$n_sample<-factor(snr_intra$n_sample,levels=ll,ordered=T)

x<-unique(snr_intra[,c("n_sample","groups")])
ll1<-x[order(x$n_sample),"groups"]
snr_intra$groups<-factor(snr_intra$groups,levels=ll1,ordered=T)

p_intra_snr<-ggplot(snr_intra,aes(x=sgr,y=snr,fill=n_sample))+
  geom_boxplot(color="black",alpha=0.9)+
  theme_few()+
  mytheme12 +
  myx90+
  scale_fill_brewer(palette = "YlGnBu",name="N library") +
  facet_grid(.~groups,scales="free",space="free_x")+
  labs(x="",y="SNR")

p_intra_snr



#---------------------------------#
###########RC_intra###############
#---------------------------------#

rc_intra<-readRDS("intraN_RC.rds")

rc_intra$sgr<-""

for ( i in 1:nrow(rc_intra)){
  x<-rc_intra$detail[i]
  
  x1<-sort(unique(sample12$sample[match( unlist(strsplit(x,";")),sample12$sample_rep)]))
  
  rc_intra$sgr[i]<-pasteto(x1)
}



rc_intra$groups<-gsub("R2","R2*",rc_intra$groups,fix=T)
rc_intra$groups<-gsub("R3","R3*",rc_intra$groups,fix=T)


ll<-names(table(rc_intra$n_sample))
rc_intra$n_sample<-factor(rc_intra$n_sample,levels=ll,ordered=T)


x<-unique(rc_intra[,c("n_rep","groups")])
ll1<-x[order(x$n_rep),"groups"]
rc_intra$groups<-factor(rc_intra$groups,levels=ll1,ordered=T)

p_intra_rc<-ggplot(rc_intra,aes(x=sgr,y=RC,fill=n_sample))+
  geom_boxplot(color="black",alpha=0.9,outlier.colour = NA)+
  theme_few()+
  mytheme12 +
  myx90+
  scale_fill_brewer(palette = "YlGnBu",name="N library") +
  facet_grid(.~groups,scales="free",space="free_x")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1,size=12))+
  labs(x="",y="RC")

p_intra_rc



#---------------------------------#
###########N ref SNR###############
#---------------------------------#

refN<-readRDS("refN_SNR.rds")
colnames(refN)[colnames(refN)=="group"]<-"groups"

refN$sgr<-""

for ( i in 1:nrow(refN)){
  x<-refN$detail[i]
  
  x1<-sort(unique(sample12$sample[match( unlist(strsplit(x,";")),sample12$sample_rep)]))
  
  refN$sgr[i]<-pasteto(x1)
}


refN$groups<-paste(refN$groups,"*",sep="")
refN$groups<-gsub("G0R0*","/",refN$groups,fix=T)
refN$groups<-gsub("G1R1*","G1R1",refN$groups,fix=T)
refN$groups<-gsub("G1R2*","G1R2",refN$groups,fix=T)
refN$sgr<-gsub("NA;","Without \n correction",refN$sgr)

x<-unique(refN[,c("n_sample","groups")])
ll1<-x[order(x$n_sample),"groups"]
refN$groups<-factor(refN$groups,levels=ll1,ordered=T)



p2<-ggplot(refN,aes(x=sgr,y=SNR,fill=factor(n_sample)))+
  geom_boxplot(color="black",alpha=0.9)+
  theme_few()+
  mytheme12 +
  myx90+
  scale_fill_brewer(palette = "YlGnBu",name="N library") +
  facet_grid(.~groups,scales="free",space="free_x")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1,size=12))+
  labs(x="",y="SNR")

p2



g<-plot_grid(p_intra_snr,p_intra_rc,p2,ncol=1,align="hv",axis="tblr")

g
