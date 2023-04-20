source("library.R")

pvalue<-function(x,group){
      obj<-try(t.test(x[group==levels(group)[1]], x[group==levels(group)[2]],var.equal = TRUE), silent=TRUE)
    if (is(obj, "try-error")){p.value =1}	else {p.value = obj$p.value}
 
    return(p.value)   
}

####################import data

logexpr<-readRDS()
ratioexpr<-readRDS()

meta<-readRDS()

sample_combn<-data.frame(
  sampleA=rep("D6",3),
  sampleB=c("D5","F7","M8")
)

sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")

usample<-unique(as.character(meta$sample))

expd_refs<-c()
for ( i in 1:length(usample)){
  
  x<-expdgene1_p[expdgene1_p$sample==usample[i],]
  expd_sample<-data.frame(table(x$gene))
  expd_sample<-cbind(expd_sample,usample[i])
  expd_refs<-rbind(expd_refs,expd_sample)
}

expd_refs<-data.frame(expd_refs)
colnames(expd_refs)<-c("gene","Freq","sample")
expd_refs$gene<-as.character(expd_refs$gene)
expd_refs$sample<-as.character(expd_refs$sample)


expd_refs_f<-expd_refs[expd_refs$Freq==13,]

alldet<-names(which(table(expd_refs_f$gene)==4))

detect_genes_pairs<-c()
for (i in 1:nrow(sample_combn)){
  s<-intersect(expd_refs_f$gene[expd_refs_f$sample==sample_combn$sampleA[i]],expd_refs_f$gene[expd_refs_f$sample==sample_combn$sampleB[i]])
  detect_genes_pairs<-rbind(detect_genes_pairs,cbind(s,as.character(sample_combn$compare[i])))
}
colnames(detect_genes_pairs)<-c("gene","compare")
detect_genes_pairs<-data.frame(detect_genes_pairs)

detect_genes_pairs$gene_compare<-paste(detect_genes_pairs$gene,detect_genes_pairs$compare)

DEGs_p<-DEGs[DEGs$batch %in% passbatch,]
DEGs_p$protocol<-sapply(strsplit(as.character(DEGs_p$batch),"_"),function(x){x[1]})
DEGs_p$gene_compare<-paste(DEGs_p$gene,DEGs_p$compare)

DEGs_p_dec<-DEGs_p[DEGs_p$gene_compare %in%  detect_genes_pairs$gene_compare,]

#####filter p<0.05

DEGs_p_f<-DEGs_p_dec[DEGs_p_dec$`P.Value`<0.05,]

DEGs_p_cal<-data.frame(table(paste(DEGs_p_f$gene,DEGs_p_f$compare)))
DEGs_p_cal$gene<-sapply(strsplit(as.character(DEGs_p_cal$Var1)," "),function(x){x[1]})
DEGs_p_cal$compare<-sapply(strsplit(as.character(DEGs_p_cal$Var1)," "),function(x){x[2]})

DEGs_p_cal_f<-DEGs_p_cal[DEGs_p_cal$Freq>=4,]
DEGs_p_cal_f_fd<-DEGs_p_cal_f[DEGs_p_cal_f$Var1 %in% detect_genes_pairs$gene_compare,]

DEGs_p_cal_f_fd<-DEGs_p_cal_f_fd[order(DEGs_p_cal_f_fd$compare),]
DEGs_p_cal_f_fd$gene_compare<-paste(DEGs_p_cal_f_fd$gene,DEGs_p_cal_f_fd$compare)



#########################filtering protocol-dependent genes
prot_fc<-c()
for ( i in 1:nrow(DEGs_p_cal_f_fd2)){
  
  c=as.character(DEGs_p_cal_f_fd2$compare)[i]
  g=as.character(DEGs_p_cal_f_fd2$gene)[i]
  m<-DEGs_p[intersect(which(DEGs_p$gene==g),which(DEGs_p$compare==c)),]
  
  if(length(which(grepl("P",m$protocol)))>=3 && length(which(grepl("R",m$protocol)))>=3){
    p<-pvalue(m$logFC,as.factor(m$protocol))
    f<-tapply(m$logFC,as.factor(m$protocol),mean)
    fc<-f[1]-f[2]
    prot_fc<-rbind(prot_fc,c(g,c,p,f,fc))
  }
    
  if(i %in% seq(0,nrow(DEGs_p_cal_f_fd2),by=200)){
print(i)
  }
  
}

colnames(prot_fc)<-c("gene","compare","pvalue","mean_polyA","mean_riboZ","FC_P2R")
prot_fc<-data.frame(prot_fc)
for ( i in 1:ncol(prot_fc)){
  prot_fc[,i]<-as.character(prot_fc[,i])
  
}

for ( i in 3:ncol(prot_fc)){
  prot_fc[,i]<-as.numeric(prot_fc[,i])
  
}

prot_fc$gene_compare<-paste(prot_fc$gene,prot_fc$compare)
ff<-prot_fc[intersect(which(prot_fc$pvalue<0.05),which(abs(prot_fc$FC_P2R)>=1)),]


#####filter 
ref_FC_f2<-DEGs_p_cal_f_fd2[!as.character(DEGs_p_cal_f_fd2$gene_compare) %in% as.character(ff$gene_compare),]




DEGs_p<-DEGs[DEGs$batch %in% passbatch,]
DEGs_p$gene_compare<-paste(DEGs_p$gene,DEGs_p$compare)

#function
rse <- function(x){sd(x, na.rm = T)/(mean(x, na.rm=T)*sqrt(length(x[!is.na(x)])))}
gene_compare=levels(as.factor(DEGs_p$gene_compare))

rse.DEG<-data.frame(
  rseFC=tapply(DEGs_p$logFC,as.factor(DEGs_p$gene_compare),function(x){rse(2^x)*100}),
  meanlogFC=tapply(DEGs_p$logFC,as.factor(DEGs_p$gene_compare),mean),
  medianp=tapply(DEGs_p$P.Value,as.factor(DEGs_p$gene_compare),median),
  gene_compare=levels(as.factor(DEGs_p$gene_compare))
)

ref_FC_f2$RSE.FC<-rse.DEG$rseFC[match(ref_FC_f2$gene_compare,rse.DEG$gene_compare)]
ref_FC_f2$meanlogFC<-rse.DEG$meanlogFC[match(ref_FC_f2$gene_compare,rse.DEG$gene_compare)]
ref_FC_f2$FC<-2^ref_FC_f2$meanlogFC
ref_FC_f2$medianp<-rse.DEG$medianp[match(ref_FC_f2$gene_compare,rse.DEG$gene_compare)]


###################################
#########reference DEG

DEGs_p<-DEGs[DEGs$batch %in% passbatch,]
DEGs_p$protocol<-sapply(strsplit(as.character(DEGs_p$batch),"_"),function(x){x[1]})
DEGs_p$gene_compare<-paste(DEGs_p$gene,DEGs_p$compare)


DEGs_p$DEG_type<-"non-DEG"
DEGs_p$DEG_type[intersect(which(DEGs_p$P.Value<0.05),which(DEGs_p$logFC>=1))]<-"up-regulate"
DEGs_p$DEG_type[intersect(which(DEGs_p$P.Value<0.05),which(DEGs_p$logFC<=(-1)))]<-"down-regulate"

########DEGs_cal 

DEG_ref_cal<-data.frame(
  N_up=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="up-regulate"))}),
  N_non=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="non-DEG"))}),
  N_down=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="down-regulate"))}))

DEG_ref_cal$Final<-"non-DEG"
DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_up>=2),which(DEG_ref_cal$N_down>=2))]<-"conflicting"

DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_up>=4),which(DEG_ref_cal$N_down==0))]<-"up-regulate"
DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_down>=4),which(DEG_ref_cal$N_up==0))]<-"down-regulate"
DEG_ref_cal$gene<-sapply(strsplit(rownames(DEG_ref_cal)," "),function(x){x[1]})
DEG_ref_cal$compare<-sapply(strsplit(rownames(DEG_ref_cal)," "),function(x){x[2]})

DEG_ref_cal$Var1<-rownames(DEG_ref_cal)

DEG_ref_cal_f<-DEG_ref_cal[DEG_ref_cal$Var1 %in%refFC_final$gene_compare, ]

refFC_final$N_up<-DEG_ref_cal$N_up[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$N_non<-DEG_ref_cal$N_non[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$N_down<-DEG_ref_cal$N_down[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$Final<-DEG_ref_cal$Final[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]

DEG_ref<-refFC_final[grep("regulate",refFC_final$Final),]













