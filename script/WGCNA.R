####################WGCNA


setwd("C:/Users/ingri/Documents/working/quartet/quartet_20200503/RNA/RNA_revision_20230327/")
source("C:/Users/ingri/Documents/working/quartet/quartet_20200503/RNA/RNA_revision_20230327/script/1_preprocess/library.R")
source("C:/Users/ingri/Documents/working/quartet/quartet_20200503/RNA/RNA_revision_20230327/script/1_preprocess/color_theme.R")

source("script/func/func_snrdb_2023.R")



library(WGCNA)

enableWGCNAThreads()
options(stringsAsFactors = FALSE);

############function
pastetogether<-function(x){
  x<-unique(x)
  y<-""
  for ( i in 1:length(x)){
    y<-paste(y,x[i],sep=";")
  }
  y<-gsub("^;","",y)
  y<-gsub(";$","",y)
  return(y)
}


##########import=====================================================


##meta

meta<-readRDS("data/metaQ_r252_20230406.rds")

ubatch<-unique(as.character(meta$batch))
usample<-as.character(unique(meta$sample))

#######expr
ratio_D6<-readRDS("data/exprMat_RatioD6_r58395c252_20210701.rds")



## batch &name
passbatch<-readRDS("expr_mat/3a_passbatch_13.rds")
name<-as.character(meta$library[meta$batch %in% passbatch])

expdgene<-readRDS("expr_mat/detect_genelist_20230406.rds")
expdgene1<-expdgene[expdgene$Num>1,]


##################################################################
#######################ref detect genes
###################### detected== all 13 batches
####################################################################
expdgene1_p<-expdgene1[expdgene1$batch %in% passbatch,]

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


#########3alldet gene
gg<-alldet

###########preprocessing
exp<-ratio_D6[rownames(ratio_D6) %in% gg,colnames(ratio_D6)%in% name]

exp<-exp[head(order(apply(exp,1,sd),decreasing = T),10000),]

saveRDS(exp,"expr_mat/5a_WGCNA_relative/ratioexp_WGCNA_r10000c156_20230411.rds")

meta_f<-meta[match(name,meta$library),]


#### 
  datExpr<-t(exp)
  datTraits<-data.frame(meta_f[,c("sample","rep")])
  rownames(datTraits)<-meta_f$library
  datTraits$relation<-substr(datTraits$sample,1,1)
  datTraits$batch<-meta_f$batch[match(rownames(datTraits),meta_f$library)]
  datTraits$age<-ifelse(datTraits$sample %in% c("F7","M8"),1,2)
  datTraits$sex<-ifelse(datTraits$sample %in% c("F7"),1,2)
  
  
  for ( i in 1:ncol(datTraits)){
    datTraits[,i]<-as.numeric(as.factor(datTraits[,i]))
  }
  
  
 
#=====================================================================================
  #
  #  Code chunk 2
  #
  #=====================================================================================
  
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  
  
  #=====================================================================================
  #
  #  Code chunk 3
  #
  #=====================================================================================

  
  
  b<-50
  
  net = blockwiseModules(datExpr, power = 6,
                         TOMType = "unsigned", minModuleSize = b,
                         reassignThreshold = 4, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = b, nThreads = 4,
                         verbose = 3)
  
  
  #=====================================================================================
  #
  #  Code chunk 4
  #
  #=====================================================================================
  
  
  ## open a graphics window
  #sizeGrWindow(12, 9)
  ## Convert labels to colors for plotting
  #mergedColors = labels2colors(net$colors)
  ## Plot the dendrogram and the module colors underneath
  #plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  #                    "Module colors",
  #                    dendroLabels = FALSE, hang = 0.03,
  #                    addGuide = TRUE, guideHang = 0.05)
  
  
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = paste("expr_mat/5a_WGCNA_relative/02-networkConstruction-auto_N",b,"_20230411.RData",sep=""))
  
  modulegenecolors<-data.frame(cbind(moduleLabels,moduleColors))
  
  
  ####Relating modules to external clinical traits and identifying important genes
  
  
# black      blue     brown     green      grey      pink       red turquoise    yellow 
#  177      1777      1508       477      2590       133       229      2368       741 
  
  
  #=====================================================================================
  #
  #  Code chunk 2
  #
  #=====================================================================================
  
  
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  save(MEs,moduleTraitCor,moduleTraitPvalue, file = paste("expr_mat/5a_WGCNA_relative/03-networkConstruction-auto_N",b,"_20230411.RData",sep=""))
  
  
  print("networkConstruction")
  
  ################
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  
  library(RColorBrewer)
  png(paste("chart/5a_WGCNA_relative/wgcna_labeledheatmap_N",b,"_20220808.png",sep=""),height=1280,width=800)
  p1<-labeledHeatmap(Matrix = moduleTraitCor,
                     xLabels = names(datTraits),
                     yLabels = names(MEs),
                     ySymbols = names(MEs),
                     colorLabels = FALSE,
                     colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                     textMatrix = textMatrix,
                     setStdMargins = FALSE,
                     cex.text =1,
                     zlim = c(-1,1),
                     main = paste("Module-trait relationships"))
  
  dev.off()
  
  print("correlation plot")
  
  #=====================================================================================
  #
  #  Code chunk 9
  #
  #=====================================================================================
  
  
  geneinfo<-readRDS("data/geneinfo_r58395_20211210.rds")
  
  probes = colnames(datExpr)
  probes2annot = match(probes, geneinfo$Gene_ID)
  # The following is the number or probes without annotation:
  sum(is.na(probes2annot))
  
  # Create the starting data frame
  geneInfo = data.frame(gene = probes,
                        geneSymbol = geneinfo$Gene_Name[probes2annot],
                        entrez_id = geneinfo$entrezgene_id[probes2annot],
                        moduleColor = moduleColors)
  
  geneInfo$moduleID<-as.numeric(as.factor(geneInfo$moduleColor))
  
  saveRDS(geneInfo,paste("expr_mat/5a_WGCNA_relative/WGCNA_geneInfo_modules_N",b,"_20230411.rds",sep=""))
  
  
  ####gene_info
  
  ##############
  moduleTraitCor<-data.frame(moduleTraitCor)
  moduleTraitCor$N<-""
  moduleTraitCor$color<-gsub("^ME","",rownames(moduleTraitCor))  
  
  for ( i in 1:nrow(moduleTraitCor)){
    col<-gsub("^ME","",rownames(moduleTraitCor)[i])  
    moduleTraitCor$N[i]<-length(which(geneInfo$moduleColor==col))
  }
  
  saveRDS(moduleTraitCor,paste("expr_mat/5a_WGCNA_relative/moduleTraitCor_N",b,"_20230411.rds",sep=""))
  
  
  
  
  mod<-moduleTraitCor
  mod<-mod[!mod$color=="grey",]
  mod<-mod[order(as.numeric(as.character(mod$N)),decreasing=T),]
  
  modules<-unique(mod$color)
  
  
  
  ########PCA
  
  
  pcs_all<-c()
  pc1comp<-c()
  for ( i in 1:nrow(mod)){
    
    gg<-geneInfo$gene[geneInfo$moduleColor==mod$color[i]]
    
    mat<-exp[rownames(exp)%in% gg,]
    
    dat_z<-t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
    dat_z<-dat_z[apply(dat_z,1,function(x){length(which(is.na(x)))==0}),]
    
    sample<-meta$sample[match(colnames(dat_z),meta$library)]
    
    pca_prcomp <- prcomp(t(dat_z ), scale=F)
    pcs <- data.frame(predict(pca_prcomp))
    pcs<-pcs[,1:2]
    pcs$sample = sample
    pcs$module=mod$color[i]
    pcs$library=rownames(pcs)
    
    pcs_all<-rbind(pcs_all,pcs)
    
    pc1comp<-rbind(pc1comp,c(mod$color[i],summary(pca_prcomp)$importance[2,1]*100,summary(pca_prcomp)$importance[2,2]*100))
  }
  
  saveRDS(pcs_all,paste("expr_mat/5a_WGCNA_relative/PCA_value_N",b,"_20230411.rds",sep=""))
  
  
  pc1comp<-data.frame(pc1comp)
  colnames(pc1comp)<-c("modue","PC1comp","PC2comp")
  
  pc1comp$PC1comp<-as.numeric(as.character(pc1comp$PC1comp))
  pc1comp$PC2comp<-as.numeric(as.character(pc1comp$PC2comp))
  
  saveRDS(pc1comp,paste("expr_mat/5a_WGCNA_relative/PCA_comp_N",b,"_20230411.rds",sep=""))
  
  
#  "turquoise" "blue"  "brown"  "yellow"  "green"   "red" "black" "pink" 
 #     45        47      47       54         47      44     49      50
  
  head(pcs_all)
  
  pcs_all$PC12<-sqrt(pcs_all$PC1^2+pcs_all$PC2^2)
  pcs_all_1<-pcs_all[pcs_all$module %in% modules,]
  pcs_all_1$module<-factor(pcs_all_1$module,levels=modules,ordered=T)
  
  
 
  
  n<-ggplot(pcs_all_1,aes(x=PC1,y=1,fill=sample))+
    geom_point(aes(fill=sample,color=sample),size=8,alpha=0.5)+
    #geom_point(color="grey",shape=1,size=8)+
    scale_color_manual(values=colors.sample.Quartet.fill,name="Group")+
    theme_few()+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y =element_blank(),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text())+
    facet_wrap(module~.,scales="free",ncol=1);n
  
  
  ggsave(paste("chart/5a_WGCNA_relative/pca_pc1_N",b,"_20230411.png",sep=""),n,width=5,height=7)
  ggsave(paste("chart/5a_WGCNA_relative/pca_pc1_N",b,"_20230411.pdf",sep=""),n,width=5,height=7)
  
  print("PCA")
  
  
  ##################################################
  ##########functional analysis based on ensembl ID GO
  ##################################################
  
  umodules<-modules
  go_enrich<-c()
  for ( i in 1:length(umodules)){
    tab<-geneInfo[geneInfo$moduleColor==umodules[i],]
    g1<-as.character(unique(tab$gene))
    ########GOterm
    ego<-data.frame(enrichGO(g1, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = 'ALL', pvalueCutoff =  0.2, pAdjustMethod = "fdr", qvalueCutoff = 0.2))
    if(nrow(ego)>0){
      go_enrich<-rbind(go_enrich,cbind(umodules[i],ego))
    }
    print(i)
  }
  
  go_enrich$geneName<-""
  for ( i in 1:nrow(go_enrich)){
   # print(i)
    g<-unlist(strsplit(go_enrich$geneID[i],"/"))
    #gname<-geneInfo$geneSymbol[match(g,geneInfo$entrez_id)]
    gname<-geneInfo$geneSymbol[match(g,geneInfo$gene)]
    go_enrich$geneName[i]<-pastetogether(gname)
  }
  colnames(go_enrich)[1]<-"modulecolor"
  
  
  write.csv(go_enrich,paste("expr_mat/5a_WGCNA_relative/go_enrich_N",b,"_20230411.csv",sep=""))
  
  