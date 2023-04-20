# process data
library(reshape2)
library(xlsx)
library(data.table)
library(stringr)

#draw plot
library(RColorBrewer)
library(matrixStats)
library(Vennerable)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(gg.gap)
library(ggthemes)
library(GGally)
library(ggalt)
#library(ggpmisc)
library(ggsci)
library(pheatmap)
library(scales)

## functional analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
####scatterplot correlation
library(smplot)
######density point
library(ggpointdensity)

### RNA-seq analysis
library(edgeR)
library(limma)


#20230405
## check version package.version ()

#R.Version() ## R version 4.1.2 (2021-11-01)

#> package.version("ggplot2")
#[1] "3.3.5"
#> package.version("pheatmap")
#[1] "1.0.12"
#> package.version("clusterProfiler")
#[1] "4.2.2"
#> package.version("edgeR")
#[1] "3.36.0"
#> package.version("limma")
#[1] "3.50.0"
#> package.version("WGCNA")
#[1] "1.71"
#> package.version("ggsci")
#[1] "2.9"
#> package.version("GGally")
#[1] "2.1.2"
