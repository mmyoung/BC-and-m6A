load("basic_data.RData")
##using the four basal-featured genes for samples clustering
library(ConsensusClusterPlus)
title= "sample_cluster"
Basal_sample_name <- clinical_info[clinical_info$PAM50Call_RNAseq=="Basal","sampleID"]
nonBasal_sample_name <- clinical_info[clinical_info$PAM50Call_RNAseq %in% c("LumA","LumB"),"sampleID"]
filter_sample_index <- colnames(m6A_expr) %in% c(Basal_sample_name,nonBasal_sample_name)
results = ConsensusClusterPlus(as.matrix(m6A_expr[c("IGF2BP2","IGF2BP3","YTHDC2","RBM15"),filter_sample_index]),
                               maxK=9,reps=50,
                               pItem=0.8,
                               pFeature=1,
                               clusterAlg="hc",
                               title = title,
                               distance="pearson",
                               seed=1262118,
                               plot="pdf")

####consensus clustering analysis with the DNA methylation of four basal-featured genes
selected_probes <- c("cg12781915","cg08939418","cg02860543",
                     "cg16466899","cg00508334","cg22826239",
                     "cg07297397","cg20265043","cg27135125",
                     "cg02302089","cg08584665"
)
title= "11probe_cluster"
filter_sample_index <- colnames(probe_filter_2) %in% c(Basal_sample_name,nonBasal_sample_name)
results = ConsensusClusterPlus(as.matrix(probe_filter_2[selected_probes,filter_sample_index]),
                               maxK=9,reps=50,pItem=0.8,pFeature=1,
                               title=title,
                               clusterAlg="hc",
                               distance="pearson",
                               seed=126211,plot="pdf")