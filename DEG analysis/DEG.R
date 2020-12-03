load(file="consensuscluster.RData")
load(file = "basic_data.RData")
cluter_1_sample <- names(results[[2]]$consensusClass[results[[2]]$consensusClass==1])
cluter_2_sample <- names(results[[2]]$consensusClass[results[[2]]$consensusClass==2])
cluster1_LumA <- intersect(cluter_1_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="LumA","sampleID"])
cluster1_LumB <- intersect(cluter_1_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="LumB","sampleID"])
cluster1_Her2 <- intersect(cluter_1_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="Her2","sampleID"])
cluster2_LumA <- intersect(cluter_2_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="LumA","sampleID"])
cluster2_LumB <- intersect(cluter_2_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="LumB","sampleID"])
cluster2_Her2 <- intersect(cluter_2_sample,clinical_info[clinical_info$PAM50Call_RNAseq=="Her2","sampleID"])
library(DESeq2)
row.names(HiSeqV2) <- HiSeqV2$sample
diff_expr_df <- HiSeqV2[,c(cluter_1_sample,cluter_2_sample)]
expr_df_rsem <- 2^diff_expr_df-1
expr_df_round <-round(expr_df_rsem,digits = 0)
expr_df_round[1:6,1:3]
res_list <- list()
for (i in c("LumA","LumB")){
  cluster1_sample_name <- get(paste("cluster1_",i,sep = ""))
  cluster2_sample_name <- get(paste("cluster2_",i,sep = ""))
  coldata_tmp <- data.frame(row.names=c(cluster1_sample_name,cluster2_sample_name), 
                            condition=c(rep("cluster1",length(cluster1_sample_name)),rep("cluster2",length(cluster2_sample_name))))
  coldata_tmp$condition <- factor(coldata_tmp$condition,levels = c("cluster1","cluster2"))
  dds_tmp <- DESeqDataSetFromMatrix(countData=expr_df_round[,c(cluster1_sample_name,cluster2_sample_name)],
                                    colData=coldata_tmp, design=~condition)
  dds_tmp <- DESeq(dds_tmp)
  res_tmp <- results(dds_tmp,contrast=c("condition","cluster1","cluster2"))
  diff_gene_df <- as.data.frame(res_tmp)[which(abs(res_tmp$log2FoldChange)>1 & res_tmp$padj<0.05),]
  res_list[[i]] <- res_tmp
  write.table(diff_gene_df,
              file = paste("cluster_DESeq2\\",i,"_DESeq2_sig_diff_gene.txt",sep=""),
              sep = "\t",quote = F,row.names = T)
  write.table(res_tmp,
              file = paste("cluster_DESeq2\\",i,"_DESeq2_res.txt",sep=""),
              sep = "\t",quote = F,row.names = T)
}