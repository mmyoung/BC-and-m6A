load("methylation_basic_data.RData")
wilcox_res <- read.table("tumorVSnormal_wilcox.txt",header=T,row.names=1,stringsAsFactors=F)
sig_diff_probes <- wilcox_res[abs(wilcox_res$tumorVSnormal_diff)>0.2 & wilcox_res$tumorVSnormal_PValue<0.05,]
dim(sig_diff_probes)
library(pheatmap)
probe_anno <- data.frame(gene=sig_diff_probes$gene,
                         group=sig_diff_probes$Group,
                         row.names = row.names(sig_diff_probes),
                         stringsAsFactors = F)
probe_anno$gene[grep(";",probe_anno$gene)] <- 
  unlist(lapply(probe_anno$gene[grep(";",probe_anno$gene)],function(x){str_split(x,";")[[1]][1]}))
probe_anno$group[grep(";",probe_anno$group)] <- 
  unlist(lapply(probe_anno$group[grep(";",probe_anno$group)],function(x){str_split(x,";")[[1]][1]}))
probe_anno$probes <- row.names(probe_anno)
probe_anno <- probe_anno[order(ordered(probe_anno$gene,levels = m6A_regulator$human_m6A_regulators)),]
row.names(probe_anno) <- probe_anno$probes 
probe_anno <- probe_anno %>% select(gene,group)
probe_anno$gene <- factor(probe_anno$gene,levels = unique(probe_anno$gene))
probe_anno$group <- factor(probe_anno$group,levels = unique(probe_anno$group))
sample_anno <- data.frame(status=c(rep("normal",length(normal_sample)),rep("tumor",length(tumor_sample))),
                          row.names=c(normal_sample,tumor_sample))
gene <- brewer.pal(6,"Set3")
names(gene) <- levels(probe_anno$gene)
group <- brewer.pal(5,"Dark2")
names(group) <- levels(probe_anno$group)
ann_colors <- list(
  status=c(normal="#81b29a",tumor="#e76f51"),
  group=group,
  gene=gene
)
pheatmap(as.matrix(probe_filter_2[row.names(probe_anno),c(normal_sample,tumor_sample)]),
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = F,
         show_rownames = T,
         scale = "row",
         annotation_row = probe_anno,
         annotation_col = sample_anno,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

##hierarchical clustering of DNA methylation probes whose SD is greater than 0.2 in all tumor samples
minus_normal_sample <- colnames(probe_filter_2)[!colnames(probe_filter_2) %in% clinical_info[clinical_info$PAM50Call_RNAseq == "Normal","sampleID"]]
probe_sd <- apply(probe_filter_2[,minus_normal_sample],1,sd)
probe_sd_sort <- sort(probe_sd,decreasing = T)
probe_sd_filter <- names(probe_sd_sort[probe_sd_sort>0.2]) 
length(probe_sd_filter)
sample_anno <- data.frame(subtype= clinical_info[match(minus_normal_sample,clinical_info$sampleID),"PAM50Call_RNAseq"],
                          sample = minus_normal_sample)
sample_anno <- sample_anno[order(sample_anno$subtype),]
sample_annot <- data.frame(subtype=sample_anno$subtype,
                           row.names = sample_anno$sample)
sample_annot$subtype <- factor(sample_annot$subtype)
probe_anno <- data.frame(gene=reflect_df[match(probe_sd_filter,row.names(reflect_df)),"gene"],
                         row.names = probe_sd_filter)
probe_anno$gene <- factor(probe_anno$gene,levels = unique(probe_anno$gene))
gene <- brewer.pal(8,"Set3")
names(gene) <- levels(probe_anno$gene)
ann_colors = list(
  subtype=c(Basal="#C5194C",Her2="#fcd74e",LumA="#6FBEC2",LumB="#ED726B"),
  gene=gene
)
pheatmap(as.matrix(probe_filter_2[probe_sd_filter,row.names(sample_annot)]),
         cluster_rows = T,
         cluster_cols = T,
         scale = "row",
         show_rownames = T,
         show_colnames = F,
         clustering_method = "complete",
         annotation_row = probe_anno,
         annotation_col = sample_annot,
         annotation_colors = ann_colors
)