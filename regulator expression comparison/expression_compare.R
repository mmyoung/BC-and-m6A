load("basic_data.RData")
wilcoxon.test_res <- read.table("TvsN_regulator_wilcox.txt",row.names=1,stringsAsFactors=F,header=T)
library(reshape2)
library(ggplot2)
library(dplyr)
merged_expr_meta <- merge(m6A_regulator_t,clinical_info[,c("sampleID","PAM50Call_RNAseq")],by="sampleID")
merged_expr_meta <- merged_expr_meta[merged_expr_meta$PAM50Call_RNAseq != "",]
melt_merged_expr_meta <- melt(merged_expr_meta)
melt_merged_expr_meta <- melt_merged_expr_meta %>% mutate(group= ifelse(PAM50Call_RNAseq=="Normal","Normal","Tumor"))
melt_merged_expr_meta$variable <- factor(melt_merged_expr_meta$variable,levels = m6A_regulator$human_m6A_regulators)
melt_merged_expr_meta$group <- factor(melt_merged_expr_meta$group,levels = c("Normal","Tumor"))
diff_res <- wilcoxon.test_res
##filtering the genes with difference significance of P<0.001 (17)
gene_3star <- row.names(diff_res[diff_res$TvsN_PValue<0.001,])
sig_data <- melt_merged_expr_meta[melt_merged_expr_meta$variable %in% gene_3star,]
gene_max_value <- summarise(group_by(sig_data,variable),max_value=max(value))
merged_data <- merge(gene_max_value,diff_res,by.x="variable",by.y=0)
star_sign <- as.character(symnum(merged_data$TvsN_PValue, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                 symbols = c("***", "**", "*", "ns"),
                                 abbr.colnames = FALSE, na = ""))
merged_data$star_sign <- star_sign
merged_data <- merged_data[match(gene_3star,merged_data$variable),]
merged_data$x_pos <- c(1:17)
merged_data$max_value <- merged_data$max_value+0.1
merged_data$color <- ifelse(merged_data$TvsN_diff<1,"#588b8b","#c8553d")
sig_data$variable <- as.character(sig_data$variable)
sig_data$variable[sig_data$variable=="BAT2"] <- "PRRC2A"
sig_data$variable[sig_data$variable=="METT10D"] <- "METTL16"
sig_data$variable <- factor(sig_data$variable,
                            levels = gsub("BAT2","PRRC2A",gsub("METT10D","METTL16",gene_3star)))
p1 <- ggplot(sig_data,aes(x=variable,y=value,fill=group))+
  geom_boxplot(width = 0.6,lwd=0.05)+
  scale_fill_manual(values=c("#81b29a","#e76f51"))+
  labs(x="",y="log2(RSEM+1)")+
  theme(axis.text.x = element_text(vjust = 0.7,hjust = 0.7,angle = 20,colour=merged_data$color),axis.title.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p2 <-  geom_text(data = merged_data,aes(x=merged_data$x_pos,y=merged_data$max_value,label=merged_data$star_sign),inherit.aes = FALSE)
p1+p2  
ggsave("gene_3star_expr_boxplot_TvsN.pdf",p1+p2,width = 10,height = 6)
##filtering the genes with difference significance of P<0.05 & >0.001 (11)
gene_non3star <- row.names(diff_res[diff_res$TvsN_PValue>0.001,])
sig_data <- melt_merged_expr_meta[melt_merged_expr_meta$variable %in% gene_non3star,]
gene_max_value <- summarise(group_by(sig_data,variable),max_value=max(value))
merged_data <- merge(gene_max_value,diff_res,by.x="variable",by.y=0)
star_sign <- as.character(symnum(merged_data$TvsN_PValue, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                 symbols = c("***", "**", "*", "ns"),
                                 abbr.colnames = FALSE, na = ""))
merged_data$star_sign <- star_sign
merged_data <- merged_data[match(gene_non3star,merged_data$variable),]
merged_data$x_pos <- c(1:11)
merged_data$max_value <- merged_data$max_value+0.1
merged_data$color <- ifelse(merged_data$star_sign=="ns","black",ifelse(merged_data$TvsN_diff>1,"#c8553d","#588b8b"))
sig_data$variable <- as.character(sig_data$variable)
sig_data$variable[sig_data$variable=="BAT2"] <- "PRRC2A"
sig_data$variable[sig_data$variable=="METT10D"] <- "METTL16"
sig_data$variable <- factor(sig_data$variable,
                            levels = gsub("BAT2","PRRC2A",gsub("METT10D","METTL16",gene_non3star)))
p1 <- ggplot(sig_data,aes(x=variable,y=value,fill=group))+
  geom_boxplot(width = 0.6,lwd=0.05)+
  scale_fill_manual(values=c("#81b29a","#e76f51"))+
  labs(x="",y="log2(RSEM+1)")+
  theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5,colour=merged_data$color),axis.title.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p2 <-  geom_text(data = merged_data,aes(x=merged_data$x_pos,y=merged_data$max_value,label=merged_data$star_sign),inherit.aes = FALSE)
p1+p2  
ggsave("gene_non3star_expr_boxplot_TvsN.pdf",p1+p2,width = 10,height = 6)

##heatmap of gene expression among four subtypes 
meta_variable_sel <- c("PAM50Call_RNAseq","sampleID",
                       "pathologic_T","pathologic_N",
                       "pathologic_M","pathologic_stage")
sample_info_filter <- clinical_info[clinical_info$sampleID %in% colnames(m6A_expr),][,meta_variable_sel] 
rownames(sample_info_filter) <- sample_info_filter$sampleID
sample_info_filter <- sample_info_filter[sample_info_filter$PAM50Call_RNAseq!="",]
ordered_gene <- c("METTL3","METTL14","WTAP","KIAA1429","RBM15","RBM15B","ZC3H13","CBLL1","METTL5","METT10D","ZCCHC4",
                  "ALKBH5","FTO",
                  "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","HNRNPA2B1","HNRNPC","RBMX","IGF2BP1","IGF2BP2","IGF2BP3","EIF3A","FMR1","BAT2","ZNF217")
ordered_regulator_expr <- m6A_expr[ordered_gene,]
library(tidyverse)
sample_info_filter_ordered <- sample_info_filter[order(sample_info_filter$PAM50Call_RNAseq),]
ordered_regulator_expr <- ordered_regulator_expr[,sample_info_filter_ordered$sampleID]
sample_anno <-  sample_info_filter_ordered %>%
  filter(PAM50Call_RNAseq!="Normal") %>%
  dplyr::select(PAM50Call_RNAseq,sampleID)%>%
  column_to_rownames(var = "sampleID")
colnames(sample_anno) <- "subgroup"
sample_anno$subgroup <- factor(sample_anno$subgroup,levels = c("Basal","Her2","LumA","LumB"))
gene_anno <- data.frame(group=c(rep("writer",11),rep("eraser",2),rep("reader",15)),
                        row.names = rownames(ordered_regulator_expr)
)
Pvalue <- list()
for (i in row.names(gene_anno)){
  y <- kruskal.test(m6A_regulator_t[row.names(sample_anno),i]~sample_anno$subgroup)
  Pvalue[i]<-y$p.value
}
star_sign <- as.character(symnum(as.numeric(Pvalue), cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                 symbols = c("***", "**", "*", "ns"),
                                 abbr.colnames = FALSE, na = ""))
gene_Pvalue <- data.frame(gene=names(Pvalue),
                          Pvalue=as.numeric(Pvalue),
                          star_sign=star_sign)
ann_colors = list(
  subgroup=c(Basal = "#C5194C",Her2 = "#fcd74e",LumA = "#6FBEC2",LumB = "#ED726B"),
  group=c(writer="#ef00d1",eraser="yellow",reader="#ff7763")
)
library(pheatmap)
heatmap_data <- ordered_regulator_expr
row.names(heatmap_data) <- paste(gene_Pvalue$gene,gene_Pvalue$star_sign)
row.names(heatmap_data) <- gsub("METT10D","METTL16",gsub("BAT2","PRRC2A",row.names(heatmap_data)))
row.names(gene_anno) <- paste(gene_Pvalue$gene,gene_Pvalue$star_sign)
row.names(gene_anno) <- gsub("METT10D","METTL16",gsub("BAT2","PRRC2A",row.names(gene_anno)))
ppp<-pheatmap(heatmap_data[,row.names(sample_anno)],
              cluster_cols = F,
              scale = "row",
              cluster_rows = F,
              #cutree_rows = 5,
              show_colnames =F,
              show_rownames = T,
              annotation_col = sample_anno ,
              annotation_row = gene_anno,
              annotation_colors = ann_colors,
              gaps_row = c(11,13),
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave(filename = "regulator_heatmap_4subtype.pdf",ppp,width = 10,height = 6)
