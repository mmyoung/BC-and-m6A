ordered_gene <- c("METTL3","METTL14","WTAP","KIAA1429","RBM15","RBM15B","ZC3H13","CBLL1","METTL5","METTL16","ZCCHC4",
                  "FTO","ALKBH5",
                  "EIF3A","FMR1","HNRNPA2B1","HNRNPC","IGF2BP1","IGF2BP2","IGF2BP3","YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","PRRC2A","RBMX","ZNF217")
expr_raw_data <- read.table("expr_cox_regression.txt",header = T,row.names = 1,stringsAsFactors = F)
expr_raw_data <- expr_raw_data[,ordered_gene]
new_data <- as.data.frame(matrix(rep(0,28*3),nrow = 3))
colnames(new_data) <- colnames(raw_data)
row.names(new_data) <- c("OS","DFS","PFS")
for (i in colnames(expr_raw_data)){
  for (j in row.names(new_data)){
    if (expr_raw_data[paste(j,".HR",sep = ""),i] > 1 && expr_raw_data[paste(j,".P",sep = ""),i] < 0.05){
      new_data[j,i] <- 1
    }
    else if (expr_raw_data[paste(j,".HR",sep = ""),i] < 1 && expr_raw_data[paste(j,".P",sep = ""),i] < 0.05){
      new_data[j,i] <- (-1)
    }
  }
}
library(pheatmap)
colnames(new_data)[colnames(new_data)=="BAT2"] <- "PRRC2A"
colnames(new_data)[colnames(new_data)=="METT10D"] <- "METTL16"
gene_anno <- data.frame(group = c(rep("writer",11),rep("eraser",2),rep("reader",15)),
                        row.names = colnames(new_data))
ann_colors = list(
  group=c(writer="#ef00d1",eraser="yellow",reader="#ff7763")
)
pheatmap(new_data,
              scale = "none",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              cluster_cols= FALSE,
              cluster_rows = FALSE,
              cellwidth = 20,
              cellheight = 20,
              annotation_col = gene_anno,
              annotation_colors = ann_colors)


##copy number variation associated with survival heatmap
ordered_gene <- c("METTL3","METTL14","WTAP","KIAA1429","RBM15","RBM15B","ZC3H13","CBLL1","METTL5","METT10D","ZCCHC4",
                  "FTO","ALKBH5",
                  "EIF3A","FMR1","HNRNPA2B1","HNRNPC","IGF2BP1","IGF2BP2","IGF2BP3","YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","BAT2","RBMX","ZNF217")
raw_data <- read.table("clipboard",header = T,row.names = 1,stringsAsFactors = F)
raw_data <- raw_data[,ordered_gene]
raw_data[1:6,1:6]
gain_new_data <- as.data.frame(matrix(rep(0,28*3),nrow = 3))
colnames(gain_new_data) <- colnames(raw_data)
row.names(gain_new_data) <- c("OS.gain","DFS.gain","PFS.gain")
for (i in colnames(raw_data)){
  for (j in row.names(gain_new_data)){
    if (raw_data[paste(j,".HR",sep = ""),i] > 1 && raw_data[paste(j,".P",sep = ""),i] < 0.05){
      gain_new_data[j,i] <- 1
    }
    else if (raw_data[paste(j,".HR",sep = ""),i] < 1 && raw_data[paste(j,".P",sep = ""),i] < 0.05){
      gain_new_data[j,i] <- (-1)
    }
  }
}
raw_data <- read.table("CNV_cox_regression.txt",header = T,row.names = 1,stringsAsFactors = F)
raw_data <- raw_data[,ordered_gene]
raw_data[1:6,1:6]
loss_new_data <- as.data.frame(matrix(rep(0,28*3),nrow = 3))
colnames(loss_new_data) <- colnames(raw_data)
row.names(loss_new_data) <- c("OS.loss","DFS.loss","PFS.loss")
for (i in colnames(raw_data)){
  for (j in row.names(loss_new_data)){
    if (raw_data[paste(j,".HR",sep = ""),i] > 1 && raw_data[paste(j,".P",sep = ""),i] < 0.05){
      loss_new_data[j,i] <- 1
    }
    else if (raw_data[paste(j,".HR",sep = ""),i] < 1 && raw_data[paste(j,".P",sep = ""),i] < 0.05){
      loss_new_data[j,i] <- (-1)
    }
  }
}
library(ComplexHeatmap)
library(circlize)
UpColor <- colorRamp2(breaks = c(-1,0, 1), colors = c("navy", "white", "firebrick3"))
DnColor <- colorRamp2(breaks = c(-1, 0,1), colors = c("navy", "white", "firebrick3"))
col_an <- HeatmapAnnotation(type = c(rep("W", 11), rep("E", 2), rep("R", 15)),
                            show_annotation_name = F,
                            col = list(type = c( "W"="#ef00d1","E"="yellow","R"="#ff7763")),
                            show_legend = F
)
DiagFunc <- 
  function(up, down){
    function(j, i, x, y, width, height, fill){ 
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey")) 
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
    }
  }
p1 <- Heatmap(gain_new_data, 
              column_title = "", 
              rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              top_annotation = col_an, 
              cell_fun = DiagFunc(up = gain_new_data, down = loss_new_data),
              column_names_rot = 0,
              column_names_centered = TRUE,
              heatmap_width = unit(10, "npc"),
              heatmap_height = unit(1, "npc")
          
)
p1
pdf("complexplot_N2.pdf",width = 27,height = 3)
draw(p1, 
     annotation_legend_side = "bottom", 
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE)
dev.off()

##DNA methylation associated with survival
methy_raw_data <- read.table("methylation_cox_regression.txt",header = T,row.names = 1,stringsAsFactors = F)
library(dplyr)
methy_raw_data$probe <- row.names(methy_raw_data)
methy_raw_data_2 <- methy_raw_data[order(ordered(methy_raw_data$gene,levels=ordered_gene)),]
methy_raw_data_2$gene <- factor(methy_raw_data_2$gene,levels = unique(methy_raw_data_2$gene))
methy_raw_data_2 <- 
  do.call(rbind,lapply(split(methy_raw_data_2,methy_raw_data_2$gene),function(x) x[order(x[,"group"]),]))
head(methy_raw_data_2)
row.names(methy_raw_data_2) <- methy_raw_data_2$probe
head(methy_raw_data_2)
new_data <- as.data.frame(matrix(rep(0,74*3),nrow = 3))
colnames(new_data) <- row.names(methy_raw_data_2)
row.names(new_data) <- c("OS","DFS","PFS")
for (i in row.names(methy_raw_data_2)){
  for (j in row.names(new_data)){
    if (methy_raw_data_2[i,paste(j,".HR",sep = "")] > 1 && methy_raw_data_2[i,paste(j,".P",sep = "")] < 0.05){
      new_data[j,i] <- 1
    }
    else if (methy_raw_data_2[i,paste(j,".HR",sep = "")] < 1 && methy_raw_data_2[i,paste(j,".P",sep = "")] < 0.05){
      new_data[j,i] <- (-1)
    }
  }
}
library(pheatmap)
probe_anno <- data.frame(gene=unlist(methy_raw_data_2$gene),
                         group=unlist(methy_raw_data_2$group),
                         row.names = row.names(methy_raw_data_2))
anno_col <- list(gene = ,group=)

p <- pheatmap(new_data,
              scale = "none",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              cluster_cols= FALSE,
              cluster_rows = FALSE,
              cellwidth = 9,
              cellheight = 9,
              annotation_col = probe_anno)