summary_res <- read.table("CNV_tumor_normal_freq.txt",
            row.names = F,sep = "\t",quote = F)
TN_plot_data <- summary_res
TN_plot_data$del_percent <- -1*TN_plot_data$del_percent
TN_plot_data$gene <- factor(TN_plot_data$gene,levels = m6A_regulator$human_m6A_regulators)

ggplot(data=TN_plot_data,aes(x=gene))+
  geom_bar(aes(x=gene,y=del_percent,fill=subtype),stat = "identity",position = "dodge")+
  geom_bar(aes(x=gene,y=amp_percent,fill=subtype),stat = "identity",position = "dodge")+
  theme_classic()+
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle = 30,margin = unit(c(0,0.5,0.5,0.5),"cm"),hjust=1, vjust=1))+
  scale_y_continuous(breaks = c(-50,0,50),labels = c("50%","0","50%"))+
  labs(x="",y="frequency",fill="")+
  scale_fill_manual(values = c("#81b29a","#e76f51"))+
  annotate("text",x=23,y=55,label="gain",size=6)+
  annotate("text",x=23,y=-55,label="loss",size=6)

##CNV frequency complexheatmap among four subtypes
summary_res <- read.table(file = "CNV_subtype_normal_freq.txt",header = T,
              sep = "\t",quote = "")
summary_res <- summary_res[summary_res$subtype!="Normal",]  
library(ComplexHeatmap)
library(circlize)
library(tidyr)
options(stringsAsFactors = FALSE)
complexheatmap_data_pre <- summary_res %>% 
  mutate(del_frequency=del_sample/total_sample,amp_frequency=amp_sample/total_sample)
head(complexheatmap_data_pre)
complexheatmap_data_pre$amp_frequency <- round(complexheatmap_data_pre$amp_frequency,digits = 2)
complexheatmap_data_pre$del_frequency <- round(complexheatmap_data_pre$del_frequency,digits = 2)
up <- complexheatmap_data_pre[,c("gene","subtype","amp_frequency")] %>%
  spread(subtype,amp_frequency) %>% 
  column_to_rownames(var = "gene")
up <- up[m6A_regulator$human_m6A_regulators,]
dn <- complexheatmap_data_pre[,c("gene","subtype","del_frequency")] %>%
  spread(subtype,del_frequency) %>%
  column_to_rownames(var = "gene")
dn <- dn[m6A_regulator$human_m6A_regulators,]

UpColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFADD","#AB221F"))
DnColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFADD","#3878C1"))

row_an <- HeatmapAnnotation(type = c(rep("W", 11), rep("E", 2), rep("R", 15)),
                            show_annotation_name = F,
                            col = list(type = c( "W" = "#FAC67A", "E" = "#51B743","R" = "#5AC9FA")),
                            show_legend = T,
                            annotation_legend_param = list(title = "m6A group",nrow = 1),
                            which = "row"
                            )
DiagFunc <- 
  function(up, down){
  function(j, i, x, y, width, height, fill){ 
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                 unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey")) 
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                 unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                 gp = gpar(fill = DnColor(dn[i, j]), col = "grey")) 
    grid.text(up[i, j], x-0.35*width, y+0.15*height)
    grid.text(dn[i, j], x+0.35*width, y-0.15*height)
  }
  }

p1 <- Heatmap(up, 
              column_title = "", 
              rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              left_annotation = row_an, 
              cell_fun = DiagFunc(up = up, down = dn),
              column_names_rot = 0,
              column_names_centered = TRUE
              )
p1

lgd <- list(Legend(title = "CNV Gain", 
                   col_fun = UpColor, 
                   at = c(0,0.5,1),  
                   direction = "horizontal"
                   ), 
            Legend(title = "CNV Loss", 
                   col_fun = DnColor,
                   at = c(0,0.5,1), 
                   direction = "horizontal"
                   ) 
            )
pdf("complexplot_N2.pdf",width = 7,height = 10)
draw(p1, annotation_legend_list = lgd, 
     annotation_legend_side = "bottom", 
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE)
dev.off()