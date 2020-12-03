load("basic_data.RData")
library(dplyr)
clear_sample <- clinical_info[(substring(clinical_info$sampleID,14,15) %in% c("01","11")) & (clinical_info$PAM50Call_RNAseq !=""),"sampleID"]
expr_t <- m6A_regulator_t[clear_sample,] %>% select(-sampleID)
sample_meta <- clinical_info[match(rownames(expr_t),clinical_info$sampleID),"PAM50Call_RNAseq"]
library(Rtsne)
tsne_expr <- Rtsne(apply(as.matrix(expr_t),2,as.numeric))
colors = c(Normal="#81b29a",Basal = "#C5194C",Her2 = "#fcd74e",LumA = "#6FBEC2",LumB = "#ED726B")
#colors = rainbow(length(unique(sample_meta)))
#names(colors) = unique(sample_meta)
pdf("tsneout_2.pdf",width = 6.8,height = 5.5)
plot(tsne_expr$Y,col=colors[sample_meta],asp=1,lwd = 1.5)
legend("topright",
       col = c("#81b29a","#C5194C","#fcd74e","#6FBEC2","#ED726B"),
       pch = 1,lwd=1.5,
       legend = c("Normal","Basal","Her2","LumA","LumB"))
dev.off()
