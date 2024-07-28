##############################################################
########################################  ####################
#######################################     ##################
#######################################      #################
###############                DMP analysis         ##########
#######################################      #################
#######################################     ##################
########################################  ####################
##############################################################
library(rtracklayer)
library(dplyr)
library(do)
library(limma)
library(ggplot2)
library(ggpubr)
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step3-output-myDMP_hbv_T_N.Rdata")
myDMP_hbv <- myDMP

load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step3-output-myDMP_hcv_T_N.Rdata")
myDMP_hcv <- myDMP
rm(myDMP)

myDMP_hbv_hyper <- myDMP_hbv$T_to_N[myDMP_hbv$T_to_N$logFC >=0,]
myDMP_hbv_hypo <- myDMP_hbv$T_to_N[myDMP_hbv$T_to_N$logFC < 0,]
myDMP_hcv_hyper <- myDMP_hcv$T_to_N[myDMP_hcv$T_to_N$logFC >=0,]
myDMP_hcv_hypo <- myDMP_hcv$T_to_N[myDMP_hcv$T_to_N$logFC < 0,]

dat_for_chi_hbv <- myDMP_hbv_hypo
dat_for_chi_hcv <- myDMP_hcv_hypo
###chi test for feature
if(T){
  table(dat_for_chi_hbv$feature)
  table(dat_for_chi_hcv$feature)
  chisq.test(c(dat_for_chi_hbv$feature,dat_for_chi_hcv$feature),
             c(rep("HBV",length(dat_for_chi_hbv$feature)),
               rep("HCV",length(dat_for_chi_hcv$feature))))
}





################################################################################
#############--------- hyper/hypo methylated percentage-----------##############
################################################################################
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation_table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_table <- as.data.frame(annotation_table)
annotation_table <- annotation_table[c("chr","pos","Name","Relation_to_Island","UCSC_RefGene_Group")]
annotation_table$UCSC_RefGene_Group[which(annotation_table$UCSC_RefGene_Group == "")] <- "IGR"
annotation_table$feature <- strsplit2(annotation_table$UCSC_RefGene_Group,";")[,1]


dmp_hbv <- cbind("Prob"=rownames(myDMP_hbv$T_to_N),
                 #"group" = rep("HBV", nrow(myDMP_hbv$T_to_N)),
                 "Meth"= ifelse(myDMP_hbv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))
dmp_hcv <- cbind("Prob"=rownames(myDMP_hcv$T_to_N),
                 #"group" = rep("HCV", nrow(myDMP_hcv$T_to_N)),
                 "Meth"= ifelse(myDMP_hcv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))



dat_for_bar_hyper <- data.frame("region_name" =rep(names(table(annotation_table$feature)),2),
                                "Totol_prob"= rep(as.vector(table(annotation_table$feature)),2),
                                "Hyper_prob"=c(as.vector(table(myDMP_hbv_hyper$feature)),as.vector(table(myDMP_hcv_hyper$feature))),
                                "group"=c( rep("HBV",7), rep("HCV",7))
                                )
dat_for_bar_hyper$percent <- round((dat_for_bar_hyper$Hyper_prob/dat_for_bar_hyper$Totol_prob)*100,2)

dat_for_bar_hypo <-  data.frame("region_name" =rep(names(table(annotation_table$feature)),2),
                                "Totol_prob"= rep(as.vector(table(annotation_table$feature)),2),
                                "Hypo_prob"=c(as.vector(table(myDMP_hbv_hypo$feature)),as.vector(table(myDMP_hcv_hypo$feature))),
                                "group"=c( rep("HBV",7), rep("HCV",7))
                                )
dat_for_bar_hypo$percent <- round((dat_for_bar_hypo$Hypo_prob/dat_for_bar_hypo$Totol_prob)*100,2)

p_5utr <- prop.test(dat_for_bar_hyper$Hyper_prob[c(3,10)], dat_for_bar_hyper$Totol_prob[c(3,10)],alternative = "two.sided")
p_igr <- prop.test(dat_for_bar_hyper$Hyper_prob[c(5,12)], dat_for_bar_hyper$Totol_prob[c(5,12)],alternative = "two.sided")
p_tss1500 <- prop.test(dat_for_bar_hyper$Hyper_prob[c(6,13)], dat_for_bar_hyper$Totol_prob[c(6,13)],alternative = "two.sided")
if(T){
  pdf("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/bar_plot_hyper.pdf")
  ggplot(dat_for_bar_hyper,aes(x= region_name, percent)) + 
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.8,aes(fill = group)) + 
    theme_bw() +
    scale_x_discrete(limits = c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","IGR")) +
    scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.5, colour = "black")) +
    labs(x = "Region", y = "Proportion (%)") +
    annotate("text", label ="Adjusted p < 0.0001", x= Inf, y= Inf, vjust = 1, hjust = 1)
  dev.off()
}
################################################################################
#################------------------bar plot--------------#######################
################################################################################
if(T){
  pdf("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/bar_plot_hypo.pdf")
  ggplot(dat_for_bar_hypo,aes(x= region_name, percent)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.8,aes(fill = group)) + 
    theme_bw() +
    scale_x_discrete(limits = c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","IGR")) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.5, colour = "black")) +
    labs(x = "Region", y = "Proportion (%)") +
    annotate("text", label ="Adjusted p < 0.0001", x= Inf, y= Inf, vjust = 1, hjust = 2.5)
  dev.off()
}
################################################################################
#################-----------------pie chart--------------#######################
################################################################################
if(T){
  numbers <- as.vector(table(myDMP_hcv_hypo$feature))
  name <- names(table(myDMP_hcv_hypo$feature))
  percent <- 100*(numbers/sum(numbers)) %>% round(.,4) 
  percent <- paste(percent,"%",sep = "")
  pdf(file = "~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/hcv/pie_t_n_feat_hypo.pdf", width = 8,height = 8)
  pie(numbers, labels = percent, col = rainbow(length(percent)))
  #legend("topleft", name, fill = rainbow(length(percent)), cex = 0.8)
  dev.off()
}

################################################################################
################################################################################
#################-----------------chi test for cgi-------#######################
################################################################################
################################################################################

if(T){
  table(dat_for_chi_hbv$cgi)
  table(dat_for_chi_hcv$cgi)
  chisq.test(c(dat_for_chi_hbv$cgi,dat_for_chi_hcv$cgi),
             c(rep("HBV",length(dat_for_chi_hbv$cgi)),
               rep("HCV",length(dat_for_chi_hcv$cgi))))
}



################################################################################
################################################################################
#################-----------------DMP comparison---------#######################
################################################################################
################################################################################
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step3-output-myDMP_hbv_T_N.Rdata")
myDMP_hbv <- myDMP

load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step3-output-myDMP_hcv_T_N.Rdata")
myDMP_hcv <- myDMP
rm(myDMP)
dmp_hbv <- cbind("Prob"=rownames(myDMP_hbv$T_to_N),
                 #"group" = rep("HBV", nrow(myDMP_hbv$T_to_N)),
                 "Meth"= ifelse(myDMP_hbv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))
dmp_hcv <- cbind("Prob"=rownames(myDMP_hcv$T_to_N),
                 #"group" = rep("HCV", nrow(myDMP_hcv$T_to_N)),
                 "Meth"= ifelse(myDMP_hcv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))

dmp_hbv <- as.data.frame(dmp_hbv)
dmp_hcv <- as.data.frame(dmp_hcv)
rownames(dmp_hbv) <- dmp_hbv$Prob
rownames(dmp_hcv) <- dmp_hcv$Prob

library(VennDiagram)
venn <- venn.diagram(list(HBV = c(rownames(dmp_hbv)), 
                          HCV = c(rownames(dmp_hcv))), 
                     filename = NULL,
                     fill = c("red", "green"), 
                     alpha = c(0.5, 0.5), 
                     cex = 1.5,cat.fontface = 4
                     #lty =2
                     )
grid.draw(venn)
pdf("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/venn.pdf")
grid.draw(venn)
dev.off()

#unique dmp analysis
dmp_hbv <- cbind("Prob"=rownames(myDMP_hbv$T_to_N),
                 #"group" = rep("HBV", nrow(myDMP_hbv$T_to_N)),
                 "Meth"= ifelse(myDMP_hbv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))
dmp_hcv <- cbind("Prob"=rownames(myDMP_hcv$T_to_N),
                 #"group" = rep("HCV", nrow(myDMP_hcv$T_to_N)),
                 "Meth"= ifelse(myDMP_hcv$T_to_N$logFC >= 0,"Hypermethylation","Hypomethylation"))
dmp_hbv <- as.data.frame(dmp_hbv)
dmp_hcv <- as.data.frame(dmp_hcv)
rownames(dmp_hbv) <- dmp_hbv$Prob
rownames(dmp_hcv) <- dmp_hcv$Prob
unique_dmp_hbv <- dmp_hbv[!rownames(dmp_hbv) %in% intersect(dmp_hbv$Prob,dmp_hcv$Prob),]
unique_dmp_hbv$virus <- "HBV"
unique_dmp_hcv <- dmp_hcv[!rownames(dmp_hcv) %in% intersect(dmp_hbv$Prob,dmp_hcv$Prob),]
unique_dmp_hcv$virus <- "HCV"
unique_dmp <- rbind(unique_dmp_hbv, unique_dmp_hcv)
write.csv(unique_dmp, "~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/unique_dmp_no_psuedo.csv")
dmp_hbv_share <- dmp_hbv[intersect(dmp_hbv$Prob,dmp_hcv$Prob),]
dmp_hcv_share <- dmp_hcv[intersect(dmp_hbv$Prob,dmp_hcv$Prob),]
dmp_hbv_share <- rename(dmp_hbv_share,Meth = "Meth_hbv")
dmp_hcv_share <- rename(dmp_hcv_share,Meth = "Meth_hcv")
dmp_share <- merge(dmp_hbv_share,dmp_hcv_share, by = "Prob")
true_dmp_share <- rbind(dmp_share[dmp_share$Meth_hbv == "Hypermethylation" & dmp_share$Meth_hcv == "Hypermethylation",],
                        dmp_share[dmp_share$Meth_hbv == "Hypomethylation" & dmp_share$Meth_hcv == "Hypomethylation",])

psuedo_dmp_share <- rbind(dmp_share[dmp_share$Meth_hbv == "Hypermethylation" & dmp_share$Meth_hcv == "Hypomethylation",],
                          dmp_share[dmp_share$Meth_hbv == "Hypomethylation" & dmp_share$Meth_hcv == "Hypermethylation",])
write.csv(psuedo_dmp_share, "~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/psuedo_dmp_share.csv")

#heatmap for psuedo dmp in common
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step2-output-myNorm_hbv.Rdata")
myDMP_hbv_heatmap <- myNorm
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step2-output-myNorm_hcv.Rdata")
myDMP_hcv_heatmap <- myNorm
rm(myNorm)

myDMP_hbv_heatmap <- as.data.frame(myDMP_hbv_heatmap)
colnames(myDMP_hbv_heatmap) <- paste("HBV",colnames(myDMP_hbv_heatmap), sep  = "-")

myDMP_hcv_heatmap <- as.data.frame(myDMP_hcv_heatmap)
colnames(myDMP_hcv_heatmap) <- paste("HCV",colnames(myDMP_hcv_heatmap), sep  = "-")
ph_dat <- cbind(myDMP_hbv_heatmap, myDMP_hcv_heatmap)
ph_dat <- ph_dat[psuedo_dmp_share$Prob,]
ph_dat <- cbind(ph_dat[substr(colnames(ph_dat),1,5) == "HBV-T"],ph_dat[substr(colnames(ph_dat),1,5) == "HBV-N"], ph_dat[substr(colnames(ph_dat),1,5) == "HCV-N"], ph_dat[substr(colnames(ph_dat),1,5) == "HCV-T"])
if(T){
  data_directory <- setwd("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/")
  rt = ph_dat
  Type =c(rep("HBV-HCC", table(substr(colnames(ph_dat),1,5))["HBV-T"]),
          rep("Non-Tumor", table(substr(colnames(ph_dat),5,5))["N"]),
          rep("HCV-HCC", table(substr(colnames(ph_dat),1,5))["HCV-T"])
  )
  #Type <- substr(colnames(rt),1,3)
  names(Type) = colnames(rt)
  Type = as.data.frame(Type)
  heat_plot_pdf <- paste0(data_directory,"/Heatmap_psuedo_dmp_share.pdf")
  pdf(heat_plot_pdf, height = 8, width = 13)
  pheatmap::pheatmap(rt,
                     annotation = Type,
                     color = colorRampPalette(c("blue","white","red"))(50),
                     scale = 'row',
                     cluster_cols = F,
                     fontsize = 9,
                     fontsize_row = 8,
                     fontsize_col = 6,
                     main = "",
                     show_colnames = F,
                     show_rownames = F)
  dev.off()
}

#pie chart for uni cpg sites
if(T){
  pie_dat <- myDMP_hbv$T_to_N[unique_dmp_hbv$Prob,]
  pie_dat <- pie_dat[pie_dat$logFC > 0,]
  numbers <- as.vector(table(pie_dat$feature))
  name <- names(table(pie_dat$feature))
  percent <- 100*(numbers/sum(numbers)) %>% round(.,4) 
  percent <- paste(percent,"%",sep = "")
  pdf(file = "~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/hbv/pie_t_n_feat_uni_hyper.pdf", width = 8,height = 8)
  pie(numbers, labels = percent, col = rainbow(length(percent)))
  legend("topleft", name, fill = rainbow(length(percent)), cex = 0.8)
  dev.off()
  rm(pie_dat)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


if(T){
  protein_coding_anno <-  subset(gene_anno,gene_type == "protein_coding")
  protein_coding_anno <- protein_coding_anno[c(1,2,3,4,5,7,10,12)]
  protein_coding_anno <- subset(protein_coding_anno,type == "gene")
  rownames(protein_coding_anno) <- protein_coding_anno$gene_name
}
meth <- myDMP$T_to_N[c("adj.P.Val","B","logFC","CHR","MAPINFO","Strand","Type","gene","feature","cgi")]
meth$probeId <- rownames(meth)

if(T){
  dist <- as.data.frame(cbind((protein_coding_anno["TMEM82","start"] - meth[meth$CHR=="1",]$MAPINFO)/1000,
                              (protein_coding_anno["TMEM82","end"] - meth[meth$CHR=="1",]$MAPINFO)/1000))
  dist$probeId <- meth[meth$CHR=="1",]$probeId
  genebody <- as.numeric(dist$V1*dist$V2 < 0)
  ifelse(dist$V1-dist$V2 < 0,relative_dist <- dist$V1,relative_dist <- dist$V2)
  
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(createKEGGdb)
setwd("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/annotation_dat")
species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
#createKEGGdb::create_kegg_db(species)
# You will get KEGG.db_1.0.tar.gz file in your working directory
library(KEGG.db)
data(gcSample)

pie_dat <- myDMP_hbv$T_to_N[psuedo_dmp_share$Prob,]
head(pie_dat$gene)
if(T){
  #overall kegg and go
  if(T){
  entrez <- bitr(unique(na.omit(pie_dat$gene)), fromType = "SYMBOL",
                 toType = c( "ENTREZID"),
                 OrgDb = org.Hs.eg.db)
  kegg <- enrichKEGG(entrez$ENTREZID,
                     organism = "hsa",
                     #keyType = "kegg",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     #universe,
                     minGSSize = 100,
                     maxGSSize = 5000,
                     qvalueCutoff = 0.05,
                     use_internal_data = T
  )
  kegg@result$p.adjust <- signif(kegg@result$p.adjust,2)
  GO_BP <-enrichGO(entrez$ENTREZID,
                   ont = "BP",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  GO_CC <-enrichGO(entrez$ENTREZID,
                   ont = "CC",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  GO_MF <-enrichGO(entrez$ENTREZID,
                   ont = "MF",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  GO_BP@result$p.adjust <- signif(GO_BP@result$p.adjust,2)
  GO_CC@result$p.adjust <- signif(GO_CC@result$p.adjust,2)
  GO_MF@result$p.adjust <- signif(GO_MF@result$p.adjust,2)
  kegg <- as.data.frame(kegg)
  GO_BP <- as.data.frame(GO_BP)
  GO_CC <- as.data.frame(GO_CC)
  GO_MF <- as.data.frame(GO_MF)
  kegg$logp <- -log10(kegg$p.adjust)
  kegg <- kegg[order(kegg$logp, decreasing = T),]
  GO_BP$logp <- -log10(GO_BP$p.adjust)
  GO_BP <- GO_BP[order(GO_BP$logp, decreasing = T),]
  GO_CC$logp <- -log10(GO_CC$p.adjust)
  GO_CC <- GO_CC[order(GO_CC$logp, decreasing = T),]
  GO_MF$logp <- -log10(GO_MF$p.adjust)
  GO_MF <- GO_MF[order(GO_MF$logp, decreasing = T),]
  GO <- rbind(GO_BP[1:10,],GO_CC[1:10,],GO_MF[1:10,])
  GO$group <- c(rep("BP",10), rep("CC",10), rep("MF",10))
  p_kegg <-ggplot(kegg[1:10,], aes(x = logp, y = Description, fill = logp)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(kegg[1:10,]$Description)) +
    theme(panel.background = element_blank(),
          axis.text.y.left = element_text(size = 14)) +
    labs(x = "-log10(adjustedPvalue)", fill = "-log(adjustedP)", y = "")
  p_bp <- ggplot(GO_BP[1:10,], aes(x = logp, y = Description, fill = logp)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(GO_BP[1:10,]$Description)) +
    theme(panel.background = element_blank(),
          axis.text.y.left = element_text(size = 12)) +
    labs(x = "-log10(adjustedPvalue)", fill = "", y = "")
  p_cc<- ggplot(GO_CC[1:10,], aes(x = logp, y = Description, fill = logp)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(GO_CC[1:10,]$Description)) +
    theme(panel.background = element_blank(),
          axis.text.y.left = element_text(size = 12)) +
    labs(x = "-log10(adjustedPvalue)", fill = "", y = "")
  p_mf <- ggplot(GO_MF[1:10,], aes(x = logp, y = Description, fill = logp)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(GO_MF[1:10,]$Description)) +
    theme(panel.background = element_blank(),
          axis.text.y.left = element_text(size = 12)) +
    labs(x = "-log10(adjustedPvalue)", fill = "", y = "") 
  p_go <- ggplot(GO, aes(x = logp, y = Description, fill = group)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(GO$Description)) +
    theme(panel.background = element_blank(),
          axis.text.y.left = element_text(size = 12)) +
    labs(x = "-log10(adjustedPvalue)", fill = "", y = "") 
  }
  ################-------------------------##################
  #hyper/hypo kegg and go
  if(T){
    pie_dat_hyper <- pie_dat[pie_dat$logFC > 0,]
    pie_dat_hypo <- pie_dat[pie_dat$logFC < 0,]
    ####----####
    #hyper
    if(T){
      entrez <- bitr(unique(na.omit(pie_dat_hyper$gene)), fromType = "SYMBOL",
                     toType = c( "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
      kegg <- enrichKEGG(entrez$ENTREZID,
                         organism = "hsa",
                         #keyType = "kegg",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         #universe,
                         minGSSize = 100,
                         maxGSSize = 5000,
                         qvalueCutoff = 0.05,
                         use_internal_data = T
      )
      kegg@result$p.adjust <- signif(kegg@result$p.adjust,2)
      GO_BP <-enrichGO(entrez$ENTREZID,
                       ont = "BP",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_CC <-enrichGO(entrez$ENTREZID,
                       ont = "CC",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_MF <-enrichGO(entrez$ENTREZID,
                       ont = "MF",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_BP@result$p.adjust <- signif(GO_BP@result$p.adjust,2)
      GO_CC@result$p.adjust <- signif(GO_CC@result$p.adjust,2)
      GO_MF@result$p.adjust <- signif(GO_MF@result$p.adjust,2)
      kegg <- as.data.frame(kegg)
      GO_BP <- as.data.frame(GO_BP)
      GO_CC <- as.data.frame(GO_CC)
      GO_MF <- as.data.frame(GO_MF)
      kegg$logp <- -log10(kegg$p.adjust)
      kegg_hyper <- kegg[order(kegg$logp, decreasing = T),]
      kegg_hyper <- kegg_hyper[1:10,]
      GO_BP$logp <- -log10(GO_BP$p.adjust)
      GO_BP <- GO_BP[order(GO_BP$logp, decreasing = T),]
      GO_CC$logp <- -log10(GO_CC$p.adjust)
      GO_CC <- GO_CC[order(GO_CC$logp, decreasing = T),]
      GO_MF$logp <- -log10(GO_MF$p.adjust)
      GO_MF <- GO_MF[order(GO_MF$logp, decreasing = T),]
      GO_hyper <- rbind(GO_BP[1:10,],GO_CC[1:10,],GO_MF[1:10,])
      GO_hyper$group <- c(rep("BP",10), rep("CC",10), rep("MF",10))
      GO_hyper$direction <- c(rep("Hyper", 30))
    }
    #hypo
    if(T){
      entrez <- bitr(unique(na.omit(pie_dat_hypo$gene)), fromType = "SYMBOL",
                     toType = c( "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
      kegg <- enrichKEGG(entrez$ENTREZID,
                         organism = "hsa",
                         #keyType = "kegg",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         #universe,
                         minGSSize = 100,
                         maxGSSize = 5000,
                         qvalueCutoff = 0.05,
                         use_internal_data = T
      )
      kegg@result$p.adjust <- signif(kegg@result$p.adjust,2)
      GO_BP <-enrichGO(entrez$ENTREZID,
                       ont = "BP",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_CC <-enrichGO(entrez$ENTREZID,
                       ont = "CC",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_MF <-enrichGO(entrez$ENTREZID,
                       ont = "MF",OrgDb = "org.Hs.eg.db", minGSSize = 500, maxGSSize = 5000,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      GO_BP@result$p.adjust <- signif(GO_BP@result$p.adjust,2)
      GO_CC@result$p.adjust <- signif(GO_CC@result$p.adjust,2)
      GO_MF@result$p.adjust <- signif(GO_MF@result$p.adjust,2)
      kegg <- as.data.frame(kegg)
      GO_BP <- as.data.frame(GO_BP)
      GO_CC <- as.data.frame(GO_CC)
      GO_MF <- as.data.frame(GO_MF)
      kegg$logp <- -log10(kegg$p.adjust)
      kegg_hypo <- kegg[order(kegg$logp, decreasing = T),]
      kegg_hypo <- kegg_hypo[1:10,]
      GO_BP$logp <- -log10(GO_BP$p.adjust)
      GO_BP <- GO_BP[order(GO_BP$logp, decreasing = T),]
      GO_CC$logp <- -log10(GO_CC$p.adjust)
      GO_CC <- GO_CC[order(GO_CC$logp, decreasing = T),]
      GO_MF$logp <- -log10(GO_MF$p.adjust)
      GO_MF <- GO_MF[order(GO_MF$logp, decreasing = T),]
      GO_hypo <- rbind(GO_BP[1:10,],GO_CC[1:10,],GO_MF[1:10,])
      GO_hypo$group <- c(rep("BP",10), rep("CC",10), rep("MF",10))
      GO_hypo$direction <- c(rep("Hypo", 30))
      #GO_hypo$logp <- GO_hypo$logp
    }
    #ploting
    if(T){
      p_go_hyper <- ggplot(GO_hyper, aes(x = logp, y = Description, fill = group)) +
        geom_bar(stat = "identity", width = 0.9) +
        scale_y_discrete(limits = rev(GO_hyper$Description)) +
        #scale_x_continuous(limits = c(-100,100)) +
        theme(panel.background = element_blank(),
              axis.text.y.left = element_text(size = 12)) +
        labs(x = "-log10(adjustedPvalue)", fill = "", y = "") 
      
      
      p_go_hypo <- ggplot(GO_hypo, aes(x = logp, y = Description, fill = group)) +
        geom_bar(stat = "identity", width = 0.9) +
        scale_y_discrete(limits = rev(GO_hypo$Description)) +
        #scale_x_continuous(limits = c(-100,100)) +
        theme(panel.background = element_blank(),
              axis.text.y.left = element_text(size = 12)) +
        labs(x = "-log10(adjustedPvalue)", fill = "", y = "") 
      
      p_kegg_hyper <- ggplot(kegg_hyper, aes(x = logp, y = Description, fill = logp)) +
        geom_bar(stat = "identity", width = 0.9) +
        scale_y_discrete(limits = rev(kegg_hyper$Description)) +
        #scale_x_continuous(limits = c(-100,100)) +
        theme(panel.background = element_blank(),
              axis.text.y.left = element_text(size = 14)) +
        labs(x = "-log10(adjustedPvalue)", fill = "-log(adjustedP)", y = "") 
      p_kegg_hypo <- ggplot(kegg_hypo, aes(x = logp, y = Description, fill = logp)) +
        geom_bar(stat = "identity", width = 0.9) +
        scale_y_discrete(limits = rev(kegg_hypo$Description)) +
        #scale_x_continuous(limits = c(-100,100)) +
        theme(panel.background = element_blank(),
              axis.text.y.left = element_text(size = 14)) +
        labs(x = "-log10(adjustedPvalue)", fill = "-log(adjustedP)", y = "")
    }
  }
  if(T){
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go_hyper.pdf", width = 10,height = 7)
    p_go_hyper
    dev.off()
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go_hypo.pdf", width = 10,height = 7)
    p_go_hypo
    dev.off()
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_kegg_hyper.pdf", width = 8,height = 7)
    p_kegg_hyper
    dev.off()
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_kegg_hypo.pdf", width = 8,height = 7)
    p_kegg_hypo
    dev.off()
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_kegg.pdf", width = 8,height = 7)
    p_kegg
    dev.off()
    pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go.pdf", width = 10,height = 7)
    p_go
    dev.off()
    }
}


#pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go_bp.pdf", width = 8,height = 7)
#p_bp
#dev.off()
#pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go_cc.pdf", width = 8,height = 7)
#p_cc
#dev.off()
#pdf("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/output/dan_myth/tumor_vs_norm/go_kegg/uni_cgi_go_mf.pdf", width = 8,height = 7)
#p_mf
#dev.off()








