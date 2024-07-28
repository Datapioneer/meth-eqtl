##############################################################
########################################  ####################
#######################################     ##################
#######################################      #################
##########       meth analysis                ################
#######################################      #################
#######################################     ##################
########################################  ####################
##############################################################

library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(tidyverse)
library(readr)
library(do)
library(magrittr)
library(reshape2)
library(impute)

####################################################################
#######-------match meth data and clinical dat.    -------##########
####################################################################
##### match meth data and clinical dat
load("~/project2/tcga_meth/rDat/TCGA_meth.RData")
load("~/project2/tcga_meth/rDat/viro_stage_info.RData")
rownames(viro_info) <- viro_info$bcr_patient_barcode
dat_meth_sig <- data.frame(assay(dat_meth_tumor))
gsub("\\.","-",colnames(dat_meth_sig)) %>% substr(.,1,12) -> colnames(dat_meth_sig)
dat_meth_sig <- intersect(colnames(dat_meth_sig), rownames(viro_info)) %>% dat_meth_sig[,.]
viro_info <- intersect(colnames(dat_meth_sig), rownames(viro_info)) %>% viro_info[.,]

a <- dat_meth_sig
b <- viro_info
a[1:4,1:4]
beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))

# 一定要保证，甲基化信号值矩阵，和表型信息，是一一对应的


####################################################################
##########-------meth analysis DMP----------------------############
####################################################################
if(T){
  library(ChAMP)
  # beta 信号值矩阵里面不能有NA值
  myLoad=champ.filter(beta = a,pd = b)
  save(myLoad,file = '~/project2/tcga_meth/rDat/step1-output.Rdata')
  #CpG.GUI(CpG = rownames(myLoad$beta), arraytype = "450K")
  champ.QC(beta = myLoad$beta,
           pheno=group_list,
           mdsPlot=TRUE,
           densityPlot=TRUE,
           dendrogram=TRUE,
           PDFplot=TRUE,
           Rplot=TRUE,
           Feature.sel="None",
           resultsDir="./CHAMP_QCimages/")
  #step 1
  myLoad  # 存储了甲基化信号矩阵和表型信息。
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm)
  group_list=myLoad$pd$history_hepato_carcinoma_risk_factors
  table(group_list)
  myDMP <- champ.DMP(beta = myNorm,pheno=group_list) #champ workflow
  head(myDMP[[1]])
  save(myDMP,file = '~/project2/tcga_meth/rDat/step3-output-myDMP.Rdata')
  
}
####################################################################
##########-------meth analysis DMR----------------------############
####################################################################
if(F){
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  #DMR.GUI(DMR=myDMR)
  
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  head(myBlock$Block)
  #Block.GUI(Block=myBlock,beta = myNorm,pheno=group_list,
  #          runDMP=TRUE,compare.group=NULL,arraytype="450K")
  
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  #myEpiMod <- champ.EpiMod(beta=myNorm,pheno=group_list)
  save(myDMR,file = '~/project2/tcga_meth/rDat/step3-output-myDMR.Rdata')
  save(myGSEA,file = '~/project2/tcga_meth/rDat/step3-output-myGSEA.Rdata')
  
}

#####################################
######------volcano plot-----########
#####################################
if(T){
  DEG <- myDMP$HCV_to_HBV
  # step 1 setting threshold
  beta_diff_cutoff <- 0.05
  DEG$result = as.factor(ifelse(DEG$adj.P.Val < 1e-3 & abs(DEG$deltaBeta) > beta_diff_cutoff,
                                ifelse(DEG$deltaBeta > 0 ,'UP','DOWN'),'NOT'))
  table(DEG$result)
  # step 2 setting title
  this_tile <- paste0('Cutoff for Different Methylation is ',round(beta_diff_cutoff,3), #round保留小数位数
                      '\nThe number of up prob is ',nrow(DEG[DEG$result =='UP',]) ,
                      '\nThe number of down prob is ',nrow(DEG[DEG$result =='DOWN',])
  )
  
  # step 3
  library(ggplot2)
  ggplot(data = DEG, aes(x = deltaBeta, y = -log10(adj.P.Val), color=result)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("Different Methylation") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## 这里要注意和之前设置的result三个因子相对应，DOWN就设为blue，NOT就设为black
  DEGs <- DEG[DEG$adj.P.Val < 1e-3,]
  DEGs <- DEGs[abs(DEGs$deltaBeta) > beta_diff_cutoff,]
  DEGs <- subset(DEG,adj.P.Val < 1e-3 & abs(deltaBeta) > beta_diff_cutoff )
  DEGs <- DEGs["result"]
  
}


################################################################
##########------DMP/R/block HBV tumor vs normal-----############
################################################################
if(T){
  dat_meth_sig <- dat_meth_sig[,viro_info$history_hepato_carcinoma_risk_factors == "HBV"]
  colnames(dat_meth_sig) <- paste("T",colnames(dat_meth_sig), sep = "-")
  dat_meth_sig_norm <- data.frame(assay(dat_meth_norm))
  gsub("\\.","-",colnames(dat_meth_sig_norm)) %>% substr(.,1,12) -> colnames(dat_meth_sig_norm)
  colnames(dat_meth_sig_norm) <- paste("N",colnames(dat_meth_sig_norm), sep = "-")
  dat_meth_sig <- cbind(dat_meth_sig, dat_meth_sig_norm)
  a <- dat_meth_sig
  b <- data.frame("group" = substr(colnames(dat_meth_sig),1,1),"tcga_id" = colnames(dat_meth_sig))
  rownames(b) <- b$tcga_id
  a[1:4,1:4]
  beta=as.matrix(a)
  beta=impute.knn(beta)
  betaData=beta$data
  betaData=betaData+0.00001
  a=betaData
  a[1:4,1:4]
  identical(colnames(a),rownames(b))
  library(ChAMP)
  # beta 信号值矩阵里面不能有NA值
  myLoad=champ.filter(beta = a,pd = b)
  save(myLoad,file = '~/project2/tcga_meth/rDat/step1-output_hbv_T_N.Rdata')
  group_list=b$group
  table(group_list)
  champ.QC(beta = myLoad$beta,
           pheno=group_list,
           mdsPlot=TRUE,
           densityPlot=TRUE,
           dendrogram=TRUE,
           PDFplot=TRUE,
           Rplot=TRUE,
           Feature.sel="None",
           resultsDir="~/project2/tcga_meth/output/tumor_vs_norm/hbv/")
  #step 1
  myLoad  # 存储了甲基化信号矩阵和表型信息。
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm)
  myDMP <- champ.DMP(beta = myNorm,pheno=group_list) #champ workflow
  head(myDMP[[1]])
  save(myDMP,file = '~/project2/tcga_meth/rDat/step3-output-myDMP_hbv_T_N.Rdata')
  
  ############----------meth analysis DMR--------################
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  save(myDMR,file = '~/project2/tcga_meth/rDat/step3-output-myDMR_hbv_T_N.Rdata')
  save(myGSEA,file = '~/project2/tcga_meth/rDat/step3-output-myGSEA_hbv_T_N.Rdata')
  
}


################################################################
##########------DMP/R/block HCV tumor vs normal-----############
################################################################
if(T){
  dat_meth_sig <- dat_meth_sig[,viro_info$history_hepato_carcinoma_risk_factors == "HCV"]
  colnames(dat_meth_sig) <- paste("T",colnames(dat_meth_sig), sep = "-")
  dat_meth_sig_norm <- data.frame(assay(dat_meth_norm))
  gsub("\\.","-",colnames(dat_meth_sig_norm)) %>% substr(.,1,12) -> colnames(dat_meth_sig_norm)
  colnames(dat_meth_sig_norm) <- paste("N",colnames(dat_meth_sig_norm), sep = "-")
  dat_meth_sig <- cbind(dat_meth_sig, dat_meth_sig_norm)
  a <- dat_meth_sig
  b <- data.frame("group" = substr(colnames(dat_meth_sig),1,1),"tcga_id" = colnames(dat_meth_sig))
  rownames(b) <- b$tcga_id
  a[1:4,1:4]
  beta=as.matrix(a)
  beta=impute.knn(beta)
  betaData=beta$data
  betaData=betaData+0.00001
  a=betaData
  a[1:4,1:4]
  identical(colnames(a),rownames(b))
  library(ChAMP)
  # beta 信号值矩阵里面不能有NA值
  myLoad=champ.filter(beta = a,pd = b)
  save(myLoad,file = '~/project2/tcga_meth/rDat/step1-output_hcv_T_N.Rdata')
  group_list=b$group
  table(group_list)
  champ.QC(beta = myLoad$beta,
           pheno=group_list,
           mdsPlot=TRUE,
           densityPlot=TRUE,
           dendrogram=TRUE,
           PDFplot=TRUE,
           Rplot=TRUE,
           Feature.sel="None",
           resultsDir="~/project2/tcga_meth/output/tumor_vs_norm/hcv/")
  #step 1
  myLoad  # 存储了甲基化信号矩阵和表型信息。
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm)
  myDMP <- champ.DMP(beta = myNorm,pheno=group_list) #champ workflow
  head(myDMP[[1]])
  save(myDMP,file = '~/project2/tcga_meth/rDat/step3-output-myDMP_hcv_T_N.Rdata')
  
  ############----------meth analysis DMR--------################
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  save(myDMR,file = '~/project2/tcga_meth/rDat/step3-output-myDMR_hcv_T_N.Rdata')
  save(myGSEA,file = '~/project2/tcga_meth/rDat/step3-output-myGSEA_hcv_T_N.Rdata')
  
}



################################################################
###########----------DMP distribution--------------#############
################################################################
DEG <- myDMP$T_to_N
beta_diff_cutoff <- 0.05
DEG$result = as.factor(ifelse(DEG$adj.P.Val < 1e-3 & abs(DEG$deltaBeta) > beta_diff_cutoff,
                              ifelse(DEG$deltaBeta > 0 ,'UP','DOWN'),'NOT'))
DEGs <- DEG[DEG$result != "UP",]

###########----------pie chart--------------#############
if(T){
  numbers <- as.vector(table(DEGs$feature))
  name <- names(table(DEGs$feature))
  percent <- 100*(numbers/sum(numbers)) %>% round(.,4) 
  percent <- paste(percent,"%",sep = "")
  pdf(file = "~/project2/tcga_meth/output/tumor_vs_norm/hcv/pie_t_n_feat_up.pdf", width = 8,height = 8)
  pie(numbers, labels = percent, col = rainbow(length(percent)))
  legend("topleft", name, fill = rainbow(length(percent)), cex = 0.8)
  dev.off()
}
DEGs_hbv %>% head()
DEGs_hcv %>% head()
DEGs <- DEGs_hbv
HBV_Hypermethylated <- DEGs_hbv[DEGs_hbv$result =="UP",]$feature %>% table() %>% proportions() %>% data.frame(.) %>% .$Freq
HCV_Hypermethylated <- DEGs_hcv[DEGs_hcv$result =="UP",]$feature %>% table() %>% proportions() %>% data.frame(.) %>% .$Freq
HBV_Hypomethylated <- DEGs_hbv[DEGs_hbv$result =="DOWN",]$feature %>% table() %>% proportions() %>% data.frame(.) %>% .$Freq
HCV_Hypomethylated <- DEGs_hcv[DEGs_hcv$result =="DOWN",]$feature %>% table() %>% proportions() %>% data.frame(.) %>% .$Freq
dat_for_bar <- cbind("HBV_Hypermethylated"= HBV_Hypermethylated,
                     "HCV_Hypermethylated" = HCV_Hypermethylated,
                     "HBV_Hypomethylated" = HBV_Hypomethylated,
                     "HCV_Hypomethylated" =HCV_Hypomethylated) %>% data.frame()
rownames(dat_for_bar) <- c("1stExon", "3'UTR",   "5'UTR" ,  "Body",    "IGR", "TSS1500", "TSS200" )
dat_for_bar
###########----------bar plot--------------#############


axis(2,c("HBV_Hypermethylated","HCV_Hypermethylated", "HBV_Hypomethylated","HCV_Hypomethylated"))

