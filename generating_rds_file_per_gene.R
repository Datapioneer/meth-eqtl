#############################################################
#############################################################
##########---------probibility calculation-----------########
#############################################################
#############################################################

############--------------------------------------Data preparation---------------------############

#############################################################
##########----------obtain meth anno data#-----------########
#############################################################
library(rtracklayer)
library(dplyr)
library(do)
library(readr)
cpgisland_info <- read.csv("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/annotation_dat/CGI.20220904")
cpgisland_info <- as.data.frame(strsplit2(cpgisland_info$Probe_ID.Knowledgebase,"\tCGI;"))
cpgisland_info <- data.frame("Probe_ID" = cpgisland_info$V1, "CpGi" = cpgisland_info$V2 )
cpgisland_info <- cpgisland_info[!duplicated(cpgisland_info$Probe_ID),]
rownames(cpgisland_info) <- cpgisland_info$Probe_ID

hhm_info <- read.csv("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/annotation_dat/ChromHMM.20220303")
hhm_info <- as.data.frame(strsplit2(hhm_info$Probe_ID.Knowledgebase,"\tChromHMM;"))
hhm_info <- data.frame("Probe_ID" = hhm_info$V1, "ChromHMM" = hhm_info$V2 )
rownames(hhm_info) <- hhm_info$Probe_ID

strand_info <- read_tsv("/Users/huangyiwen/Documents/My project/ultimate plan/newline/HCC_bioinformatics/annotation_dat/HM450.hg38.manifest.gencode.v36.tsv.gz")
strand_info <- as.data.frame(strand_info)
strand_info <- subset(strand_info, select = c("probeID", "CpG_chrm", "CpG_beg","CpG_end","probe_strand"))
strand_info <- rename(strand_info, Probe_ID = probeID)
rownames(strand_info) <- strand_info$Probe_ID

meth_anno <- merge(strand_info, hhm_info, by = "Probe_ID", all = FALSE )
meth_anno <- merge(meth_anno, cpgisland_info, by = "Probe_ID" , all = FALSE)
meth_anno$pos <- meth_anno$CpG_beg
rownames(meth_anno) <- meth_anno$Probe_ID

rm(strand_info,cpgisland_info, hhm_info)
#############################################################
##########----------obtain meth beta data------------########
#############################################################
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step1-output_hbv_T_N.Rdata")
meth_beta_hbv <- myLoad
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/step1-output_hcv_T_N.Rdata")
meth_beta_hcv <- myLoad
rm(myLoad)

meth_beta_hbv <- t(meth_beta_hbv$beta[,meth_beta_hbv$pd$group =="T"])  
rownames(meth_beta_hbv) <- gsub("T-","",rownames(meth_beta_hbv))

meth_beta_hcv <- t(meth_beta_hcv$beta[,meth_beta_hcv$pd$group =="T"])  
rownames(meth_beta_hcv) <- gsub("T-","",rownames(meth_beta_hcv))

glmnet.cv
#############################################################
##########----------obtain gene annotation ----------########
#############################################################
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/bulk_rna_seq/expr_extraction.RData")
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/bulk_rna_seq/DEGs.RData")
expr_raw_hbv <- expr_raw[,substr(colnames(expr_raw),1,3) == "HBV"]
expr_raw_hcv <- expr_raw[,substr(colnames(expr_raw),1,3) == "HCV"]
dim(expr_raw_hbv)
dim(expr_raw_hcv)
rm(expr_raw)
gene_annotation <- rtracklayer::import("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/annotation_dat/gencode.v36.annotation.gtf")
gene_anno <- as.data.frame(gene_annotation)
gene_anno <- subset(gene_anno, gene_type == "protein_coding" & type == "gene")
gene_anno <- gene_anno$gene_id %in% rownames(DEG_HBV) %>%gene_anno[.,]
gene_anno$gene_id <- substr(gene_anno$gene_id,1,15)

rm(gene_annotation,DEG_HBV,DEG_HCV,DEGs_HBV,DEGs_HCV)
save.image("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/data_prep_cpg_prob.Rdata")


############--------------------------------------Data assemble---------------------############
load("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/data_prep_cpg_prob.Rdata")
######----------####
colnames(expr_raw_hbv) <- substr(colnames(expr_raw_hbv),5,16)
expr_raw_hbv <- expr_raw_hbv[,intersect(colnames(expr_raw_hbv),rownames(meth_beta_hbv))]
meth_beta_hbv <- meth_beta_hbv[intersect(colnames(expr_raw_hbv),rownames(meth_beta_hbv)),]
rownames(expr_raw_hbv) <- substr(rownames(expr_raw_hbv),1,15)
expr_cpm_hbv <- cpm(expr_raw_hbv, log = TRUE)

colnames(expr_raw_hcv) <- substr(colnames(expr_raw_hcv),5,16)
expr_raw_hcv <- expr_raw_hcv[,intersect(colnames(expr_raw_hcv),rownames(meth_beta_hcv))]
meth_beta_hcv <- meth_beta_hcv[intersect(colnames(expr_raw_hcv),rownames(meth_beta_hcv)),]
rownames(expr_raw_hcv) <- substr(rownames(expr_raw_hbv),1,15)
expr_cpm_hcv <- cpm(expr_raw_hcv, log = TRUE)

######----------######
rownames(meth_anno) <- meth_anno$Probe_ID
meth_anno <- meth_anno[intersect(meth_anno$Probe_ID, colnames(meth_beta_hbv)),]
meth_beta_hbv <- meth_beta_hbv[,intersect(meth_anno$Probe_ID, colnames(meth_beta_hbv))]
meth_beta_hcv <- meth_beta_hcv[,intersect(meth_anno$Probe_ID, colnames(meth_beta_hcv))]



for (i in 1:nrow(gene_anno)) {
  ##########-----------gen info---------##########
  gene <- gene_anno[i,"gene_id"]  #check for i = 885. 546
  geneChr <- as.character(gene_anno[i,"seqnames"]) %>% gsub("chr","",.)
  strand <- as.character(gene_anno[i,"strand"])
  startbp <- gene_anno[i,"start"]
  endbp <- gene_anno[i,"end"]
  geneLen <- gene_anno[i,"width"]/1000
  
  
  ##########-----------cpg info---------##########
  info <- subset(meth_anno, CpG_chrm == paste("chr", geneChr, sep = ""))
  if(strand == "+"){
    info <- subset(info, (startbp-info$pos )/1000 <= 500 & (endbp-info$pos)/1000 >= -500)
  }else{
    info <- subset(info, (endbp-info$pos )/1000 <= 500 & (startbp-info$pos)/1000 >= -500)
  }
  info$genebody <- ifelse((startbp-info$pos)*(endbp-info$pos)< 0, info$genebody <- 1, info$genebody <- 0)
  
  if(strand == "+"){
    info$dist1 <- (startbp-info$pos)/1000
    info$dist2 <- (endbp-info$pos)/1000
    info$dist <- ifelse(info$genebody == 1, 
                        info$dist <- (info$pos-startbp)/(geneLen*1000),
                        ifelse(abs(info$dist1)-abs(info$dist2) > 0,
                               info$dist2,
                               info$dist1))
  }else{
    info$dist1 <- (info$pos-startbp)/1000
    info$dist2 <- (info$pos-endbp)/1000
    info$dist <- ifelse(info$genebody == 1, 
                        info$dist <- (startbp-info$pos)/(geneLen*1000),
                        ifelse(abs(info$dist1)-abs(info$dist2) > 0,
                               info$dist2,
                               info$dist1))
  }
  info <- subset(info, select = c("Probe_ID", "pos", "genebody", "dist","ChromHMM","cgi"))
  info <- rename(info, probeid = Probe_ID, bpLocation = pos)
  if(nrow(info) > 1){
    info$gene <- gene
    
    ##########-----------exprdat---------##########
    exprdata <- expr_cpm_hcv[gene,]
    ##########-----------methdat---------##########
    methdata <- meth_beta_hcv[,info$probeid]
    #methdata <- as.data.frame(methdata)
    ##########-------spearman------------##########
    corr <- c()
    for (j in 1:ncol(methdata)) {
      corr <- c(corr,cor.test(methdata[,j],exprdata, method = "spearman")$estimate)
    }
    info$corr <- corr
    
    
    
    ###########----lasso------###########
    lasso_fit <- cv.glmnet(x=methdata, y=exprdata, alpha=1, family="gaussian", nfolds = 3)
    lasso_fit <-  glmnet(x=methdata, y=exprdata, alpha=1, family="gaussian", lambda = lasso_fit$lambda.min)
    info$lasso <- as.vector(lasso_fit[["beta"]])
    
    
    filepath <- paste0("~/Documents/My project/ultimate plan/newline/HCC_bioinformatics/rDat/meth/prob/",gene,".Rds", collapse = "")
    #list(info = info, gne = gene, geneChr = geneChr, strand = strand, startbp = startbp, endbp = endbp, geneLen = geneLen, methdata = methdata, exprdata = exprdata) %>%
    save(info, gene, geneChr, strand, startbp, endbp, geneLen, methdata, exprdata, file = filepath)
    
    
  }else{}
}


