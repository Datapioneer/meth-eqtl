a <- as.data.frame(gene_annotation)
a <- subset(a, gene_type == "protein_coding" & type == "transcript")
b <- a$gene_id %in% names(table(a$gene_id) [table(a$gene_id) == 1]) %>% a[.,] # single isoform

c <- a$gene_id %in% names(table(a$gene_id) [table(a$gene_id) > 1]) %>% a[.,] # multiple isoforms with same tss and tes

d <- data.frame()
for (i in 1:length(unique(c$gene_id))) {
  transcripts <- c[c$gene_id %in% unique(c$gene_id)[i], ]
  same_start <- c()
  same_end <- c()
  for (j in 2:nrow(transcripts)) {
    same_start <- c(same_start, identical(transcripts[1, "start"], transcripts[j,"start"]))
    same_end <-  c(same_start, identical(transcripts[1, "end"], transcripts[j,"end"]))
  }
  if(FALSE %in% same_start | FALSE %in% same_end){
    d <- d
  }else{
    d <- rbind(d, transcripts)
  }
   
}
gene_ens <- c(b$gene_id,d$gene_id)

