library(matchSCore2)
library(reshape2)

tissues <- c("Lung", "Bone marrow", "Skeletal muscle", "Spleen", "Tongue", "Trachea")

###Collect markers
#Human
for (tis in tissues){
  Z_exp <- as.matrix(read.csv(paste0("Data/Expression/",tis,"/",tis,"_Z_scores.txt"), row.names=1, sep=""))
  markers <- list()
  for (ct in colnames(Z_exp)){
    markers[ct] <- list(names(which(Z_exp[,ct] >= 2)))
  }
  save(markers, file = paste0("Output/Markers/",tis,".Robj"))
}

#Mouse
for (tis in tissues){
  Z_exp <- as.matrix(read.csv(paste0("Data/Mouse/Expression/",tis,"/",tis,"_Z_scores.txt"), row.names=1, sep=""))
  markers <- list()
  for (ct in colnames(Z_exp)){
    markers[ct] <- list(names(which(Z_exp[,ct] >= 2)))
  }
  save(markers, file = paste0("Output/Mouse/Markers/",tis,".Robj"))
}

#Compare similarities (using matchSCore2) and create files per each
#tissue-pairwise comparison (For human inter-tissue similarities)
for (tis1 in tissues){
  load(paste0("Output/Markers/",tis1,".Robj"))
  tis1_mark <- markers
  for (tis2 in tissues){
    if (tis1 != tis2){
      load(paste0("Output/Markers/",tis2,".Robj"))
      tis2_mark <- markers
      m_score <- matchSCore2(tis1_mark, tis2_mark, ylab = tis1, xlab = tis2)
      ranks <- melt(m_score$JI.mat)
      ranks <- ranks[order(ranks[,"value"], decreasing = T),]
      colnames(ranks) <- c(tis1, tis2, "value")
      write.table(ranks, file = paste0("Output/Cell type similarities/",tis1,"_",tis2,"_ranks.tsv"), sep = "\t")
      }}}

#Dictionary of human-mouse genes (See HM_dictionary.py file)
gene_dic <- read.delim("Output\\Mouse\\gene_dic.tsv", header=FALSE, row.names = 1)

#Compare similarities (using matchSCore2) and create files per each tissue-pairs
#comparison (For human-mouse tissue similarities)
for (tis in tissues){
  load(paste0("Output/Mouse/Markers/",tis,".Robj"))
  markers_m <- list()
  for (ct in names(markers)){
    orthologs <- c()
    for (gene in markers[[ct]]){
      if (gene %in% rownames(gene_dic)){
        additionals <- strsplit(gene_dic[gene,1], ",")
        orthologs <- c(orthologs, additionals)
      }
    }
    orthologs <- unique(orthologs)
    markers_m[[ct]] <- orthologs
  }
  load(paste0("Output/Markers/",tis,".Robj"))
  m_score <- matchSCore2(markers, markers_m, xlab = "Mouse", ylab = "Human")
  ranks <- melt(m_score$JI.mat)
  ranks <- ranks[order(ranks[,"value"], decreasing = T),]
  colnames(ranks) <- c("Human", "Mouse", "value")
  write.table(ranks, file = paste0("Output/Mouse/Cell type similarities/", tis,"_ranks.tsv"), sep = "\t")
}
