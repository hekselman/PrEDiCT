library(Seurat)
library(Matrix)

tissues <- c("Limb_Muscle", "Lung", "Marrow",
             "Spleen", "Tongue", "Trachea")

#Before running this code, make sure to set the active identity of
#the mouse cells to the new annotations as supplied in the
#supplementary data of the study
for (tis in tissues){
  load(file = paste("Data/Mouse/Robjects/", tis, "_2Merging.Robj", sep = ""))
  if (file.exists(paste("Data/Mouse/Robjects/", tis, "_3BatchEffect.Robj", sep = ""))){
    load(file = paste("Data/Mouse/Robjects/", tis, "_3BatchEffect.Robj", sep = ""))
  }
  temp <- matrix(nrow = length(row.names(tm@assays$RNA@data)), ncol = length(levels(tm@active.ident)))
  rownames(temp) <- row.names(tm@assays$RNA@data)
  colnames(temp) <- levels(tm@active.ident)
  
  #Per cluster, calculate percentage of cells expressing each gene (>=0.05)
  for (i in 1:length(levels(tm@active.ident))){
    identToUse <- levels(tm@active.ident)[i]
    print(identToUse)
    clusterobj <- subset(tm, idents = identToUse)
    clusterobj@assays$RNA@data@x[clusterobj@assays$RNA@data@x >= 0.05] <- 1
    clusterobj@assays$RNA@data@x[clusterobj@assays$RNA@data@x < 0.05] <- 0
    temp_table <- Matrix(clusterobj@assays$RNA@data, sparse = T)
    tempmatrix <- as.matrix((rowSums(temp_table) / length(colnames(x=clusterobj)) * 100), nrow = length(row.names(temp_table)))
    temp[,identToUse] <- tempmatrix[,1]
  }
  wrtr <- write.table(temp, file = paste("Data/Mouse/Expression/", tis, "/",tis,"_PCT.txt",sep=""))  #Writes a file for percentage of cell expressing the gene
  pct10 <- temp
  pct10 <- pct10[!apply(pct10, 1, function(x) all(x < 10)), ] #Removes rows that all values < 10
  pct10genes <- row.names(pct10)
  LogAvgEx <- AverageExpression(tm, show.progress = T, features = pct10genes, return.seurat = F)
  expthres <- median(unlist(LogAvgEx$RNA)) #Median of average expression (Serves as the expression threshold)
    for (gene in row.names(LogAvgEx$RNA)){
    for (cluster in (levels(tm@active.ident)))
      if (pct10[gene, cluster] < 10){
        LogAvgEx$RNA[gene, cluster] <- 0
      }
  }
  LogAvgEx$RNA <- LogAvgEx$RNA[!apply(LogAvgEx$RNA, 1, function(x) all(x < expthres )),] #Removes genes expressed lower than the threshold
  thresgenes <- row.names(LogAvgEx$RNA)
  LogAvgEx <- AverageExpression(tm, show.progress = T, features = thresgenes, return.seurat = F)
  wrtr <- write.table(LogAvgEx$RNA, file = paste("Data/Mouse/Expression/", tis, "/", tis, "_FinalLogAvgEx.txt", sep= "")) #Writes a file with final table of average expression
  
  #Calculation of Z scores of average gene expression
  pref.exp <- matrix(0, nrow = length(row.names(LogAvgEx$RNA)), ncol = length(levels(tm@active.ident)))
  pref.exp <-  data.frame(pref.exp)
  row.names(pref.exp) <- row.names(LogAvgEx$RNA)
  colnames(pref.exp) <- levels(tm@active.ident)
  for (gene in row.names(pref.exp)){
    gene_mean <- mean(as.numeric(LogAvgEx$RNA[gene,]))
    gene_sd <- sd(as.numeric(LogAvgEx$RNA[gene,]))
    for (cluster in c(1:length(levels(tm@active.ident)))){
      pref.exp[gene, cluster] <- ((LogAvgEx$RNA[gene, cluster] - gene_mean) / gene_sd)
    }}
  wrtr <- write.table(pref.exp, file = paste("Data/Mouse/Expression/", tis, "/", tis, "_Z_scores.txt", sep= "")) #Writes a file for Z-score transformed values
}
