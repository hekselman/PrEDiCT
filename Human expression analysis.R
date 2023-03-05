library(Seurat)
library(dplyr)
library(Matrix)
library(SeuratDisk)
library(ggplot2)

#Convert Tabula Sapiens dataset to Seurat object (Full data should be downloaded from Tabula Sapi)
tissues <- c("Bone_Marrow", "Lung", "Muscle", "Spleen", "Tongue", "Trachea")
for (tis in tissues){
  Convert(paste("TS_",tis,".h5ad",sep=""), paste("TS_",tis,".h5seurat",sep=""))
  TS <- LoadH5Seurat(paste("TS_", tis,".h5seurat",sep=""), assays="RNA")
  save(TS, file = paste("Data/Robjects/", tis, "_1Data.Robj",sep=""))
}

#Draw the data acquired from 10X technology
tissues <- c("Bone_Marrow", "Lung", "Muscle", "Spleen", "Tongue", "Trachea")
for (tis in tissues){
  load(file = paste("Data/Robjects/", tis, "_1Data.Robj",sep=""))
  TS <- SplitObject(TS, split.by = "method")
  TS <- TS$`10X`
  print(tis)
  print(table(as.factor(TS@meta.data$donor)))
  TS <- SetIdent(TS, value = "cell_ontology_class")
  save(TS, file = paste("Data/Robjects/", tis, "_2tenX.Robj",sep=""))
}

####Unify the identities of cell type names across samples
#Lung
load(file = "Data/Robjects/Lung_2tenX.Robj")
levels(as.factor(TS@active.ident))
TS <- RenameIdents(TS, "cd4-positive alpha-beta t cell" = "cd4-positive, alpha-beta t cell")
TS <- RenameIdents(TS, "cd8-positive alpha-beta t cell" = "cd8-positive, alpha-beta t cell")
save(TS, file = "Data/Robjects/Lung_2tenX.Robj")
#Trachea
load(file = "Data/Robjects/Trachea_2tenX.Robj")
levels(as.factor(TS@active.ident))
TS <- RenameIdents(TS, "tracheal goblet cell" = "goblet cell")
save(TS, file = "Data/Robjects/Trachea_2tenX.Robj")


tissues <- c("Bone_Marrow", "Lung", "Muscle", "Spleen", "Tongue", "Trachea")
for (tis in tissues){
  load(file = paste("Data/Robjects/", tis, "_2tenX.Robj",sep=""))
  TS <- NormalizeData(TS)
  temp <- matrix(nrow = length(row.names(TS@assays$RNA@data)), ncol = length(levels(TS@active.ident)))
  rownames(temp) <- row.names(TS@assays$RNA@data)
  colnames(temp) <- levels(TS@active.ident)
  
  #Per cluster, calculate percentage of cells expressing each gene (>=0.05)
  for (i in 1:length(levels(TS@active.ident))){
    identToUse <- levels(TS@active.ident)[i]
    print(identToUse)
    clusterobj <- subset(TS, idents = identToUse)
    clusterobj@assays$RNA@data@x[clusterobj@assays$RNA@data@x >= 0.05] <- 1
    clusterobj@assays$RNA@data@x[clusterobj@assays$RNA@data@x < 0.05] <- 0
    tempmatrix <- as.matrix((rowSums(clusterobj@assays$RNA@data, sparseResult=T) / length(colnames(x=clusterobj)) * 100), nrow = length(row.names(clusterobj@assays$RNA@data)))
    temp[,identToUse] <- tempmatrix[,1]
  }
  wrtr <- write.table(temp, file = paste("Data/Expression/",tis,"_PCT.txt",sep=""))   #Writes a file for percentage of cell expressing the gene
  pct10 <- temp
  pct10 <- pct10[!apply(pct10, 1, function(x) all(x < 10)), ] #Removes rows that all values < 10
  pct10genes <- row.names(pct10)
  LogAvgEx <- AverageExpression(TS, show.progress = T, features = pct10genes, return.seurat = F)
  expthres <- median(unlist(LogAvgEx$RNA)) #Median of average expression (Serves as the expression threshold)
  for (gene in row.names(LogAvgEx$RNA)){
    for (cluster in (levels(TS@active.ident)))
      if (pct10[gene, cluster] < 10){
        LogAvgEx$RNA[gene, cluster] <- 0
      }
  }
  LogAvgEx$RNA <- LogAvgEx$RNA[!apply(LogAvgEx$RNA, 1, function(x) all(x < expthres )),] #Removes genes expressed lower than the threshold
  thresgenes <- row.names(LogAvgEx$RNA)
  LogAvgEx <- AverageExpression(TS, show.progress = T, features = thresgenes, return.seurat = F)
  wrtr <- write.table(LogAvgEx$RNA, file = paste("Data/Expression/",tis,"_FinalLogAvgEx.txt",sep="")) #Writes a file with final table of average expression

  #Calculate Z scores
  pref.exp <- matrix(0, nrow = length(row.names(LogAvgEx$RNA)), ncol = length(levels(TS@active.ident)))
  pref.exp <-  data.frame(pref.exp)
  row.names(pref.exp) <- row.names(LogAvgEx$RNA)
  colnames(pref.exp) <- levels(TS@active.ident)
  for (gene in row.names(pref.exp)){
    mean_g <- mean(as.numeric(LogAvgEx$RNA[gene,]))
    sd_g <- sd(as.numeric(LogAvgEx$RNA[gene,]))
    for (cluster in c(1:length(levels(TS@active.ident)))){
      pref.exp[gene, cluster] <- ((LogAvgEx$RNA[gene, cluster] - mean_g) / sd_g)
    }}
  wrtr <- write.table(pref.exp, file = paste("Data/Expression/", tis, "_Z_scores.txt",sep ="")) #save Z-score transformed values
}