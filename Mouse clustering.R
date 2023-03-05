library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(ranger)
library(sgof)

counts <- readRDS(file = "Data/Mouse/TM_droplet_mat.rds") #Counts file from Tabula Muris

###Choose a tissue to run (metadata files from Tabula Muris)
#Limb_Muscle
metadata <- read.csv("Data/Mouse/Metadata/Limb_Muscle.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
table(tm@meta.data$mouse.id)
#Lung
metadata <- read.csv("Data/Mouse/Metadata/Lung.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
tm <- subset(tm, subset = mouse.id == "3-F-56" | mouse.id == "3-M-5/6")
table(tm@meta.data$mouse.id)
#Marrow
metadata <- read.csv("Data/Mouse/Metadata/Marrow.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
table(tm@meta.data$mouse.id)
#Spleen
metadata <- read.csv("Data/Mouse/Metadata/Spleen.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
table(tm@meta.data$mouse.id)
#Tongue
metadata <- read.csv("Data/Mouse/Metadata/Tongue.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
table(tm@meta.data$mouse.id)
#Trachea
metadata <- read.csv("Data/Mouse/Metadata/Trachea.csv", row.names = 1, header = T)
tm <- CreateSeuratObject(counts = counts, meta.data = metadata)
table(tm@meta.data$mouse.id)

###Run per tissue
tm@meta.data <- droplevels(tm@meta.data)
tm <- subset(tm, cells = rownames(tm@meta.data[!is.na(tm@meta.data$cell_ontology_id),]))
rm(counts, metadata)
tm <- NormalizeData(tm)
tm <- FindVariableFeatures(tm, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.05,3))
VariableFeaturePlot(tm)
tm <- ScaleData(tm, vars.to.regress = "mouse.id")
tm <- RunPCA(tm, features = VariableFeatures(object = tm), npcs = 70)
tm <- JackStraw(tm, dims = 70)
tm <- ScoreJackStraw(tm, dims = 1:70)
ElbowPlot(tm, ndims = 70)
JackStrawPlot(tm, dims = 1:70)

#Numebr of significant principal components of each tissue
PCs <- 53 #Limb_Muscle
PCs <- 32 #Lung
PCs <- 25 #Marrow
PCs <- 37 #Spleen
PCs <- 45 #Tongue
PCs <- 62 #Trachea

###Run per tissue
tm <- RunUMAP(tm, dims = 1:PCs)
tm <- FindNeighbors(tm, dims = 1:PCs)
tm <- FindClusters(tm)
DimPlot(tm, reduction = "umap")
tm <- BuildClusterTree(tm, reorder = T, reorder.numeric = T, dims = 1:PCs)
PlotClusterTree(tm)

###Save the relevant tissue
save(tm, file = "Data/Mouse/Robjects/Limb_Muscle_1BasicCluster.Robj") #Limb_Muscle
save(tm, file = "Data/Mouse/Robjects/Lung_1BasicCluster.Robj") #Lung
save(tm, file = "Data/Mouse/Robjects/Marrow_1BasicCluster.Robj") #Marrow
save(tm, file = "Data/Mouse/Robjects/Spleen_1BasicCluster.Robj") #Spleen
save(tm, file = "Data/Mouse/Robjects/Tongue_1BasicCluster.Robj") #Tongue
save(tm, file = "Data/Mouse/Robjects/Trachea_1BasicCluster.Robj") #Trachea

#Functions for assessing nodes and splits of hierarchical tree
AssessSplit <- function(
  object,
  node,
  cluster1,
  cluster2,
  genes.training = NULL,
  print.output = TRUE,
  ...
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (!missing(x = node)) {
    if (!missing(x = cluster1) || !missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (!node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, ident = group1)
  group2.cells <- WhichCells(object = object, ident = group2)
  assess.data <- subset(
    x = object,
    cells = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group1.cells,
    value = "g1"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group2.cells,
    value = "g2"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    # training.genes = assess.data@var.genes,
    training.genes = genes.training,
    training.classes = assess.data@active.ident,
    ...
  )
  oobe <- rfc$prediction.error
  if (print.output) {
    message(paste0("Out of Bag Error: ", round(x = oobe, digits = 4) * 100, "%"))
  }
  return(oobe)
}

AssessNodes <- function(
  object,
  node.list,
  all.below = FALSE,
  genes.training = NULL
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (missing(x = node.list)) {
    node.list <- GetAllInternalNodes(tree = tree)
  } else {
    possible.nodes <- GetAllInternalNodes(tree = tree)
    if (any(!node.list %in% possible.nodes)) {
      stop(paste(
        node.list[!(node.list %in% possible.nodes)],
        "not valid internal nodes"
      ))
    }
    if (length(x = node.list == 1) && all.below) {
      node.list <- c(node.list, DFT(tree = tree, node = node.list))
    }
  }
  oobe <- pbsapply(
    X = node.list,
    FUN = function(x) {
      return(AssessSplit(
        object = object,
        node = x,
        genes.training = genes.training,
        print.output = FALSE,
        verbose = FALSE
      ))
    }
  )
  return(data.frame(node = node.list, oobe))
}

SetIfNull <- function(x, default) {
  if (is.null(x = x)) {
    return(default)
  } else {
    return(x)
  }
}

GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if(! only.children){
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <-c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (! only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}

BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  PackageCheck('ranger')
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = object@assays$RNA@data)
  )
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = as.matrix(object@assays$RNA@data[training.genes, ])
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}

MergeNode <- function(object, node.use, rebuild.tree = FALSE, ...) {
  object.tree <- object@tools$BuildClusterTree
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- intersect(x = node.children, y = levels(x = object@active.ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells = children.cells,
      value = min(node.children)
    )
  }
  if (rebuild.tree) {
    object <- BuildClusterTree(object = object, ...)
  }
  return(object)
}

###Load the relevant tissue
load("Data/Mouse/Robjects/Limb_Muscle_1BasicCluster.Robj") #Limb_Muscle
load("Data/Mouse/Robjects/Lung_1BasicCluster.Robj") #Lung
load("Data/Mouse/Robjects/Marrow_1BasicCluster.Robj") #Marrow
load("Data/Mouse/Robjects/Spleen_1BasicCluster.Robj") #Spleen
load("Data/Mouse/Robjects/Tongue_1BasicCluster.Robj") #Tongue
load("Data/Mouse/Robjects/Trachea_1BasicCluster.Robj") #Trachea

###Run per tissue
node.scores <- AssessNodes(tm, genes.training = VariableFeatures(tm))
node.scores <- node.scores[order(node.scores$oobe,decreasing = T),]
node.scores

###Number of splits to merge per tissue
nodes.merge <- node.scores[1:2,] #Limb_Muscle
nodes.merge <- node.scores[1:4,] #Lung
nodes.merge <- node.scores[1:1,] #Marrow
nodes.merge <- node.scores[1:2,] #Spleen
nodes.merge <- node.scores[1:1,] #Spleen 2nd round
nodes.merge <- node.scores[1:1,] #Tongue
nodes.merge <- node.scores[1:4,] #Trachea
nodes.merge <- node.scores[1:1,] #Trachea 2nd round

###Run per tissue
nodes.to.merge <- sort(nodes.merge$node, decreasing = T)
for (n in nodes.to.merge){
  tm <- MergeNode(tm, n)}
tm <- BuildClusterTree(tm, reorder = T, reorder.numeric = T, dims = 1:PCs)
PlotClusterTree(tm)
#Repeat AssessNodes until stable (p<0.05)

###Save the relevant tissue
save(tm, file = "Data/Mouse/Robjects/Limb_Muscle_2Merging.Robj") #Limb_Muscle
save(tm, file = "Data/Mouse/Robjects/Lung_2Merging.Robj") #Lung
save(tm, file = "Data/Mouse/Robjects/Marrow_2Merging.Robj") #Marrow
save(tm, file = "Data/Mouse/Robjects/Spleen_2Merging.Robj") #Spleen
save(tm, file = "Data/Mouse/Robjects/Tongue_2Merging.Robj") #Tongue
save(tm, file = "Data/Mouse/Robjects/Trachea_2Merging.Robj") #Trachea

tissues <- c("Limb_Muscle", "Lung", "Marrow",
             "Spleen", "Tongue", "Trachea")

#Check whether tissues have significant differences in their distribution of cells between batches (BH corrected)
chis <- c()
for (tis in tissues){
  load(file = paste("Data/Mouse/Robjects/", tis, "_2Merging.Robj", sep = ""))
  batch <- table(tm@meta.data[c("mouse.id","tree.ident")])
  chi<-chisq.test(batch)
  chi<-chi$p.value
  chis <- c(chis,chi)
}
chis_fdr <- BH(chis, alpha = 0.05)
print(chis_fdr$Adjusted.pvalues)

#Check which clusters are not evenly distributed between batches (BH corrected)
#Writes a table with a post hoc analysis for the chi2 test above
for (tis in tissues){
  load(file = paste("Data/Mouse/Robjects/", tis, "_2Merging.Robj", sep = ""))
  batch <- table(tm@meta.data[c("mouse.id","tree.ident")])
  pairs <- c()
  chis <- c()
  prop <- c()
  for (row in 1:nrow(batch)){
    for (col in 1:ncol(batch)){
      entry <- batch[row, col]
      rest_col <- sum(batch[,col]) - entry
      rest_row <- sum(batch[row,]) - entry
      rest <- sum(batch) - rest_col - rest_row - entry
      batch_compare <- data.frame(rows = c(entry,rest_row), cols = c(rest_col,rest))
      chi <- chisq.test(batch_compare)
      chi <- chi$p.value
      pairs <- c(pairs,paste0(row, "_", col, sep=""))
      chis <- c(chis,chi)
      prop <- c(prop,entry / sum(batch[row,]))
    }
  }
  correction <- data.frame(pair = pairs, pval = chis, prop = prop)
  correction <- correction[order(correction$pval, decreasing = F),]
  p_batch_fdr <- BH(correction$pval, alpha = 0.05)
  correction["fdr"] <- p_batch_fdr$Adjusted.pvalues
  write.table(correction, file = paste("Output/Mouse/Batches/", tis, ".tsv",sep=""), sep = "\t")
}

#Merge similar clusters that came from different batches (found only in lung and tongue)
#Lung
load(file = "Data/Mouse/Robjects/Lung_2Merging.Robj")
PCs <- 32
tm <- MergeNode(tm, 25)
tm <- BuildClusterTree(tm, reorder = T, reorder.numeric = T, dims = 1:PCs)
PlotClusterTree(tm)
node.scores <- AssessNodes(tm, genes.training = VariableFeatures(tm))
node.scores <- node.scores[order(node.scores$oobe,decreasing = T),]
node.scores
batch <- table(tm@meta.data[c("mouse.id","tree.ident")])
pairs <- c()
chis <- c()
prop <- c()
for (row in 1:nrow(batch)){
  for (col in 1:ncol(batch)){
    entry <- batch[row, col]
    rest_col <- sum(batch[,col]) - entry
    rest_row <- sum(batch[row,]) - entry
    rest <- sum(batch) - rest_col - rest_row - entry
    batch_compare <- data.frame(rows = c(entry,rest_row), cols = c(rest_col,rest))
    chi <- chisq.test(batch_compare)
    chi <- chi$p.value
    pairs <- c(pairs,paste0(row, "_", col, sep=""))
    chis <- c(chis,chi)
    prop <- c(prop,entry / sum(batch[row,]))
  }
}
correction <- data.frame(pair = pairs, pval = chis, prop = prop)
correction <- correction[order(correction$pval, decreasing = F),]
p_batch_fdr <- BH(correction$pval, alpha = 0.05)
correction["fdr"] <- p_batch_fdr$Adjusted.pvalues
write.table(correction, file = paste("Output/Mouse/Batches/Lung0.tsv",sep=""), sep = "\t")
save(tm, file = "Data/Mouse/Robjects/Lung_3BatchEffect.Robj")

#Tongue
load(file = "Data/Mouse/Robjects/Tongue_2Merging.Robj")
PCs <- 45
tm <- MergeNode(tm, 19)
tm <- MergeNode(tm, 26)
tm <- BuildClusterTree(tm, reorder = T, reorder.numeric = T, dims = 1:PCs)
PlotClusterTree(tm)
node.scores <- AssessNodes(tm, genes.training = VariableFeatures(tm))
node.scores <- node.scores[order(node.scores$oobe,decreasing = T),]
node.scores
batch <- table(tm@meta.data[c("mouse.id","tree.ident")])
pairs <- c()
chis <- c()
prop <- c()
for (row in 1:nrow(batch)){
  for (col in 1:ncol(batch)){
    entry <- batch[row, col]
    rest_col <- sum(batch[,col]) - entry
    rest_row <- sum(batch[row,]) - entry
    rest <- sum(batch) - rest_col - rest_row - entry
    batch_compare <- data.frame(rows = c(entry,rest_row), cols = c(rest_col,rest))
    chi <- chisq.test(batch_compare)
    chi <- chi$p.value
    pairs <- c(pairs,paste0(row, "_", col, sep=""))
    chis <- c(chis,chi)
    prop <- c(prop,entry / sum(batch[row,]))
  }
}
correction <- data.frame(pair = pairs, pval = chis, prop = prop)
correction <- correction[order(correction$pval, decreasing = F),]
p_batch_fdr <- BH(correction$pval, alpha = 0.05)
correction["fdr"] <- p_batch_fdr$Adjusted.pvalues
write.table(correction, file = paste("Output/Mouse/Batches/Tongue0.tsv",sep=""), sep = "\t")
save(tm, file = "Output/Mouse/Robjects/Tongue_3BatchEffect.Robj")
