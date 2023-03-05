### Success function
success_test <- function(pref_cts, simi_cts){
  if (length(intersect(pref_cts, simi_cts)) > 0){
    run_success <- 1
  } else {
    run_success <- 0
  }
  return(run_success)
}

#General
tissues <- c("Lung", "Bone marrow", "Skeletal muscle",
             "Spleen", "Tongue", "Trachea")
p_vals <- c()


### Testing similarities between affected cell types across human tissues
tissues_bar <- c()
for (tis1 in tissues) {
  for (tis2 in tissues) {
    if (tis1 != tis2) {
      tissues_bar <- c(tissues_bar, paste0(tis1,"_",tis2))
      ranks <- read.delim(paste0("Output/Cell type similarities/",tis1,"_",tis2,"_ranks.tsv"), row.names=1)
      ranks <- ranks[c(1:floor(length(rownames(ranks))/10)),]
      ranks <- ranks[which(ranks[,"value"] >= 0.05),]
      score1 <- read.delim(paste0("Output/PrEDiCT score/",tis1,".tsv"))
      score2 <- read.delim(paste0("Output/PrEDiCT score/",tis2,".tsv"))
      ### Random ###
      repetitions <- 1000
      total_success <- c()
      for (r in c(1:repetitions)){
        run_success <- c()
        dis_num <- 0
        for (dis in unique(score1$Disease.Phenotypic.series.accession..OMIM.)){
          pref_cts1 <- score1[score1$PrEDiCT.score >= 2 & score1$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
          if (length(pref_cts1) > 0){
            dis_num <- dis_num + 1
            simi_cts2 <- ranks[ranks[,gsub(" ",".",tis1)] %in% pref_cts1,gsub(" ",".",tis2)]
            pref_cts2 <- score2[score2$PrEDiCT.score >= 2 & score2$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
            random_cts2 <- sample(unique(score2$Cell.type), length(pref_cts2))
            run_success <- c(run_success, success_test(pref_cts = random_cts2, simi_cts = simi_cts2))
          }
        }
        num_success <- sum(run_success)
        total_success <- c(total_success, num_success)
      }
      ### Real data ###
      run_success <- c()
      for (dis in unique(score1$Disease.Phenotypic.series.accession..OMIM.)){
        pref_cts1 <- score1[score1$PrEDiCT.score >= 2 & score1$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
        if (length(pref_cts1) > 0){
          simi_cts2 <- ranks[ranks[,gsub(" ",".",tis1)] %in% pref_cts1,gsub(" ",".",(tis2))]
          pref_cts2 <- score2[score2$PrEDiCT.score >= 2 & score2$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
          run_success <- c(run_success, success_test(pref_cts = pref_cts2, simi_cts = simi_cts2))
        }
      }
      num_success <- sum(run_success)
      real_success <- num_success
      
      #### Check P value
      p <- sum(total_success >= real_success)/repetitions
      p_vals <- c(p_vals, p)
      print(paste0(tis1, "_", tis2))
      #write.table(data.frame(real_success,total_success), file = paste0("Permutation_success\\Success_", tis1, "_", tis2, ".tsv"), sep = "\t")
    }
  }
}
print(FDR(pvalues = data.frame(tis = tissues_bar,ps = p_vals)))


### Testing similarities between affected cell types in human versus mouse tissues
for (tis in tissues){
  ranks <- read.delim(paste0("Output/Mouse/Cell type similarities/",tis,"_ranks.tsv"))
  ranks <- ranks[c(1:floor(length(rownames(ranks))/10)),]
  ranks <- ranks[which(ranks[,"value"] >= 0.05),]
  human_score <- read.delim(paste0("Output/PrEDiCT score/",tis,".tsv"))
  mouse_score <- read.delim(paste0("Output/Mouse/PrEDiCT score/",tis,".tsv"))
  ### Random ###
  repetitions <- 1000
  total_success <- c()
  for (r in c(1:repetitions)){
    run_success <- c()
    for (dis in unique(human_score$Disease.Phenotypic.series.accession..OMIM.)){
      human_pref_cts <- human_score[human_score$PrEDiCT.score >= 2 & human_score$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
      if (length(human_pref_cts) > 0){
        mouse_ortho_cts <- ranks[ranks$Human %in% human_pref_cts,"Mouse"]
        mouse_pref_cts <- mouse_score[mouse_score$PrEDiCT.score >= 2 & mouse_score$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
        mouse_random_cts <- sample(unique(mouse_score$Cell.type), length(mouse_pref_cts))
        run_success <- c(run_success, success_test(pref_cts = mouse_random_cts, simi_cts = mouse_ortho_cts))
      }
    }
    num_success <- sum(run_success)
    total_success <- c(total_success, num_success)
  }
  ### Real data ###
  run_success <- c()
  dis_num <- 0
  for (dis in unique(human_score$Disease.Phenotypic.series.accession..OMIM.)){
    human_pref_cts <- human_score[human_score$PrEDiCT.score >= 2 & human_score$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
    if (length(human_pref_cts) > 0){
      dis_num <- dis_num + 1
      mouse_ortho_cts <- ranks[ranks$Human %in% human_pref_cts,"Mouse"]
      mouse_pref_cts <- mouse_score[mouse_score$PrEDiCT.score >= 2 & mouse_score$Disease.Phenotypic.series.accession..OMIM. == dis,"Cell.type"]
      run_success <- c(run_success, success_test(pref_cts = mouse_pref_cts, simi_cts = mouse_ortho_cts))
    }
  }
  num_success <- sum(run_success)
  real_success <- num_success
  
  #### Check P value
  p <- sum(total_success >= real_success)/repetitions
  p_vals <- c(p_vals, p)
  print(tis)
}
print(FDR(pvalues = data.frame(tis = tissues,ps = p_vals)))
