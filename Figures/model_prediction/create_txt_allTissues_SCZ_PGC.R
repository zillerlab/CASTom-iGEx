# combine all tissues results (SCZ)
library(qvalue)

setwd('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ')
tissues_name <- read.table('Tissues_PGC', h=F, stringsAsFactors = F)$V1

df_tscore <- df_pathR <- df_pathGO <- df_pathWiki <- list()

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr.RData', t)))
  df_tscore[[i]] <- tmp$tscore[[1]]
  df_tscore[[i]]$tissue <- t
  
  # add genes in the pathway and if there is an improvment in significance
  df_pathR[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR[[i]]$genes_path <- NA
  df_pathR[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathR[[i]])){
    df_pathR[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_reactome[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathR[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_reactome[[1]][[j]]$tscore[,8] > df_pathR[[i]][j,13])
  }
  df_pathR[[i]]$tissue <- t
  
  df_pathGO[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO[[i]]$genes_path <- NA
  df_pathGO[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathGO[[i]])){
    df_pathGO[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_GO[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathGO[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_GO[[1]][[j]]$tscore[,8] > df_pathGO[[i]][j,15])
  }
  df_pathGO[[i]]$tissue <- t
  
  
  tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr_customPath_WikiPath2019Human.RData', t)))
  df_pathWiki[[i]] <- tmp$pathScore[[1]]
  df_pathWiki[[i]]$genes_path <- NA
  df_pathWiki[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathWiki[[i]])){
    df_pathWiki[[i]]$genes_path[j] <- paste0(tmp$info_pathScore[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathWiki[[i]]$improvement_sign[j] <- all(tmp$info_pathScore[[1]][[j]]$tscore[,8] > df_pathWiki[[i]][j,13])
  }
  df_pathWiki[[i]]$tissue <- t
}

df_tscore <- do.call(rbind, df_tscore)
df_pathR <- do.call(rbind, df_pathR)
df_pathGO <- do.call(rbind, df_pathGO)
df_pathWiki <- do.call(rbind, df_pathWiki)


# add overall corrected pval
df_tscore$Dx_pval_BHcorr_overall <- p.adjust(df_tscore[, 8], method = 'BH')
df_pathR$Dx_pval_BHcorr_overall <- p.adjust(df_pathR[, 13], method = 'BH')
df_pathGO$Dx_pval_BHcorr_overall <- p.adjust(df_pathGO[, 15], method = 'BH')
df_pathWiki$Dx_pval_BHcorr_overall <- p.adjust(df_pathWiki[, 13], method = 'BH')


# create a function to remove pathway with 1 gene and recompute pvalues
recompte_path <- function(tissues_name, res, id_pval){
  tmp <- lapply(tissues_name, function(x) res[res$tissue == x & res$ngenes_tscore>1,])
  for(i in 1:length(tmp)){
    tmp[[i]][, id_pval+1] <- qvalue(tmp[[i]][, id_pval])$qvalue
    tmp[[i]][, id_pval+2] <- p.adjust(tmp[[i]][, id_pval], method = 'BH')
  }
  tmp <- do.call(rbind, tmp)
  tmp[, id_pval+9] <- p.adjust(tmp[, id_pval], method = 'BH')
  colnames(tmp)[ncol(tmp)] <- 'Dx_pval_BHcorr_overall'
  return(tmp)
}

df_pathR_red <- recompte_path(res = df_pathR, tissues_name = tissues_name, id_pval = 13)
df_pathGO_red <- recompte_path(res = df_pathGO, tissues_name = tissues_name, id_pval = 15)
df_pathWiki_red <- recompte_path(res = df_pathWiki, tissues_name = tissues_name, id_pval = 13)


### save results
write.table(x = df_tscore, file = 'OUTPUT_all/tscore_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR, file = 'OUTPUT_all/path_Reactome_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO, file = 'OUTPUT_all/path_GO_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki, file = 'OUTPUT_all/customPath_WikiPath2019Human_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_red, file = 'OUTPUT_all/path_Reactome_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_red, file = 'OUTPUT_all/path_GO_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki_red, file = 'OUTPUT_all/customPath_WikiPath2019Human_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)

# CMC pathways
t <- tissues_name[1]
tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr_customPath_CMC_GeneSets.RData', t)))
df_pathCMC <- tmp$pathScore[[1]]
df_pathCMC$genes_path <- NA
df_pathCMC$improvement_sign  <- NA
for(j in 1:nrow(df_pathCMC)){
  df_pathCMC$genes_path[j] <- paste0(tmp$info_pathScore[[1]][[j]]$tscore$external_gene_name, collapse = ',')
  df_pathCMC$improvement_sign[j] <- all(tmp$info_pathScore[[1]][[j]]$tscore[,8] > df_pathCMC[j,13])
}
df_pathCMC$tissue <- t

df_pathCMC_red <- recompte_path(res = df_pathCMC, tissues_name = t, id_pval = 13)
write.table(x = df_pathCMC, file = 'DLPC_CMC/customPath_CMC_GeneSets_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathCMC_red, file = 'DLPC_CMC/customPath_CMC_GeneSets_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)


