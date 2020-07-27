# combine all tissues results (Asthma)
library(qvalue)

setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/')
tissues_name <- c('Lung', 'Whole_Blood')

df_tscore <- df_pathR <- df_pathGO <- list()

for(i in 1:length(tissues_name)){
  t <- tissues_name[i]
  tmp <- get(load(sprintf('predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/pval_Asthma_pheno_covCorr.RData', t)))
  df_tscore[[i]] <- tmp$tscore[[1]]
  df_tscore[[i]]$tissue <-  t
  
  df_pathR[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR[[i]]$tissue  <- t
  
  df_pathGO[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO[[i]]$tissue <- t
  
}

df_tscore <- do.call(rbind, df_tscore)
df_pathR <- do.call(rbind, df_pathR)
df_pathGO <- do.call(rbind, df_pathGO)

# add overall corrected pval
df_tscore$Asthma_BHcorr_overall <- p.adjust(df_tscore[, 8], method = 'BH')
df_pathR$Asthma_BHcorr_overall <- p.adjust(df_pathR[, 13], method = 'BH')
df_pathGO$Asthma_BHcorr_overall <- p.adjust(df_pathGO[, 15], method = 'BH')


# create a function to remove pathway with 1 gene and recompute pvalues
recompte_path <- function(tissues_name, res, id_pval){
  tmp <- lapply(tissues_name, function(x) res[res$tissue == x & res$ngenes_tscore>1,])
  for(i in 1:length(tmp)){
    tmp[[i]][, id_pval+1] <- qvalue(tmp[[i]][, id_pval])$qvalue
    tmp[[i]][, id_pval+2] <- p.adjust(tmp[[i]][, id_pval], method = 'BH')
  }
  tmp <- do.call(rbind, tmp)
  tmp[, id_pval+4] <- p.adjust(tmp[, id_pval], method = 'BH')
  
  return(tmp)
}

df_pathR_red <- recompte_path(res = df_pathR, tissues_name = tissues_name, id_pval = 13)
df_pathGO_red <- recompte_path(res = df_pathGO, tissues_name = tissues_name, id_pval = 15)


### save results
write.table(x = df_tscore, file = 'predict_UKBB/AllTissues/200kb/noGWAS/Asthma_pheno/tscore_pval_Asthma_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR, file = 'predict_UKBB/AllTissues/200kb/noGWAS/Asthma_pheno/path_Reactome_pval_Asthma_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO, file = 'predict_UKBB/AllTissues/200kb/noGWAS/Asthma_pheno/path_GO_pval_Asthma_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_red, file = 'predict_UKBB/AllTissues/200kb/noGWAS/Asthma_pheno/path_Reactome_pval_Asthma_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_red, file = 'predict_UKBB/AllTissues/200kb/noGWAS/Asthma_pheno/path_GO_pval_Asthma_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)




