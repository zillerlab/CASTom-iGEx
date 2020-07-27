# combine all tissues results (Asthma)
library(qvalue)

setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
tissues_name <- c('DLPC_CMC', 'Whole_Blood')
pheno_name <- c('LifetimeMDD', 'MDDRecur')


for(p in pheno_name){
  
  df_tscore <- df_pathR <- df_pathGO <- list()
  
  for(i in 1:length(tissues_name)){
    
    t <- tissues_name[i]
    if(t == 'DLPC_CMC'){
      file <- sprintf('OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/%s_pheno/pval_%s_pheno_covCorr.RData', p,p)
    }else{
      file <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/%s_pheno/pval_%s_pheno_covCorr.RData',t, p, p)
    }
    
    tmp <- get(load(file))
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
  df_tscore <- cbind(df_tscore, p.adjust(df_tscore[, 8], method = 'BH'))
  df_pathR<- cbind(df_pathR, p.adjust(df_pathR[, 13], method = 'BH'))
  df_pathGO <- cbind(df_pathGO, p.adjust(df_pathGO[, 15], method = 'BH'))
  colnames(df_tscore)[ncol(df_tscore)] <- sprintf('%s_BHcorr_overall', p)
  colnames(df_pathR)[ncol(df_pathR)] <- sprintf('%s_BHcorr_overall', p)
  colnames(df_pathGO)[ncol(df_pathGO)] <- sprintf('%s_BHcorr_overall', p)
  
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
  write.table(x = df_tscore, file = sprintf('OUTPUT_all/%s_pheno/tscore_pval_%s_covCorr.txt', p, p), col.names=T, row.names=F, sep = '\t', quote = F)
  write.table(x = df_pathR, file = sprintf('OUTPUT_all/%s_pheno/path_Reactome_pval_%s_covCorr.txt', p, p), col.names=T, row.names=F, sep = '\t', quote = F)
  write.table(x = df_pathGO, file = sprintf('OUTPUT_all/%s_pheno/path_GO_pval_%s_covCorr.txt', p, p), col.names=T, row.names=F, sep = '\t', quote = F)
  
  write.table(x = df_pathR_red, file = sprintf('OUTPUT_all/%s_pheno/path_Reactome_pval_%s_covCorr_filt.txt', p, p), col.names=T, row.names=F, sep = '\t', quote = F)
  write.table(x = df_pathGO_red, file = sprintf('OUTPUT_all/%s_pheno/path_GO_pval_%s_covCorr_filt.txt', p, p), col.names=T, row.names=F, sep = '\t', quote = F)
  
}
