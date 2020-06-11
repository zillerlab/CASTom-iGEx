# combine all tissues results (CAD)
library(qvalue)

setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/')
tissues_name <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage','Heart_Left_Ventricle', 'Liver')

df_tscore_HARD <- df_pathR_HARD <- df_pathGO_HARD <- list()
df_tscore_SOFT <- df_pathR_SOFT <- df_pathGO_SOFT <- list()
  
for(i in 1:length(tissues_name)){
  t <- tissues_name[i]
  tmp <- get(load(sprintf('predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData', t)))
  df_tscore_HARD[[i]] <- tmp$tscore[[1]]
  df_tscore_SOFT[[i]] <- tmp$tscore[[2]]
  df_tscore_HARD[[i]]$tissue <-  df_tscore_SOFT[[i]]$tissue <- t
  
  df_pathR_HARD[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR_SOFT[[i]] <- tmp$pathScore_reactome[[2]]
  df_pathR_HARD[[i]]$tissue <-  df_pathR_SOFT[[i]]$tissue <- t
  
  df_pathGO_HARD[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO_SOFT[[i]] <- tmp$pathScore_GO[[2]]
  df_pathGO_HARD[[i]]$tissue <-  df_pathGO_SOFT[[i]]$tissue <- t
}

df_tscore_HARD <- do.call(rbind, df_tscore_HARD)
df_tscore_SOFT <- do.call(rbind, df_tscore_SOFT)
df_pathR_HARD <- do.call(rbind, df_pathR_HARD)
df_pathR_SOFT <- do.call(rbind, df_pathR_SOFT)
df_pathGO_HARD <- do.call(rbind, df_pathGO_HARD)
df_pathGO_SOFT <- do.call(rbind, df_pathGO_SOFT)

# add overall corrected pval
df_tscore_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_tscore_HARD[, 8], method = 'BH')
df_tscore_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_tscore_SOFT[, 8], method = 'BH')
df_pathR_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_pathR_HARD[, 13], method = 'BH')
df_pathR_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_pathR_SOFT[, 13], method = 'BH')
df_pathGO_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_pathGO_HARD[, 15], method = 'BH')
df_pathGO_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_pathGO_SOFT[, 15], method = 'BH')


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
df_pathR_HARD_red <- recompte_path(res = df_pathR_HARD, tissues_name = tissues_name, id_pval = 13)
df_pathR_SOFT_red <- recompte_path(res = df_pathR_SOFT, tissues_name = tissues_name, id_pval = 13)
df_pathGO_HARD_red <- recompte_path(res = df_pathGO_HARD, tissues_name = tissues_name, id_pval = 15)
df_pathGO_SOFT_red <- recompte_path(res = df_pathGO_SOFT, tissues_name = tissues_name, id_pval = 15)


### save results
write.table(x = df_tscore_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/tscore_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_tscore_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/tscore_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_HARD_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_HARD_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_SOFT_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_SOFT_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_HARD_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_HARD_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_SOFT_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_SOFT_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)









