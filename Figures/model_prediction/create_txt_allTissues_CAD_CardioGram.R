# combine all tissues results (CAD cardiogram)
library(qvalue)

setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/')
tissues_name <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage','Heart_Left_Ventricle', 'Liver', 'Whole_Blood')

df_tscore <- df_pathR <- df_pathGO <- df_pathWiki <- list()

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  
  tmp <- get(load(sprintf('predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/pval_Dx_pheno_covCorr.RData', t)))
  df_tscore[[i]] <- tmp$tscore[[1]]
  df_tscore[[i]]$tissue <-  t
  
  # add genes in the pathway and if there is an improvment in significance
  df_pathR[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR[[i]]$genes_path <- NA
  df_pathR[[i]]$improvement_sign <- NA
  for(j in 1:nrow(df_pathR[[i]])){
    df_pathR[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_reactome[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathR[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_reactome[[1]][[j]]$tscore[,8] > df_pathR[[i]][j,13])
  }
  df_pathR[[i]]$tissue <- t
  
  df_pathGO[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO[[i]]$genes_path <- NA
  df_pathGO[[i]]$improvement_sign <- NA
  for(j in 1:nrow(df_pathGO[[i]])){
    df_pathGO[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_GO[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathGO[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_GO[[1]][[j]]$tscore[,8] > df_pathGO[[i]][j,15])
  }
  df_pathGO[[i]]$tissue <- t
  
  tmp <- get(load(sprintf('predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/pval_Dx_pheno_covCorr_customPath_WikiPath2019Human.RData', t)))
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
df_tscore$CAD_HARD_BHcorr_overall <- p.adjust(df_tscore[, 8], method = 'BH')
df_pathR$CAD_HARD_BHcorr_overall <- p.adjust(df_pathR[, 13], method = 'BH')
df_pathGO$CAD_HARD_BHcorr_overall <- p.adjust(df_pathGO[, 15], method = 'BH')
df_pathWiki$CAD_HARD_BHcorr_overall <- p.adjust(df_pathWiki[, 13], method = 'BH')


# create a function to remove pathway with 1 gene and recompute pvalues
recompte_path <- function(tissues_name, res, id_pval){
  tmp <- lapply(tissues_name, function(x) res[res$tissue == x & res$ngenes_tscore>1,])
  for(i in 1:length(tmp)){
    tmp[[i]][, id_pval+1] <- qvalue(tmp[[i]][, id_pval])$qvalue
    tmp[[i]][, id_pval+2] <- p.adjust(tmp[[i]][, id_pval], method = 'BH')
  }
  tmp <- do.call(rbind, tmp)
  tmp[, id_pval+9] <- p.adjust(tmp[, id_pval], method = 'BH')
  return(tmp)
}
df_pathR_red <- recompte_path(res = df_pathR, tissues_name = tissues_name, id_pval = 13)
df_pathGO_red <- recompte_path(res = df_pathGO, tissues_name = tissues_name, id_pval = 15)
df_pathWiki_red <- recompte_path(res = df_pathWiki, tissues_name = tissues_name, id_pval = 13)


### save results
write.table(x = df_tscore, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/tscore_pval_Dx_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_Reactome_pval_Dx_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_GO_pval_Dx_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_WikiPath2019Human_pval_Dx_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_Reactome_pval_Dx_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_GO_pval_Dx_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_WikiPath2019Human_pval_Dx_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)







