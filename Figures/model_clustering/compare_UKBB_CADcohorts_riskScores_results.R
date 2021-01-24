options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(apcluster))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')

# compare UKBB and meta-analysis results
setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/')

tissues <- c('Liver', 'Heart_Left_Ventricle')

for(tissue in tissues){
  
  meta_fold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
  ukbb_fold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue)
  
  METAriskScore_analysis_file <- sprintf('%sriskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis.RData', meta_fold)
  UKBBriskScore_analysis_file <- c(sprintf('%swithMedication_riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', ukbb_fold), 
                                   sprintf('%swithoutMedication_riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', ukbb_fold))
  
  UKBB_analysis_file <- c(sprintf('%srescaleCont_withMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', ukbb_fold), 
                                   sprintf('%srescaleCont_withoutMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', ukbb_fold))
  
  outFold <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB_CADSchukert/CAD_HARD_clustering/'
  color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
  
  # load ukbb
  original_res <- list()
  for(i in 1:length(UKBB_analysis_file)){
    
    tmp <- get(load(UKBB_analysis_file[i]))
    original_res[[i]] <- tmp$bin_reg
    
    if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
      tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
      tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
      tmp$phenoInfo$pheno_type <- tmp_name
      tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
    }
    original_res[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(original_res[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
  }
  original_res <- do.call(rbind, original_res)
  if(length(UKBB_analysis_file)>1){
    comp <- unique(original_res$comp)
    tmp <- list()
    for(i in 1:length(comp)){
      tmp[[i]] <- original_res[original_res$comp == comp[i],]
      tmp[[i]]$pval_corr <- p.adjust(tmp[[i]]$pvalue, method = 'BH')
    }
    original_res <- do.call(rbind, tmp)
  }
  
  # load ukbb
  ukbb_res <- list()
  for(i in 1:length(UKBBriskScore_analysis_file)){
    tmp <- get(load(UKBBriskScore_analysis_file[i]))
    ukbb_res[[i]] <- tmp$bin_reg
    
    if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
      tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
      tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
      tmp$phenoInfo$pheno_type <- tmp_name
      tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
    }
    ukbb_res[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(ukbb_res[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
  }
  ukbb_res <- do.call(rbind, ukbb_res)
  
  # load meta
  tmp <- get(load(METAriskScore_analysis_file))
  meta_res <- tmp$meta_analysis
  
  if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
    tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
    tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
    tmp$phenoInfo$pheno_type <- tmp_name
    tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
  }
  meta_res$pheno_type <- tmp$phenoInfo$pheno_type[match(meta_res$pheno_id,tmp$phenoInfo$pheno_id)]
  
  meta_res_tot <- meta_res
  ukbb_res_tot <- ukbb_res
  
  common_p <- intersect(unique(meta_res$pheno_id), unique(ukbb_res$pheno_id))
  meta_res <- meta_res[meta_res$pheno_id %in% common_p, ]
  ukbb_res <- ukbb_res[ukbb_res$pheno_id %in% common_p, ]
  
  comp <- unique(meta_res$comp)
  
  # for each group plot beta +/ CI
  df <- list()
  for(i in 1:length(comp)){
    
    tmp1 <- meta_res[meta_res$comp == comp[i], ]
    tmp2 <- ukbb_res[ukbb_res$comp == comp[i], ]
    tmp_pheno <- intersect(tmp1$pheno_id, tmp2$pheno_id)
    df[[i]] <- data.frame(beta_meta = tmp1$beta[match(tmp_pheno, tmp1$pheno_id)], beta_ukbb = tmp2$beta[match(tmp_pheno, tmp2$pheno_id)], 
                          CI_low_meta = tmp1$CI_low[match(tmp_pheno, tmp1$pheno_id)], CI_low_ukbb = tmp2$CI_low[match(tmp_pheno, tmp2$pheno_id)], 
                          CI_up_meta = tmp1$CI_up[match(tmp_pheno, tmp1$pheno_id)], CI_up_ukbb = tmp2$CI_up[match(tmp_pheno, tmp2$pheno_id)])
    df[[i]]$comp <- strsplit(comp[i], split = '_vs_all')[[1]][1]
    
  }
  
  df <- do.call(rbind, df)
  gr_name <- unique(df$comp)
  # plot
  df$comp <- factor(df$comp)
  
  df_ed <- data.frame(comp = gr_name, ed = sapply(gr_name, function(x) sqrt(sum((df$beta_meta[df$comp == x]- df$beta_ukbb[df$comp == x])^2))))
  df_ed$label = paste0('euclidian distance: ', as.character(round(df_ed$ed, digits = 2)))
  df_ed$comp <- factor(df_ed$comp)
  
  pl <- ggplot(data = df, aes(x = beta_ukbb, y = beta_meta))+
    geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2)+
    geom_point(alpha = 0.5, size = 1)+
    geom_errorbar(data = df, aes(x = beta_ukbb, y = beta_meta, ymin = CI_low_meta, ymax = CI_up_meta), alpha = 0.5) + 
    geom_errorbar(data = df, aes(x = beta_ukbb, y = beta_meta, xmin = CI_low_ukbb, xmax = CI_up_ukbb), alpha = 0.5) +
    facet_wrap(.~comp, ncol = length(comp), scales = 'free')+
    geom_text(data = df_ed, aes(x = -Inf, y = Inf, hjust = -0.01, vjust = 1.5, label = label), size = 3.5)+
    xlab('UKBB RS')+ ylab('meta-analysis German cohorts RS')+
    theme_bw()+ 
    theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
          legend.text = element_text(size = 6), legend.title = element_blank())
  ggsave(filename =  sprintf('%scl%s_compare_UKBB_meta_riskScore_endophenotype.png', outFold, tissue), plot = pl, width = 13, height = 3.5)
  ggsave(filename =  sprintf('%scl%s_compare_UKBB_meta_riskScore_endophenotype.pdf', outFold, tissue), plot = pl, width = 13, height = 3.5)
  
  
  #######
  original_res <- original_res[!is.na(original_res$beta),]
  original_res <- original_res[!is.infinite(original_res$CI_up) & !is.infinite(original_res$CI_low),]
  
  common_p <- intersect(unique(original_res$pheno_id), unique(ukbb_res_tot$pheno_id))
  original_res <- original_res[original_res$pheno_id %in% common_p, ]
  ukbb_res_tot <- ukbb_res_tot[ukbb_res_tot$pheno_id %in% common_p, ]
  
  comp <- unique(meta_res$comp)
  
  # for each group plot beta +/ CI
  df <- list()
  for(i in 1:length(comp)){
    
    tmp1 <- original_res[original_res$comp == comp[i], ]
    id_b <- tmp1$type_pheno != 'CONTINUOUS'
    tmp1$CI_low[id_b] <- log(tmp1$CI_low[id_b])
    tmp1$CI_up[id_b] <- log(tmp1$CI_up[id_b])
    
    tmp2 <- ukbb_res[ukbb_res$comp == comp[i], ]
    tmp_pheno <- intersect(tmp1$pheno_id, tmp2$pheno_id)
    df[[i]] <- data.frame(beta_o = tmp1$beta[match(tmp_pheno, tmp1$pheno_id)], beta_ukbb = tmp2$beta[match(tmp_pheno, tmp2$pheno_id)], 
                          CI_low_o = tmp1$CI_low[match(tmp_pheno, tmp1$pheno_id)], CI_low_ukbb = tmp2$CI_low[match(tmp_pheno, tmp2$pheno_id)], 
                          CI_up_o = tmp1$CI_up[match(tmp_pheno, tmp1$pheno_id)], CI_up_ukbb = tmp2$CI_up[match(tmp_pheno, tmp2$pheno_id)])
    df[[i]]$comp <- strsplit(comp[i], split = '_vs_all')[[1]][1]
    
  }
  
  df <- do.call(rbind, df)
  gr_name <- unique(df$comp)
  # plot
  df$comp <- factor(df$comp)
  
  df_ed <- data.frame(comp = gr_name, ed = sapply(gr_name, function(x) sqrt(sum((df$beta_o[df$comp == x] - df$beta_ukbb[df$comp == x])^2, na.rm = T))))
  df_ed$label = paste0('euclidian distance: ', as.character(round(df_ed$ed, digits = 2)))
  df_ed$comp <- factor(df_ed$comp)
  
  pl <- ggplot(data = df, aes(x = beta_ukbb, y = beta_o))+
    geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2)+
    geom_hline(yintercept = 0, color = 'red', linetype = 2)+
    geom_vline(xintercept = 0, color = 'red', linetype = 2)+
    geom_point(alpha = 0.5, size = 1)+
    geom_errorbar(data = df, aes(x = beta_ukbb, y = beta_o, ymin = CI_low_o, ymax = CI_up_o), alpha = 0.5) + 
    geom_errorbar(data = df, aes(x = beta_ukbb, y = beta_o, xmin = CI_low_ukbb, xmax = CI_up_ukbb), alpha = 0.5) +
    facet_wrap(.~comp, ncol = length(comp), scales = 'free')+
    geom_text(data = df_ed, aes(x = -Inf, y = Inf, hjust = -0.01, vjust = 1.5, label = label), size = 3.5)+
    xlab('UKBB RS')+ ylab('UKBB endophenotype')+
    theme_bw()+ 
    theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
          legend.text = element_text(size = 6), legend.title = element_blank())
  ggsave(filename =  sprintf('%scl%s_compare_UKBB_riskScore_endophenotype_original.png', outFold, tissue), plot = pl, width = 13, height = 3.5)
  ggsave(filename =  sprintf('%scl%s_compare_UKBB_riskScore_endophenotype_original.pdf', outFold, tissue), plot = pl, width = 13, height = 3.5)
  
  
}



