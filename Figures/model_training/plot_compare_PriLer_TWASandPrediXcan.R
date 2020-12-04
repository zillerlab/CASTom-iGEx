# plots comparison PriLer vs prediXcan and TWAS
# on denbi
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(latex2exp))

# each tissue can have a different model: noGWAS, CAD_GWAS, PGC_GWAS
tissues_model <- read.csv('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/final_model_gtex.csv', h=F, stringsAsFactors = F)
colnames(tissues_model) <- c('tissue', 'type')

res_CMC <- '/mnt/lucia/PriLer_TRAIN_PLOTS/res_CMC/'
res_GTEx <- '/mnt/lucia/PriLer_TRAIN_PLOTS/res_GTEx/'

df_tot_TWAS_v7 <- read.table(sprintf('%s/compare_PriLer_TWAS.txt', res_CMC), h=T, stringsAsFactors = F)
df_tot_TWAS_v7$tissue <- 'DLPC_CMC'
# remove some columns to have in common with gtex
df_tot_TWAS_v7 <- df_tot_TWAS_v7[, !colnames(df_tot_TWAS_v7) %in% c('top1.nsnps', 'top1.r2', 'top1.pv','bslmm.nsnps', 'bslmm.r2', 'bslmm.pv', 'id')]

for(t in tissues_model$tissue){
  
  print(t)
  tmp <- read.table(sprintf('%s/%s_compare_PriLer_TWAS_v7.txt', res_GTEx, t), h=T, stringsAsFactors = F)
  tmp$tissue <- t
  id <- sapply(colnames(df_tot_TWAS_v7), function(x) which(colnames(tmp) == x))
  tmp <- tmp[, id]
  df_tot_TWAS_v7 <- rbind(df_tot_TWAS_v7, tmp)
  
}

df_tot_TWAS_v6p <- read.table(sprintf('%s/compare_PriLer_TWAS.txt', res_CMC), h=T, stringsAsFactors = F)
df_tot_TWAS_v6p$tissue <- 'DLPC_CMC'
# remove some columns to have in common with gtex
df_tot_TWAS_v6p <- df_tot_TWAS_v6p[, !colnames(df_tot_TWAS_v6p) %in% c('top1.nsnps', 'top1.r2', 'top1.pv', 'id')]

for(t in tissues_model$tissue){
  
  print(t)
  tmp <- read.table(sprintf('%s/%s_compare_PriLer_TWAS_v6p.txt', res_GTEx, t), h=T, stringsAsFactors = F)
  tmp$tissue <- t
  id <- sapply(colnames(df_tot_TWAS_v6p), function(x) which(colnames(tmp) == x))
  tmp <- tmp[, id]
  df_tot_TWAS_v6p <- rbind(df_tot_TWAS_v6p, tmp)
  
}


df_tot_predixOld <- read.table(sprintf('%s/compare_PriLer_prediXcan.txt', res_CMC), h=T, stringsAsFactors = F)
df_tot_predixOld <- df_tot_predixOld[,!colnames(df_tot_predixOld) %in% 'genename']
df_tot_predixOld$tissue <- 'DLPC_CMC'

for(t in tissues_model$tissue){
  
  print(t)
  tmp <- read.table(sprintf('%s/%s_compare_PriLer_prediXcan_v6p.txt', res_GTEx, t), h=T, stringsAsFactors = F)
  tmp$tissue <- t
  id <- sapply(colnames(df_tot_predixOld), function(x) which(colnames(tmp) == x))
  tmp <- tmp[, id]
  df_tot_predixOld <- rbind(df_tot_predixOld, tmp)
  
}


df_tot_predixNew <- NULL

for(t in tissues_model$tissue){
  
  print(t)
  tmp <- read.table(sprintf('%s/%s_compare_PriLer_prediXcan_v7.txt', res_GTEx, t), h=T, stringsAsFactors = F)
  tmp$tissue <- t
  df_tot_predixNew <- rbind(df_tot_predixNew, tmp)
  
}
# load number of samples for each method 
CMC_sampleSize <- read.table('/mnt/lucia/PriLer_TRAIN_PLOTS/nSamples_methods.txt', h=F, stringsAsFactors = F)
PriLer_sampleSize <- read.table('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/Covariates/n_samples_allTissues.txt', h=F, stringsAsFactors = F)
PriLer_sampleSize <- rbind(data.frame(V1='DLPC_CMC', V2 = CMC_sampleSize$V2[CMC_sampleSize$V1 == 'PriLer']), PriLer_sampleSize)

TWAS_v6p_sampleSize <- read.table('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/TWAS/GTEx_v6p/n_samples_allTissues.txt', h=F, stringsAsFactors = F,  sep = '\t')
TWAS_v6p_sampleSize <-  rbind(data.frame(V1='DLPC_CMC', V2 = CMC_sampleSize$V2[CMC_sampleSize$V1 == 'TWAS']), TWAS_v6p_sampleSize)
TWAS_v7_sampleSize <- read.table('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/TWAS/GTEx_v7/n_samples_allTissues.txt', h=F, stringsAsFactors = F,  sep = '\t')
TWAS_v7_sampleSize <-  rbind(data.frame(V1='DLPC_CMC', V2 = CMC_sampleSize$V2[CMC_sampleSize$V1 == 'TWAS']), TWAS_v7_sampleSize)

prediXcanOld_sampleSize <- read.table('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/prediXcan/GTEx_v6p/n_samples_allTissues.txt', h=F, stringsAsFactors = F)
prediXcanOld_sampleSize <-  rbind(data.frame(V1='DLPC_CMC', V2 = CMC_sampleSize$V2[CMC_sampleSize$V1 == 'prediXcan']), prediXcanOld_sampleSize)

prediXcanNew_sampleSize <- read.table('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/prediXcan/GTEx_v7/n_samples_allTissues.txt', h=F, stringsAsFactors = F)

######################################################################################################
### function to produce plot and save summary (generalize for each dataset)
order_tissues <- rev(c(1:7, 22:23, 8:15, 19:21, 32, 30, 17:18,16, 24:29, 31, 33:34))
color_tissues <- rev(c('#000080', rep('#999900',2),'#3CB371', rep('#CD5C5C', 3), rep('#DC143C', 2), rep('#4169E1', 8), rep('#8B4513', 7), '#FF8C00', '#8A2BE2','#708090','#A52A2A',
                       '#3CB371', rep('#FA8072', 2), rep('#3CB371', 2), '#FF8C00'))
color_other_TWAS_v6p <- "#FF7F0EFF"
color_other_TWAS_v7 <- "#FF7F50"
color_other_prediXcan_v6p <- "#D62728FF"
color_other_prediXcan_v7 <-  "#9467BDFF"

plots_summary_comp <- function(df_tot, PriLer_eval, other_eval, other_nsnps, other_nsnps_window = NA, subset_tissue = NA, 
                               other_name, order_tissues, color_tissues, other_name_save, eval_title_plot, color_other, sampleSize_file_other, height_val = 10, width_val = 6.2){
  
  # remove test NA
  df_tot_noNA <- df_tot[!is.na(df_tot[, PriLer_eval]) & !is.na(df_tot$n_snps),]
  if(!is.na(subset_tissue)){
    df_tot_noNA <- df_tot_noNA[df_tot_noNA$tissue %in% subset_tissue, ]
  }
  df_pl <- data.frame(corr2 = c(df_tot_noNA[, PriLer_eval]^2, df_tot_noNA[, other_eval]), nsnps = c(df_tot_noNA$n_snps, df_tot_noNA[,other_nsnps]))
  if(!is.na(other_nsnps_window)){
    df_pl$nsnps_window <- c(df_tot_noNA$n_snps_gene_window, df_tot_noNA[,other_nsnps_window])
  }
  df_pl$tissue <- rep(df_tot_noNA$tissue,2)
  df_pl$method <-  c(rep('PriLer',nrow(df_tot_noNA)), rep(other_name, nrow(df_tot_noNA)))
  
  df_pl$log_nspns <- log2(df_pl$nsnps+1) 
  if(!is.na(other_nsnps_window)){df_pl$log_nspns_window <- log2(df_pl$nsnps_window+1)}
  names_t <- unique(df_pl$tissue)
  df_pl$tissue <- factor(df_pl$tissue, levels = names_t[order_tissues])
  df_pl$method <- factor(df_pl$method, levels = rev(c('PriLer', other_name)))
  
  # test for different distributions (corr2)
  test_res <- lapply(levels(df_pl$tissue), function(x) wilcox.test(df_pl$corr2[df_pl$tissue == x & df_pl$method == 'PriLer'],
                                                                   df_pl$corr2[df_pl$tissue == x & df_pl$method == other_name], paired = T))
  
  sign_diff_method <- sapply(test_res, function(x) x$p.value)
  # diff_method <-  sapply(test_res, function(x) x$estimate)
  diff_method <- sapply(levels(df_pl$tissue), function(x) median(df_pl$corr2[df_pl$tissue == x & df_pl$method == 'PriLer']-
                                                                   df_pl$corr2[df_pl$tissue == x & df_pl$method == other_name]))
  
  df_sign <- data.frame(sign = sign_diff_method, est = diff_method, tissue = levels(df_pl$tissue))
  df_sign$new_names <- df_sign$tissue
  df_sign$new_names[df_sign$sign<=0.01 & df_sign$sign>0.001] <- paste0('*', df_sign$new_names[df_sign$sign<=0.01 & df_sign$sign>0.001])
  df_sign$new_names[df_sign$sign<=0.001 & df_sign$sign>0.00001]<- paste0('**', df_sign$new_names[df_sign$sign<=0.001 & df_sign$sign>0.00001])
  df_sign$new_names[df_sign$sign<=0.00001]<- paste0('***', df_sign$new_names[df_sign$sign<=0.00001])
  df_sign$face <- 'plain'
  df_sign$face[df_sign$est>0] <- 'bold.italic'
  df_sign$color <- color_tissues
  
  df_diff <- data.frame(tissue = df_sign$tissue, w.test_corr2_pval = df_sign$sign, w.test_corr2_median = df_sign$est, cor_corr2_pval = 0, cor_corr2_est = 0, n_genes = 0)
  df_diff$n_genes <- sapply(df_diff$tissue, function(x) length(which(df_tot_noNA$tissue == x)))
  df_diff$cor_corr2_pval <- sapply(df_diff$tissue, function(x) cor.test(df_tot_noNA[df_tot_noNA$tissue == x, PriLer_eval]^2, df_tot_noNA[df_tot_noNA$tissue == x, other_eval])$p.value)
  df_diff$cor_corr2_est <- sapply(df_diff$tissue, function(x) cor.test(df_tot_noNA[df_tot_noNA$tissue == x, PriLer_eval]^2, df_tot_noNA[df_tot_noNA$tissue == x, other_eval])$est)
  
  df_diff <- rbind(df_diff, data.frame(tissue = 'All', 
                                       w.test_corr2_pval =  wilcox.test(df_pl$corr2[df_pl$method == 'PriLer'],df_pl$corr2[df_pl$method == other_name], paired = T)$p.value, 
                                       w.test_corr2_median =  median(df_pl$corr2[df_pl$method == 'PriLer'] - df_pl$corr2[df_pl$method == other_name]), 
                                       cor_corr2_pval = cor.test(df_tot_noNA[, PriLer_eval]^2, df_tot_noNA[, other_eval])$p.value, 
                                       cor_corr2_est = cor.test(df_tot_noNA[, PriLer_eval]^2, df_tot_noNA[, other_eval])$est, 
                                       n_genes = nrow(df_tot_noNA)))
  # save
  write.table(file = sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/Corr_pred_perTissue_PriLer_%s.txt', other_name_save), x = df_diff, col.names = T, row.names = F, sep = '\t', quote = F)
  
  pl_box <- ggplot(df_pl, aes(x = tissue, y = corr2, fill = method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.6) + ggtitle(eval_title_plot)+
    stat_summary(  position = position_dodge(width = 0.75), fun.y=mean, colour="darkred", geom="point", shape=18, size=1,show.legend = FALSE)+
    xlab('') +  theme_classic() + ylab(TeX('$corr^2$')) + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
                                                                axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
                                                                axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
                                                                axis.text.y = element_text(color = df_sign$color, face = df_sign$face))+
    scale_x_discrete(labels = df_sign$new_names)+
    scale_fill_manual(values = rev(c("#1F77B4FF",color_other)))+ coord_flip()
  
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_corr2_PriLer_%s.png', other_name_save)
  ggsave(file_name, pl_box, width = width_val, height = height_val, units = "in", dpi=500)
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_corr2_PriLer_%s.pdf', other_name_save)
  ggsave(file_name, pl_box, width = width_val, height = height_val, units = "in", dpi=500)
  
  # test for different distributions (nsnps)
  test_res_nsnps <- lapply(levels(df_pl$tissue), function(x) wilcox.test(df_pl$log_nspns[df_pl$tissue == x & df_pl$method == 'PriLer'],
                                                                         df_pl$log_nspns[df_pl$tissue == x & df_pl$method == other_name], paired = T))
  
  sign_diff_method_nsnps <- sapply(test_res_nsnps, function(x) x$p.value)
  diff_method_nsnps <-  sapply(levels(df_pl$tissue), function(x) median(df_pl$log_nspns[df_pl$tissue == x & df_pl$method == 'PriLer']-
                                                                          df_pl$log_nspns[df_pl$tissue == x & df_pl$method == other_name]))
  df_sign_nsnps <- data.frame(sign = sign_diff_method_nsnps, est = diff_method_nsnps, tissue = levels(df_pl$tissue))
  df_sign_nsnps$new_names <- df_sign_nsnps$tissue
  df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.01 & df_sign_nsnps$sign>0.001] <- paste0('*', df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.01 & df_sign_nsnps$sign>0.001])
  df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.001 & df_sign_nsnps$sign>0.00001]<- paste0('**', df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.001 & df_sign_nsnps$sign>0.00001])
  df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.00001]<- paste0('***', df_sign_nsnps$new_names[df_sign_nsnps$sign<=0.00001])
  df_sign_nsnps$face <- 'plain'
  df_sign_nsnps$face[df_sign_nsnps$est<0] <- 'bold.italic'
  df_sign_nsnps$color <- color_tissues
  
  pl_box_nsnps <- ggplot(df_pl, aes(x = tissue, y = log_nspns, fill = method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.6) + ggtitle(TeX('N. reg-SNPs final model'))+
    stat_summary(  position = position_dodge(width = 0.75), fun.y=mean, colour="darkred", geom="point", shape=18, size=1,show.legend = FALSE)+
    xlab('') +  theme_classic() + ylab(TeX('log2(n.snps+1)')) + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
                                                                      axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
                                                                      axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
                                                                      axis.text.y = element_text(color = df_sign_nsnps$color, face = df_sign_nsnps$face))+
    scale_x_discrete(labels = df_sign_nsnps$new_names)+
    scale_fill_manual(values = rev(c("#1F77B4FF",color_other)))+ coord_flip()
  
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsnps_PriLer_%s.png', other_name_save)
  ggsave(file_name, pl_box_nsnps, width = width_val, height = height_val, units = "in", dpi=500)
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsnps_PriLer_%s.pdf', other_name_save)
  ggsave(file_name, pl_box_nsnps, width = width_val, height = height_val, units = "in", dpi=500)
  
  if(!is.na(other_nsnps_window)){
    test_res_wnsnps <- lapply(levels(df_pl$tissue), function(x) wilcox.test(df_pl$log_nspns_window[df_pl$tissue == x & df_pl$method == 'PriLer'],
                                                                            df_pl$log_nspns_window[df_pl$tissue == x & df_pl$method == other_name], paired = T))
    
    sign_diff_method_wnsnps <- sapply(test_res_wnsnps, function(x) x$p.value)
    diff_method_wnsnps <-  sapply(levels(df_pl$tissue), function(x) median(df_pl$log_nspns_window[df_pl$tissue == x & df_pl$method == 'PriLer']-
                                                                             df_pl$log_nspns_window[df_pl$tissue == x & df_pl$method == other_name]))
    
    df_sign_wnsnps <- data.frame(sign = sign_diff_method_wnsnps, est = diff_method_wnsnps, tissue = levels(df_pl$tissue))
    df_sign_wnsnps$new_names <- df_sign_wnsnps$tissue
    df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.01 & df_sign_wnsnps$sign>0.001] <- paste0('*', df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.01 & df_sign_wnsnps$sign>0.001])
    df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.001 & df_sign_wnsnps$sign>0.00001]<- paste0('**', df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.001 & df_sign_wnsnps$sign>0.00001])
    df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.00001]<- paste0('***', df_sign_wnsnps$new_names[df_sign_wnsnps$sign<=0.00001])
    df_sign_wnsnps$face <- 'plain'
    df_sign_wnsnps$face[df_sign_wnsnps$est<0] <- 'bold.italic'
    df_sign_wnsnps$color <- color_tissues
    
    pl_box_wnsnps <- ggplot(df_pl, aes(x = tissue, y = log_nspns_window, fill = method)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.6) + ggtitle(TeX('N. cis-SNPs in gene window'))+
      stat_summary(  position = position_dodge(width = 0.75), fun.y=mean, colour="darkred", geom="point", shape=18, size=1,show.legend = FALSE)+
      xlab('') +  theme_classic() + ylab(TeX('log2(n.snps+1)')) + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
                                                                        axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
                                                                        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
                                                                        axis.text.y = element_text(color = df_sign_wnsnps$color, face = df_sign_wnsnps$face))+
      scale_x_discrete(labels = df_sign_wnsnps$new_names)+
      scale_fill_manual(values = rev(c("#1F77B4FF",color_other)))+ coord_flip()
    
    file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsnps_window_PriLer_%s.png', other_name_save)
    ggsave(file_name, pl_box_wnsnps, width = width_val, height = height_val, units = "in", dpi=500)
    file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsnps_window_PriLer_%s.pdf', other_name_save)
    ggsave(file_name, pl_box_wnsnps, width = width_val, height = height_val, units = "in", dpi=500)
    
  }
  
  #########################
  ### number of samples ###
  #########################
  common_tissues <- intersect(PriLer_sampleSize$V1, sampleSize_file_other$V1)
  tmp_PriLer <- PriLer_sampleSize[PriLer_sampleSize$V1 %in% common_tissues, ]
  tmp_other <- sampleSize_file_other[sampleSize_file_other$V1 %in% common_tissues, ]
  df_sample <- data.frame(tissue = c( tmp_PriLer$V1, tmp_other$V1), 
                          method = c(rep('PriLer', nrow(tmp_PriLer)), rep(other_name, nrow(tmp_other))),
                          n_samples = c(tmp_PriLer$V2, tmp_other$V2))
  if(!is.na(subset_tissue)){df_sample <- df_sample[df_sample$tissue %in% subset_tissue,]}
  df_sample$tissue <- factor(df_sample$tissue, levels =  levels(df_pl$tissue))
  df_sample$method <- factor(df_sample$method, levels = levels(df_pl$method))
  
  pl_bar <- ggplot(df_sample, aes(x = tissue, y = n_samples, fill = method)) +
    geom_bar(stat="identity", position=position_dodge(), alpha = 0.6) + ggtitle(TeX('Sample Size'))+
    xlab('') +  theme_classic() + ylab(TeX('n. samples')) + theme(legend.position = 'none', legend.direction  = 'vertical', plot.title = element_text(hjust = 0.5), 
                                                                  axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
                                                                  axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
                                                                  axis.text.y = element_text(color = color_tissues))+
    scale_fill_manual(values = rev(c("#1F77B4FF", color_other)))+ labs(color = "")+ coord_flip()
  
  
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsamples_PriLer_%s.png', other_name_save)
  ggsave(file_name, pl_bar, width = width_val, height = height_val, units = "in", dpi=500)
  file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_nsamples_PriLer_%s.pdf', other_name_save)
  ggsave(file_name, pl_bar, width = width_val, height = height_val, units = "in", dpi=500)
  
  ######################
  ## overall increase ##
  ######################
  if(is.na(subset_tissue)){
    df_tab <- data.frame(method = rep(c('PriLer', other_name),2), n_genes = c(table(df_pl$corr2[df_pl$method == 'PriLer'] >= df_pl$corr2[df_pl$method == other_name]),
                                                                              table(df_pl$nsnps[df_pl$method == 'PriLer'] >= df_pl$nsnps[df_pl$method == other_name])),
                         type = c(rep('CV',2), rep('n.reg-SNPs',2)))
    df_tab$type <- factor(df_tab$type)
    group_name <- c(TeX('CV $corr^2$'),'n.reg-SNPs')
    df_tab$method <- factor(df_tab$method, levels = c('PriLer', other_name))
    
    pl_bar <- ggplot(df_tab, aes(x = type, y = n_genes, fill = method)) +
      geom_bar(stat="identity", alpha = 0.6, width = 0.7) + ggtitle('34 tissues combined')+
      xlab('') +  theme_classic() + ylab(TeX('n. genes')) + theme(legend.position = 'bottom', legend.title = element_blank(), legend.direction  = 'vertical', plot.title = element_text(hjust = 0.5), 
                                                                  axis.text = element_text(size = 12), axis.title = element_text(size = 10), 
                                                                  axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))+
      scale_fill_manual(labels = parse(text = c(TeX(paste0('PriLer $<$', other_name)), TeX(paste0('PriLer $\\geq$', other_name)))), values = rev(c("#1F77B4FF", color_other)))+
      scale_x_discrete(labels=parse(text = group_name))+ coord_flip()
    
    file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_increaseAllTissues_PriLer_%s.png', other_name_save)
    ggsave(file_name, pl_bar, width = 5, height = 3.5, units = "in", dpi=500)
    file_name <- sprintf('/mnt/lucia/PriLer_TRAIN_PLOTS/comparison_methods_increaseAllTissues_PriLer_%s.pdf', other_name_save)
    ggsave(file_name, pl_bar, width = 5, height = 3.5, units = "in", dpi=500)
  }
}


#### TWAS (GTEx v6p)
plots_summary_comp(df_tot = df_tot_TWAS_v6p, PriLer_eval = 'test_comb_cor', other_eval = 'best.r2', other_nsnps = 'best.nsnps', other_nsnps_window = 'nsnps', 
                   other_name = 'TWAS (GTEx v6p)', other_name_save = 'TWAS_v6p', order_tissues = order_tissues, color_tissues = color_tissues, 
                   eval_title_plot =  TeX('CV $corr^2$ on test folds combined'), color_other = color_other_TWAS_v6p, sampleSize_file_other = TWAS_v6p_sampleSize)
# subset tissues:
subset_tissues <- c('DLPC_CMC', 'Adipose_Subcutaneous', 'Heart_Left_Ventricle', 'Whole_Blood')
plots_summary_comp(df_tot = df_tot_TWAS_v6p, subset_tissue = subset_tissues, PriLer_eval = 'test_comb_cor', other_eval = 'best.r2', other_nsnps = 'best.nsnps', other_nsnps_window = 'nsnps', 
                   other_name = 'TWAS (GTEx v6p)', other_name_save = 'TWAS_v6p_subsetTissues', order_tissues = rev(1:4), color_tissues = color_tissues[c(1, 26, 33,34)], 
                   eval_title_plot =  TeX('CV $corr^2$ on test folds combined'), color_other = color_other_TWAS_v6p, sampleSize_file_other = TWAS_v6p_sampleSize, width_val = 5, height_val = 3.5)

#### TWAS (GTEx v7)
plots_summary_comp(df_tot = df_tot_TWAS_v7, PriLer_eval = 'test_comb_cor', other_eval = 'best.r2', other_nsnps = 'best.nsnps', other_nsnps_window = 'nsnps', 
                   other_name = 'TWAS (GTEx v7)', other_name_save = 'TWAS_v7', order_tissues = order_tissues, color_tissues = color_tissues, 
                   eval_title_plot =  TeX('CV $corr^2$ on test folds combined'), color_other = color_other_TWAS_v7, sampleSize_file_other = TWAS_v7_sampleSize)

#### prediXcan (GTEx v6p)
plots_summary_comp(df_tot = df_tot_predixOld, PriLer_eval = 'test_comb_cor', other_eval = 'pred.perf.R2', other_nsnps = 'n.snps.in.model',
                   other_name = 'prediXcan (GTEx v6p)', other_name_save = 'prediXcan_v6p', order_tissues = order_tissues, color_tissues = color_tissues, 
                   eval_title_plot =  TeX('CV $corr^2$ on test folds combined'), color_other = color_other_prediXcan_v6p, sampleSize_file_other = prediXcanOld_sampleSize)
# subset tissues:
subset_tissues <- c('DLPC_CMC', 'Adipose_Subcutaneous', 'Heart_Left_Ventricle', 'Whole_Blood')
plots_summary_comp(df_tot = df_tot_predixOld, subset_tissue = subset_tissues, PriLer_eval = 'test_comb_cor', other_eval = 'pred.perf.R2', other_nsnps = 'n.snps.in.model', 
                   other_name = 'prediXcan (GTEx v6p)', other_name_save = 'prediXcan_v6p_subsetTissues', order_tissues = rev(1:4), color_tissues = color_tissues[c(1, 26, 33,34)], 
                   eval_title_plot =  TeX('CV $corr^2$ on test folds combined'), color_other = color_other_prediXcan_v6p, sampleSize_file_other = prediXcanOld_sampleSize, width_val = 5, height_val = 3.5)

#### prediXcan (GTEx v7)
plots_summary_comp(df_tot = df_tot_predixNew, PriLer_eval = 'test_cor', other_eval = 'pred.perf.R2', other_nsnps = 'n.snps.in.model', other_nsnps_window = 'n_snps_in_window',
                   other_name = 'prediXcan (GTEx v7)', other_name_save = 'prediXcan_v7', order_tissues = order_tissues[-length(order_tissues)]-1, color_tissues = color_tissues[-length(color_tissues)], 
                   eval_title_plot =  TeX('mean CV $corr^2$ on test folds'), color_other = color_other_prediXcan_v7, sampleSize_file_other = prediXcanNew_sampleSize)


