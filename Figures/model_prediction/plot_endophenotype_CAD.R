# plot OR/beta for endophenotypes association

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

# plot OR and beta regression coeff, for some of them plot violin plots

phenoFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/'
pheno_ann <- read.delim(sprintf('%scolor_pheno_type_UKBB.txt', phenoFold), header = T, stringsAsFactors = F)
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4'), pheno_type = c('ICD9-10_OPCS4', 'Medications')))
pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'

#######################################################################################################
tissue <- 'Liver'
fold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
withMed_glmFile <- sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)

foldPred <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
glmPredFile <- sprintf('%stscore_zscaled_clusterCases_phenoAssociationGLM_prediction_modelUKBB.RData', foldPred)


res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df$OR_or_beta <- NA
df$CI_low <- NA
df$CI_up <- NA
for(i in 1:nrow(df)){
  
  if(df$type_pheno[i] == 'CONTINUOUS'){
    df$OR_or_beta[i] <- df$beta[i]
    df$CI_low[i] <- df$OR_or_beta[i] + qnorm(0.025)*df$se_beta[i]
    df$CI_up[i] <- df$OR_or_beta[i] + qnorm(0.975)*df$se_beta[i]
  }else{
    df$OR_or_beta[i] <- exp(df$beta[i])
    df$CI_low[i] <- exp(df$beta[i] + qnorm(0.025)*df$se_beta[i])
    df$CI_up[i] <- exp(df$beta[i] + qnorm(0.975)*df$se_beta[i])
  }
}

# make plots
id_keep <- unique(df$pheno_id[df$pvalue <= 0.001])
df_red <- df[df$pheno_id %in% id_keep, ]

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1)+
  facet_wrap(comp~.,  nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red1$color)+
  scale_y_continuous(trans='log2')+
  coord_flip()

pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0)+
  facet_wrap(comp~., nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red2$color)+
  coord_flip()
len_h <- length(unique(df_red$pheno_id[df_red$type_pheno == 'CONTINUOUS'])) + length(unique(df_red$pheno_id[df_red$type_pheno != 'CONTINUOUS']))

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.3))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%sstscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'pdf')

### with Medication
res <- get(load(withMed_glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df$OR_or_beta <- NA
df$CI_low <- NA
df$CI_up <- NA
for(i in 1:nrow(df)){
  
  if(df$type_pheno[i] == 'CONTINUOUS'){
    df$OR_or_beta[i] <- df$beta[i]
    df$CI_low[i] <- df$OR_or_beta[i] + qnorm(0.025)*df$se_beta[i]
    df$CI_up[i] <- df$OR_or_beta[i] + qnorm(0.975)*df$se_beta[i]
  }else{
    df$OR_or_beta[i] <- exp(df$beta[i])
    df$CI_low[i] <- exp(df$beta[i] + qnorm(0.025)*df$se_beta[i])
    df$CI_up[i] <- exp(df$beta[i] + qnorm(0.975)*df$se_beta[i])
  }
}

# make plots
id_keep <- unique(df$pheno_id[df$pvalue <= 0.001])
df_red <- df[df$pheno_id %in% id_keep, ]

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1)+
  facet_wrap(comp~.,  nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red1$color)+
  scale_y_continuous(trans='log2')+
  coord_flip()

pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0)+
  facet_wrap(comp~., nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red2$color)+
  coord_flip()
len_h <- length(unique(df_red$pheno_id[df_red$type_pheno == 'CONTINUOUS'])) + length(unique(df_red$pheno_id[df_red$type_pheno != 'CONTINUOUS']))

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.3))

ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'pdf')


### plot specific phenotypes ###
# LDL and HDL
library(tibble)
library(scales)

pheno_id <- c('30780', '30760')
for(i in 1:length(pheno_id)){

  tmp_pheno <- data.frame(value = res$phenoDat[, pheno_id[i]], gr = paste0('gr',res$cl$gr))
  tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res$cl$gr))))
  stat_tmp <- res$bin_reg[res$bin_reg$pheno_id == pheno_id[i], colnames(res$bin_reg) %in% c('comp', 'pheno_id', 'Field', 'pvalue')]
  stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  stat_tmp <- stat_tmp[grepl('gr1',stat_tmp$group1) | grepl('gr1',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
  stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))

  stat_tmp$y.position <- mean(sapply(sort(unique(res$cl$gr)), function(x) max(tmp_pheno$value[tmp_pheno$gr == paste0('gr',x)], na.rm = T))) + 
  c(0.7, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7)
  stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
  stat_tmp <- as_tibble(stat_tmp)        
  title_name <- unique(df$Field[df$pheno_id == pheno_id[i]])
  
  pl <- ggviolin(tmp_pheno, x = "gr", y = "value", fill = "gr", alpha = 0.7, add = "boxplot", add.params = list(fill = "white"))+ 
    stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)

  pl <- ggpar(pl, legend = "none", xlab = '', ylab = title_name)
  file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise', fold)
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')

}

# load results from prediction:
res_pred <- get(load(glmPredFile))
pheno_id <- c('LDL', 'HDL')
cohort_name <- paste0('German', 1:5)

for(id in c(1,5)){
  for(i in 1:length(pheno_id)){
    
    tmp_pheno <- data.frame(value = res_pred$phenoDat[[id]][, pheno_id[i]], gr = paste0('gr',res_pred$cl[[id]]$gr))
    tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res$cl$gr))))
    stat_tmp <- res_pred$bin_reg[[id]][res_pred$bin_reg[[id]]$pheno_id == pheno_id[i], colnames(res_pred$bin_reg[[id]]) %in% c('comp', 'pheno_id', 'pvalue')]
    stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
    stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
    stat_tmp <- stat_tmp[grepl('gr1',stat_tmp$group1) | grepl('gr1',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
    stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))
    
    if(pheno_id[i] == 'LDL'){
      vect_add <- 11+16*(0:(nrow(stat_tmp)-1))
    }
    if(pheno_id[i] == 'HDL'){
      vect_add <- 7+7*(0:(nrow(stat_tmp)-1))
    }
    
    stat_tmp$y.position <- mean(sapply(sort(unique(res_pred$cl[[id]]$gr)), function(x) max(tmp_pheno$value[tmp_pheno$gr == paste0('gr',x)], na.rm = T))) + 
      vect_add
    stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
    stat_tmp <- as_tibble(stat_tmp)        
    title_name <- pheno_id[i]
    
    pl <- ggviolin(tmp_pheno, x = "gr", y = "value", fill = "gr", alpha = 0.7, add = "boxplot", add.params = list(fill = "white"))+ 
      stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)
    
    pl <- ggpar(pl, legend = "none", xlab = '', ylab = title_name)
    file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_prediction_modelUKBB_%s', foldPred, cohort_name[id])
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')
    
  }
}



#######################################################################################################
tissue <- 'Adipose_Visceral_Omentum'
fold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
foldPred <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
glmPredFile <- sprintf('%stscore_zscaled_clusterCases_phenoAssociationGLM_prediction_modelUKBB.RData', foldPred)

res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df$OR_or_beta <- NA
df$CI_low <- NA
df$CI_up <- NA
for(i in 1:nrow(df)){
  
  if(df$type_pheno[i] == 'CONTINUOUS'){
    df$OR_or_beta[i] <- df$beta[i]
    df$CI_low[i] <- df$OR_or_beta[i] + qnorm(0.025)*df$se_beta[i]
    df$CI_up[i] <- df$OR_or_beta[i] + qnorm(0.975)*df$se_beta[i]
  }else{
    df$OR_or_beta[i] <- exp(df$beta[i])
    df$CI_low[i] <- exp(df$beta[i] + qnorm(0.025)*df$se_beta[i])
    df$CI_up[i] <- exp(df$beta[i] + qnorm(0.975)*df$se_beta[i])
  }
}

# make plots
id_keep <- unique(df$pheno_id[df$pvalue <= 0.001])
df_red <- df[df$pheno_id %in% id_keep, ]

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1)+
  facet_wrap(comp~.,  nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red1$color)+
  scale_y_continuous(trans='log2')+
  coord_flip()

pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0)+
  facet_wrap(comp~., nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red2$color)+
  coord_flip()
len_h <- length(unique(df_red$pheno_id[df_red$type_pheno == 'CONTINUOUS'])) + length(unique(df_red$pheno_id[df_red$type_pheno != 'CONTINUOUS']))

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.3))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%sstscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1.5, plot = tot_pl, device = 'pdf')


### plot specific phenotypes ###
# ever smoked

pheno_id <- c('20160')
for(i in 1:length(pheno_id)){
  
  tmp_pheno <- data.frame(value = res$phenoDat[, pheno_id[i]], gr = paste0('gr',res$cl$gr))
  tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res$cl$gr))))
  stat_tmp <- res$bin_reg[res$bin_reg$pheno_id == pheno_id[i], colnames(res$bin_reg) %in% c('comp', 'pheno_id', 'Field', 'pvalue')]
  stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  stat_tmp <- stat_tmp[grepl('gr4',stat_tmp$group1) | grepl('gr4',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
  stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))
  
  stat_tmp$y.position <- mean(sapply(sort(unique(res$cl$gr)), function(x) max(tmp_pheno$value[tmp_pheno$gr == paste0('gr',x)], na.rm = T))) + 
    0.1+ 0.1*(0:(nrow(stat_tmp)-1))
  stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
  stat_tmp <- as_tibble(stat_tmp)        
  title_name <- unique(df$Field[df$pheno_id == pheno_id[i]])
  
  class_tmp <- unique(sort(tmp_pheno$value))
  new_df <- data.frame(perc = as.vector(table(tmp_pheno$gr, tmp_pheno$value)/rowSums(table(tmp_pheno$gr, tmp_pheno$value))), class = unlist(lapply(class_tmp, function(x) rep(x, length(unique(tmp_pheno$gr))))))
  new_df$gr <- rep(unique(sort(tmp_pheno$gr)), length(class_tmp))
  new_df$gr <- factor(new_df$gr, levels = unique(sort(tmp_pheno$gr)))
  new_df$class <- factor(new_df$class, levels = class_tmp)

  pl <- ggbarplot(new_df, x = "gr", y = "perc", fill = "class", palette = "Paired",
                  label = F, lab.col = "white", lab.pos = "in")+ 
    stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)
 
  pl <- ggpar(pl, legend = "top", xlab = '', ylab = paste0(title_name, '\npercentage'))
  file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise', fold)
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')
  
}


# load results from prediction:
res_pred <- get(load(glmPredFile))

cohort_name <- paste0('German', 1:5)

for(id in c(1,5)){
  
  if(id==1){
    pheno_id <- c('Smoker_ever', 'Smoker')
  }
  if(id==5){
    pheno_id <- c('Ever_Smoked', 'Previous_Smoker')
  }
  
  
  for(i in 1:length(pheno_id)){
    
    tmp_pheno <- data.frame(value = res_pred$phenoDat[[id]][, pheno_id[i]], gr = paste0('gr',res_pred$cl[[id]]$gr))
    tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res_pred$cl[[id]]$gr))))
    stat_tmp <- res_pred$bin_reg[[id]][res_pred$bin_reg[[id]]$pheno_id == pheno_id[i], colnames(res_pred$bin_reg[[id]]) %in% c('comp', 'pheno_id', 'pvalue')]
    stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
    stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
    stat_tmp <- stat_tmp[grepl('gr4',stat_tmp$group1) | grepl('gr4',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
    stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))
    
    vect_add <- 0.1+ 0.1*(0:(nrow(stat_tmp)-1))
    
    stat_tmp$y.position <- 1 + vect_add
    stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
    stat_tmp <- as_tibble(stat_tmp)        
    title_name <- pheno_id[i]
    
    class_tmp <- unique(sort(tmp_pheno$value))
    new_df <- data.frame(perc = as.vector(table(tmp_pheno$gr, tmp_pheno$value)/rowSums(table(tmp_pheno$gr, tmp_pheno$value))), class = unlist(lapply(class_tmp, function(x) rep(x, length(unique(tmp_pheno$gr))))))
    new_df$gr <- rep(unique(sort(tmp_pheno$gr)), length(class_tmp))
    new_df$gr <- factor(new_df$gr, levels = unique(sort(tmp_pheno$gr)))
    new_df$class <- factor(new_df$class, levels = class_tmp)
    
    pl <- ggbarplot(new_df, x = "gr", y = "perc", fill = "class", palette = "Paired",
                    label = F, lab.col = "white", lab.pos = "in")+ 
      stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)
    
    pl <- ggpar(pl, legend = "top", xlab = '', ylab = paste0(title_name, '\npercentage'))
    
    file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_prediction_modelUKBB_%s', foldPred, cohort_name[id])
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')
    
  }
}


#######################################################################################################
tissue <- 'Adrenal_Gland'
fold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
foldPred <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
glmPredFile <- sprintf('%stscore_zscaled_clusterCases_phenoAssociationGLM_prediction_modelUKBB.RData', foldPred)

res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df$OR_or_beta <- NA
df$CI_low <- NA
df$CI_up <- NA
for(i in 1:nrow(df)){
  
  if(df$type_pheno[i] == 'CONTINUOUS'){
    df$OR_or_beta[i] <- df$beta[i]
    df$CI_low[i] <- df$OR_or_beta[i] + qnorm(0.025)*df$se_beta[i]
    df$CI_up[i] <- df$OR_or_beta[i] + qnorm(0.975)*df$se_beta[i]
  }else{
    df$OR_or_beta[i] <- exp(df$beta[i])
    df$CI_low[i] <- exp(df$beta[i] + qnorm(0.025)*df$se_beta[i])
    df$CI_up[i] <- exp(df$beta[i] + qnorm(0.975)*df$se_beta[i])
  }
}

# make plots
id_keep <- unique(df$pheno_id[df$pvalue <= 0.001])
df_red <- df[df$pheno_id %in% id_keep, ]

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1)+
  facet_wrap(comp~.,  nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red1$color)+
  scale_y_continuous(trans='log2')+
  coord_flip()

pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_beta, color = pheno_type, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0)+
  facet_wrap(comp~., nrow = 1, strip.position="top")+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red2$color)+
  coord_flip()
len_h <- length(unique(df_red$pheno_id[df_red$type_pheno == 'CONTINUOUS'])) + length(unique(df_red$pheno_id[df_red$type_pheno != 'CONTINUOUS']))

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.3))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+3, height = len_h*0.2+2.5, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%sstscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+3, height = len_h*0.2+2.5, plot = tot_pl, device = 'pdf')

pheno_id <- c('20003_1140883066')
for(i in 1:length(pheno_id)){
  
  tmp_pheno <- data.frame(value = res$phenoDat[, pheno_id[i]], gr = paste0('gr',res$cl$gr))
  tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res$cl$gr))))
  stat_tmp <- res$bin_reg[res$bin_reg$pheno_id == pheno_id[i], colnames(res$bin_reg) %in% c('comp', 'pheno_id', 'Field', 'pvalue')]
  stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  stat_tmp <- stat_tmp[grepl('gr4',stat_tmp$group1) | grepl('gr4',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
  stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))
  
  stat_tmp$y.position <- mean(sapply(sort(unique(res$cl$gr)), function(x) max(tmp_pheno$value[tmp_pheno$gr == paste0('gr',x)], na.rm = T))) + 
    0.1+ 0.1*(0:(nrow(stat_tmp)-1))
  stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
  stat_tmp <- as_tibble(stat_tmp)        
  title_name <- ifelse(!is.na(unique(df$meaning[df$pheno_id == pheno_id[i]])), 
                       paste(unique(df$Field[df$pheno_id == pheno_id[i]]), unique(df$meaning[df$pheno_id == pheno_id[i]]), sep = ': '), unique(df$Field[df$pheno_id == pheno_id[i]]))  
  
  class_tmp <- unique(sort(tmp_pheno$value))
  new_df <- data.frame(perc = as.vector(table(tmp_pheno$gr, tmp_pheno$value)/rowSums(table(tmp_pheno$gr, tmp_pheno$value))), class = unlist(lapply(class_tmp, function(x) rep(x, length(unique(tmp_pheno$gr))))))
  new_df$gr <- rep(unique(sort(tmp_pheno$gr)), length(class_tmp))
  new_df$gr <- factor(new_df$gr, levels = unique(sort(tmp_pheno$gr)))
  new_df$class <- factor(new_df$class, levels = class_tmp)
  
  pl <- ggbarplot(new_df, x = "gr", y = "perc", fill = "class", palette = "Paired",
                  label = F, lab.col = "white", lab.pos = "in")+ 
    stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)
  
  pl <- ggpar(pl, legend = "top", xlab = '', ylab = paste0(title_name, '\npercentage'))
  file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise', fold)
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
  ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')
  
}


# load results from prediction:
res_pred <- get(load(glmPredFile))

cohort_name <- paste0('German', 1:5)

for(id in c(1,5)){
  
  if(id==1){
    pheno_id <- c('Diabetes')
  }
  if(id==5){
    pheno_id <- c('DIAB')
  }
  
  
  for(i in 1:length(pheno_id)){
    
    tmp_pheno <- data.frame(value = res_pred$phenoDat[[id]][, pheno_id[i]], gr = paste0('gr',res_pred$cl[[id]]$gr))
    tmp_pheno$gr <- factor(tmp_pheno$gr, levels = paste0('gr', sort(unique(res_pred$cl[[id]]$gr))))
    stat_tmp <- res_pred$bin_reg[[id]][res_pred$bin_reg[[id]]$pheno_id == pheno_id[i], colnames(res_pred$bin_reg[[id]]) %in% c('comp', 'pheno_id', 'pvalue')]
    stat_tmp$group1 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
    stat_tmp$group2 <- sapply(stat_tmp$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
    stat_tmp <- stat_tmp[grepl('gr4',stat_tmp$group1) | grepl('gr4',stat_tmp$group2) | grepl('gr3',stat_tmp$group1) | grepl('gr3',stat_tmp$group2), ]
    stat_tmp <- cbind(stat_tmp, data.frame(.y. = rep('value', nrow(stat_tmp))))
    
    vect_add <- 0.1+ 0.1*(0:(nrow(stat_tmp)-1))
    
    stat_tmp$y.position <- 1 + vect_add
    stat_tmp$p.format <- scientific(stat_tmp$pvalue, digits = 2)
    stat_tmp <- as_tibble(stat_tmp)        
    title_name <- pheno_id[i]
    
    class_tmp <- unique(sort(tmp_pheno$value))
    new_df <- data.frame(perc = as.vector(table(tmp_pheno$gr, tmp_pheno$value)/rowSums(table(tmp_pheno$gr, tmp_pheno$value))), class = unlist(lapply(class_tmp, function(x) rep(x, length(unique(tmp_pheno$gr))))))
    new_df$gr <- rep(unique(sort(tmp_pheno$gr)), length(class_tmp))
    new_df$gr <- factor(new_df$gr, levels = unique(sort(tmp_pheno$gr)))
    new_df$class <- factor(new_df$class, levels = class_tmp)
    
    pl <- ggbarplot(new_df, x = "gr", y = "perc", fill = "class", palette = "Paired",
                    label = F, lab.col = "white", lab.pos = "in")+ 
      stat_pvalue_manual(stat_tmp, label = "p.format", xmin = "group1", xmax = "group2", bracket.nudge.y = 0.5, label.size = 3)
    
    pl <- ggpar(pl, legend = "top", xlab = '', ylab = paste0(title_name, '\npercentage'))
    
    file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_prediction_modelUKBB_%s', foldPred, cohort_name[id])
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i],'.png'), width = 3, height = 5, dpi = 320, device = 'png')
    ggsave(pl, filename = paste0(file_name, '_boxplot_pheno', pheno_id[i], '.pdf'), width = 3, height = 5, device = 'pdf')
    
  }
}

##############################################################################################################################################


