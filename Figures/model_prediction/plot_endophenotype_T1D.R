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
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4', 'chartreuse4', 'chocolate2'), pheno_type = c('ICD9-10_OPCS4', 'Medications', 'Alcohol', 'Medication')))
pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'

#######################################################################################################
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
tissue <- 'Adipose_Subcutaneous'
fold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
withMed_glmFile <- sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)

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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
# remove comparison without significant pval
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
  
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')

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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
# remove comparison without significant pval
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
  
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')


#######################################################################################################
tissue <- 'Liver'
fold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
withMed_glmFile <- sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)

res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df <- df[!df$pheno_type %in% 'Alcohol',]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
df_red <- df_red[!is.na(df_red$pvalue),]
# remove comparison without significant pval
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  # print(i)
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')

### with Medication
res <- get(load(withMed_glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
# df <- df[!df$pheno_type %in% 'Alcohol',]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
# remove comparison without significant pval
df_red <- df_red[!is.na(df_red$pvalue),]
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
  
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')


#######################################################################################################
tissue <- 'Pancreas'
fold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
withMed_glmFile <- sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)

res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df <- df[!df$pheno_type %in% c('Alcohol', 'Smoking'),]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
df_red <- df_red[!is.na(df_red$pvalue),]
# remove comparison without significant pval
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  # print(i)
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')

### with Medication
res <- get(load(withMed_glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
# df <- df[!df$pheno_type %in% 'Alcohol',]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
# remove comparison without significant pval
df_red <- df_red[!is.na(df_red$pvalue),]
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
  
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')


#######################################################################################################
tissue <- 'Whole_Blood'
fold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissue)
glmFile <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)
withMed_glmFile <- sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', fold)

res <- get(load(glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
df <- df[!df$pheno_type %in% c('Alcohol', 'Smoking'),]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
df_red <- df_red[!is.na(df_red$pvalue),]
# remove comparison without significant pval
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  # print(i)
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')

### with Medication
res <- get(load(withMed_glmFile))

# compute odds ratio and 95% CI
df <- res$bin_reg
df$pheno_type <- res$phenoInfo$pheno_type[match(df$pheno_id,res$phenoInfo$pheno_id)]
# df <- df[!df$pheno_type %in% 'Alcohol',]
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
id_keep <- unique(df$pheno_id[df$pval_corr_overall <= 0.01])
df_red <- df[df$pheno_id %in% id_keep, ]
# remove comparison without significant pval
df_red <- df_red[!is.na(df_red$pvalue),]
name_c <- unique(df_red$comp)
for(i in 1:length(name_c)){
  tmp <- df_red[df_red$comp == name_c[i],]
  if(all(tmp$pval_corr_overall>0.01)){
    df_red <- df_red[!df_red$comp %in% name_c[i], ]
  }
  
}

df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr_overall<=0.01] <- 'yes'
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

tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, 1.6))

ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%swithMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', fold), width = len_w+2, height = len_h*0.2+1, plot = tot_pl, device = 'pdf')


#################################################
## venn diagram intersection significant genes ##
library('gdata')
library('VennDiagram')

table_genes <- read.xls('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/T1D_clustering/41588_2015_BFng3245_MOESM391_ESM.xls', h=T)
id_genes <- unique(table_genes$genes)
id_genes <- id_genes[id_genes != '']

# load genes results
tscore <- read.table('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/AllTissues/200kb/noGWAS/T1D_pheno/tscore_pval_T1D_covCorr.txt', h=T, stringsAsFactors=F)
tissues <- unique(tscore$tissue)

train_fold <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/'
train_fold <- paste0(train_fold, tissues, '/200kb/noGWAS/')   
# gene location
tscore$start_position <- NA
tscore$chrom <- NA

for(i in 1:length(train_fold)){
  
  tmp <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[i]), h=T,stringsAsFactors = F)
  tmp <- tmp[match(tscore$ensembl_gene_id[tscore$tissue == tissues[i]], tmp$ensembl_gene_id),]
  tscore$start_position[tscore$tissue == tissues[i]] <- tmp$start_position
  tscore$chrom[tscore$tissue == tissues[i]] <- tmp$chrom
  
}

HLA_reg <- c(28000000, 34000000)
tscore_red <- tscore
tscore_red <- tscore_red[!(tscore_red$chrom %in% 'chr6' & tscore_red$start_position <=HLA_reg[2] & tscore_red$start_position >= HLA_reg[1]) , ]
tscore_sign <- tscore_red[tscore_red$T1D_BHcorr_overall <= 0.05, ]
id_genes_new <- unique(tscore_sign$external_gene_name)


x <- list('Previously known Genes'=id_genes , 'PriLer Genes'= id_genes_new)
png('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/AllTissues/200kb/noGWAS/T1D_pheno/VennDiag_list_genes.png', units = 'in',
    width = 5, height = 5, res = 200)
v0 <- venn.diagram( x, filename=NULL,
                    fill = c("red", "blue"),
                    alpha = c(0.5, 0.5), cat.cex = 1, cex=1, cat.pos = 0)

overlaps <- calculate.overlap(x)

# extract indexes of overlaps from list names
indx <- as.numeric(substr(names(overlaps),2,2))
v0[[7]]$label  <- paste(overlaps[[3]], collapse="\n")  

grid.newpage()
grid.draw(v0)
dev.off()

pdf('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/AllTissues/200kb/noGWAS/T1D_pheno/VennDiag_list_genes.pdf', 
    width = 5, height = 5, compress = F)
v0 <- venn.diagram( x, filename=NULL,
                    fill = c("red", "blue"),
                    alpha = c(0.5, 0.5), cat.cex = 1, cex=1, cat.pos = 0)

overlaps <- calculate.overlap(x)

# extract indexes of overlaps from list names
indx <- as.numeric(substr(names(overlaps),2,2))
v0[[7]]$label  <- paste(overlaps[[3]], collapse="\n")  

grid.newpage()
grid.draw(v0)
dev.off()
