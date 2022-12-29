#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Combine endophenotype res and plot")
parser$add_argument("--endopFile", type = "character", nargs = '*', help = "file to be loaded (endophenotypes association result)")
parser$add_argument("--colorFile", type = "character", help = "")
parser$add_argument("--type_cluster_data", type = "character", help = "")
parser$add_argument("--type_input", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character", help = "")
parser$add_argument("--forest_plot", type = "logical",default = F,  help = "")
parser$add_argument("--pval_pheno", type = "double", default = 0.0001, help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
endopFile <- args$endopFile
colorFile <- args$colorFile
type_cluster_data <- args$type_cluster_data
type_cluster <- args$type_cluster
type_input <- args$type_input
outFold <- args$outFold
forest_plot <- args$forest_plot
pval_pheno <- args$pval_pheno

#################################################################
# type_input <- 'corrPCs_zscaled'
# type_cluster <- 'Cases'
# type_cluster_data <- 'tscore'
# forest_plot <- T
# pval_pheno <- 0.001
# fold='OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/'
# outFold <- fold
# endopFile <- c(sprintf('%srescaleCont_withMedication_tscore_corrPCs_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', fold),
#                sprintf('%srescaleCont_withoutMedication_tscore_corrPCs_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', fold))
# colorFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
################################################################


res_pheno <- list()
for(i in 1:length(endopFile)){
  
  tmp <- get(load(endopFile[[i]]))
  res_pheno[[i]] <- tmp$bin_reg
  if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
    tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
    tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
    tmp$phenoInfo$pheno_type <- tmp_name
    tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
  }
  if(any(c('129', '130') %in% tmp$phenoInfo$pheno_id)){
    tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_id %in% c('129', '130')] <- 'Early_life_factors'
  }
  res_pheno[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(res_pheno[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
}

res_pheno <- do.call(rbind, res_pheno)
res_pheno <- res_pheno %>% group_by(comp) %>% 
  mutate(pval_corr = p.adjust(pvalue, method = 'BH')) %>% ungroup()

# save results
write.table(x = res_pheno, 
            file = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', outFold, type_cluster_data , type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

##############################
#### endophenotype plots #####
##############################
if(forest_plot){
  
  pheno_ann <- read.delim(colorFile, header = T, stringsAsFactors = F)

  id_keep <- unique(res_pheno$pheno_id[res_pheno$pvalue <= pval_pheno | res_pheno$pval_corr <= 0.05])
  id_keep <- id_keep[!is.na(id_keep)]
  df_red <- res_pheno %>% filter(pheno_id %in% id_keep) %>% 
    mutate(new_id = ifelse(is.na(meaning), paste(Field),paste(Field, meaning, sep = '\n'))) %>%
    mutate(sign = ifelse(pval_corr <= 0.05, 'yes', 'no')) %>%
    mutate(type_res = ifelse(type_pheno == 'CONTINUOUS', 'beta', 'OR')) %>%
    mutate(OR_or_Beta = ifelse(se_beta > 100, NA, OR_or_Beta)) %>%
    mutate(CI_low = ifelse(se_beta > 100, NA, CI_low)) %>% 
    mutate(CI_up = ifelse(se_beta > 100, NA, CI_up))
  
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  
  df_red_for_ann <- df_red[!duplicated(df_red$new_id),]
  pheno_ann_red1 <- pheno_ann[match(df_red_for_ann$pheno_type[df_red_for_ann$type_pheno != 'CONTINUOUS'], pheno_ann$pheno_type), ]
  pheno_ann_red2 <- pheno_ann[match(df_red_for_ann$pheno_type[df_red_for_ann$type_pheno == 'CONTINUOUS'], pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red$comp))
  len_h <- length(unique(df_red$pheno_id))
  # change labels 
  labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  names(labs_new) <- as.character(unique(df_red$comp))
  
  P <- length(unique(df_red$comp))
  gr_color <- pal_d3(palette = 'category20')(P)
  
  if(any(df_red$type_pheno != 'CONTINUOUS')){
    
    pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
      geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
      theme_bw()+ 
      ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
      facet_wrap(comp~.,  nrow = 1, strip.position="top", labeller = labeller(comp = labs_new))+
      theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7, colour = pheno_ann_red1$color),
            strip.text = element_text(size=8, color = 'white', face = 'bold'))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=pheno_ann_red1$color)+
      scale_y_continuous(trans='log', labels = scales::number_format(accuracy = 0.01))+
      coord_flip()
    
    pl_OR <- ggplot_gtable(ggplot_build(pl_OR))
    stripr <- which(grepl('strip-t', pl_OR$layout$name))
    fills <- gr_color
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', pl_OR$grobs[[i]]$grobs[[1]]$childrenOrder))
      pl_OR$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  
  if(any(df_red$type_pheno == 'CONTINUOUS')){
    
    pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
      geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
      theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red2$color), 
            strip.text = element_text(size=8, color = 'white', face = 'bold'))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=pheno_ann_red2$color)+
      coord_flip()
    ratio_OR_beta <- sum(df_red$type_pheno == 'CONTINUOUS')/sum(df_red$type_pheno != 'CONTINUOUS')
    
    pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
    stripr <- which(grepl('strip-t', pl_beta$layout$name))
    fills <- gr_color
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
      pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  
  if(any(df_red$type_pheno == 'CONTINUOUS') & any(df_red$type_pheno != 'CONTINUOUS')){
    tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, ratio_OR_beta))
  }else{
    if(any(df_red$type_pheno == 'CONTINUOUS')){
      tot_pl <- pl_beta
    }else{
      tot_pl <- pl_OR
    }
  }
  
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.png', outFold, type_cluster_data , type_input, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.pdf', outFold, type_cluster_data , type_input, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'pdf')
  
}


