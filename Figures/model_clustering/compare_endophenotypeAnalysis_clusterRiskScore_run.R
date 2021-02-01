# charachterize cluster results: prove on CAD 

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
suppressPackageStartupMessages(library(pROC))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="compare cluster risk score and endophenotype analysis")
parser$add_argument("--riskScore_analysis_file", type = "character", nargs = '*', help = "")
parser$add_argument("--endopheno_analysis_file", type = "character", help = "")
parser$add_argument("--pval_FDR_pheno", type = "double", default = 0.05, help = "")
parser$add_argument("--pval_pheno_show", type = "double", default = 0.001, help = "")
parser$add_argument("--thr_plot", type = "double", default = 1e-10, help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--R2_pheno_rs_file", type = "character", help = "")
parser$add_argument("--meta_analysis", type = "logical",default = F,  help = "")
parser$add_argument("--measureGoodness_thr", type = "double", default = 3000, help = "")
parser$add_argument("--outFold", type = "character", help = "")

args <- parser$parse_args()
riskScore_analysis_file <- args$riskScore_analysis_file
endopheno_analysis_file <- args$endopheno_analysis_file
pval_FDR_pheno <- args$pval_FDR_pheno
pval_pheno_show <- args$pval_pheno_show
color_pheno_file <- args$color_pheno_file
pheno_name <- args$pheno_name
thr_plot <- args$thr_plot
meta_analysis <- args$meta_analysis
R2_pheno_rs_file <- args$R2_pheno_rs_file
measureGoodness_thr <- args$measureGoodness_thr
outFold <- args$outFold

#####################################################################################################
# setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/')
# tissue <- 'Colon_Sigmoid'
# riskScore_analysis_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis.RData', tissue)
# endopheno_analysis_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', tissue)
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
# pval_FDR_pheno <- 0.05
# pval_pheno_show <- 0.001
# thr_plot <- 1e-50
# meta_analysis <- T
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# R2_pheno_rs_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_corrThr0.5_relatedPhenotypes_R2_risk_score_phenotype.txt', tissue)
# measureGoodness_thr <- 3000
####################################################################################################

pheno_ann <- read.delim(color_pheno_file, header = T, stringsAsFactors = F)
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4', 'brown', 'brown','chartreuse4', 'brown'), pheno_type = c('ICD9-10_OPCS4', 'Medications', 'Medication',
                                                                                                                                    'Medical_conditions', 'Alcohol', 'Asthma_related_drugs')))
pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'
pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
pheno_ann$color[pheno_ann$pheno_type %in% c('ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')] <- 'grey40'

endo_res <- read.delim(endopheno_analysis_file, h=T, stringsAsFactors = F, sep = '\t')
endo_res_tot <- endo_res[!is.na(endo_res$pvalue), ]
endo_res <- endo_res_tot[endo_res_tot$pval_corr <= pval_FDR_pheno | endo_res_tot$pvalue <= pval_pheno_show,  ]

rs_res <- list()

for(i in 1:length(riskScore_analysis_file)){
  tmp <- get(load(riskScore_analysis_file[i]))
  if(meta_analysis){
    rs_res[[i]] <- tmp$meta_analysis
  }else{
    rs_res[[i]] <- tmp$bin_reg
  }
  
  if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
    tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
    tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
    tmp$phenoInfo$pheno_type <- tmp_name
    tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
  }
  rs_res[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(rs_res[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
}
rs_res <- do.call(rbind, rs_res)
if(length(riskScore_analysis_file)>1){
  comp <- unique(rs_res$comp)
  tmp <- list()
  for(i in 1:length(comp)){
    tmp[[i]] <- rs_res[rs_res$comp == comp[i],]
    tmp[[i]]$pval_corr <- p.adjust(tmp[[i]]$pvalue, method = 'BH')
  }
  rs_res <- do.call(rbind, tmp)
  # # save results
  # write.table(x = rs_res, 
  #             file = sprintf('%s%s_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', outFold, type_cluster_data , type_cluster), 
  #             col.names = T, row.names = F, sep = '\t', quote = F)
}


# subset with endophenotypes associations
pheno_common <- intersect(rs_res$pheno_id, endo_res$pheno_id)
rs_res_red <- rs_res[rs_res$pheno_id %in% pheno_common, ]

P <- length(unique(rs_res_red$comp))
gr_color <- pal_d3(palette = 'category20')(P)

df_red <- rs_res_red
df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr <= 0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
df_red$type_res <- 'beta'
pheno_ann_red <- pheno_ann[match(df_red$pheno_type, pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
len_h <- length(unique(df_red$pheno_id))
# change labels 
labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
names(labs_new) <- as.character(unique(df_red$comp))

pl_beta <-  ggplot(df_red, aes(x = new_id, y = OR_or_Beta, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
        strip.text = element_text(size=8, color = 'white', face = 'bold'))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red$color)+
  coord_flip()

pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
stripr <- which(grepl('strip-t', pl_beta$layout$name))
fills <- gr_color
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
  pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

ggsave(filename = sprintf('%sriskScores_tscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_beta.png', outFold, 'Cases'), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'png')
ggsave(filename = sprintf('%sriskScore_tscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_beta.pdf', outFold, 'Cases'), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'pdf')


######################
# plot associations
rs_res <- rs_res[!is.na(rs_res$pval_corr),]
pheno_sign <- unique(rs_res$pheno_id[rs_res$pval_corr <= thr_plot])

if(length(pheno_sign)>0){
rs_res_red <- rs_res[rs_res$pheno_id %in% pheno_sign, ]

df_red <- rs_res_red
df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr <= 0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
df_red$type_res <- 'beta'
pheno_ann_red <- pheno_ann[match(df_red$pheno_type, pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
len_h <- length(unique(df_red$pheno_id))
# change labels 
labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
names(labs_new) <- as.character(unique(df_red$comp))

pl_beta <-  ggplot(df_red, aes(x = new_id, y = OR_or_Beta, shape = sign))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
        strip.text = element_text(size=8, color = 'white', face = 'bold'))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=pheno_ann_red$color)+
  coord_flip()

pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
stripr <- which(grepl('strip-t', pl_beta$layout$name))
fills <- gr_color
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
  pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

ggsave(filename = sprintf('%sthr%s_riskScores_tscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_beta.png', outFold, as.character(thr_plot), 'Cases'), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'png')
ggsave(filename = sprintf('%sthr%s_riskScore_tscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_beta.pdf', outFold, as.character(thr_plot), 'Cases'), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'pdf')
}
# compute overall concordance rate: significant associations for risk scores
thr_pval <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-10, 1e-50,  1e-75, 1e-100, 1e-150)
df_repr <- list()
comp_name <- as.character(unique(rs_res$comp))

for(j in 1:length(thr_pval)){
  print(j)
  df_repr[[j]] <- data.frame(comp = comp_name, n_pheno = NA, n_same_sign = NA, fraction_rep = NA, pval = NA, thr_corr = thr_pval[j])
  
  for(i in 1:length(comp_name)){
    
    id <- rs_res$pvalue <= thr_pval[j] & rs_res$comp %in% comp_name[i]
    pheno_c <- intersect(rs_res$pheno_id[id], endo_res_tot$pheno_id)
    tmp <- rs_res[id, ]
    zstat <- tmp$z[match(pheno_c, tmp$pheno_id)]
    beta_endo <- endo_res_tot[endo_res_tot$comp %in% comp_name[i], ]
    beta_endo <- beta_endo$z[match(pheno_c, beta_endo$pheno_id)]
    df_repr[[j]]$n_pheno[i] <- length(pheno_c)
    df_repr[[j]]$n_same_sign[i] <- sum(sign(zstat*beta_endo) == 1)
    if(length(pheno_c) > 0){
      df_repr[[j]]$fraction_rep[i] <- sum(sign(zstat*beta_endo) == 1)/length(pheno_c)
      df_repr[[j]]$pval[i] <- binom.test(x = df_repr[[j]]$n_same_sign[i], n = df_repr[[j]]$n_pheno[i])$p.value
    }
    
  }
  
}

df_repr <- do.call(rbind, df_repr)

# plot distribution and save table
write.table(file = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_concordance.txt', outFold, pheno_name), x = df_repr, quote = F, sep = '\t', row.names = F, col.names = T)

# df_repr$new_thr <- -log10(df_repr$thr_corr_FDR)
df_repr$new_thr <- factor(df_repr$thr_corr, levels = thr_pval)
df_repr$new_pval <- -log10(df_repr$pval)
df_repr$comp <- factor(df_repr$comp, levels = comp_name)


pl <- ggplot(data = df_repr, aes(x = new_thr, y = fraction_rep, color = comp, group = comp, size = new_pval))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('risk score p-value threshold')+ylab('Fraction of concordance')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  # scale_x_continuous(breaks=-log10(thr_pval), labels = thr_pval)+
  scale_size_continuous(breaks = seq(min(df_repr$new_pval, na.rm = T), max(df_repr$new_pval, na.rm = T), length.out = 5), 
                        labels = formatC(10^(-seq(min(df_repr$new_pval, na.rm = T), max(df_repr$new_pval, na.rm = T), length.out = 5)), format = "e", digits = 2))+
  guides(color=guide_legend(title = ''), size = guide_legend(title = 'binomial test\np-value'))+
  scale_color_d3()

ggsave(filename = sprintf('%sriskScores_clusterCases_%s_groupCorrelation_relatedPheno_concordance.png', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sriskScores_clusterCases_%s_groupCorrelation_relatedPheno_concordance.pdf', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500, compress = F)

##################################################################################
# evaluate wrt R2 
R2_pheno_rs <- read.delim(R2_pheno_rs_file, h=T, stringsAsFactors = F, sep = '\t')

common_pheno <- intersect(intersect(R2_pheno_rs$pheno_id, unique(endo_res_tot$pheno_id)), rs_res$pheno_id)
R2_pheno_rs <- R2_pheno_rs[match(common_pheno, R2_pheno_rs$pheno_id),]
endo_res_comp <- endo_res_tot[endo_res_tot$pheno_id %in% common_pheno,]
rs_res_comp <- rs_res[rs_res$pheno_id %in% common_pheno,]

comp_name <- unique(rs_res_comp$comp)
eval_sign_w <- list()
for(i in 1:length(comp_name)){
  
  tmp_rs <- rs_res_comp[rs_res_comp$comp %in% comp_name[i], ]
  tmp_rs <- tmp_rs[match(common_pheno, tmp_rs$pheno_id),]
  tmp_endo <- endo_res_comp[endo_res_comp$comp %in% comp_name[i] , ]
  tmp_endo <- tmp_endo[match(common_pheno, tmp_endo$pheno_id), ]
  tmp_eval <- data.frame(Fstat = R2_pheno_rs$Fstat_risk, R2 = R2_pheno_rs$R2_risk, R2_nsamples = R2_pheno_rs$nsamples,
                         beta_rs = tmp_rs$beta, beta_endo = tmp_endo$beta,
                         pheno_id = common_pheno, rs_sign = F, endo_sign = F)
  
  tmp_eval$rs_sign[tmp_rs$pval_corr <= 0.05] <- T
  tmp_eval$endo_sign[tmp_endo$pval_corr <= 0.05] <- T
  # tmp_eval$measure <- tmp_eval$R2*abs(tmp_eval$beta_rs)*sqrt(tmp_eval$R2_nsamples)
  tmp_eval$measure <- tmp_eval$Fstat*abs(tmp_eval$beta_rs)
  tmp_eval$same_sign <- 0
  tmp_eval$same_sign[sign(tmp_eval$beta_rs * tmp_eval$beta_endo) == 1] <- 1
  tmp_eval$comp <- comp_name[i]
  eval_sign_w[[i]] <- tmp_eval
  
}
eval_sign_w <- do.call(rbind, eval_sign_w)
write.table(file = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_measureGoodnessPred.txt', outFold, pheno_name), x = eval_sign_w, quote = F, sep = '\t', row.names = F, col.names = T)

# consider only significant differences in rs
eval_sign_w <- eval_sign_w[eval_sign_w$rs_sign, ]
tmp_eval <- eval_sign_w
tmp_eval$comp <- factor(tmp_eval$comp)
# plot Fstat, beta
pl <- ggplot(data =  tmp_eval, aes(x = Fstat, y = beta_rs, 
                                      color = measure, shape = comp))+
  geom_point(alpha = 0.7, size = 1)+
  xlab('F-stat')+ylab('beta risk-score')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = 'F-stat * |beta|'))+
  scale_color_gradient(low = 'blue', high = 'red', breaks = seq(round(min(eval_sign_w$measure)), round(max(eval_sign_w$measure)), length.out = 10), 
                       limits = c(min(eval_sign_w$measure), max(eval_sign_w$measure)))

ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_distBetaFstat.png', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_distBetaFstat.pdf', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500, compress = F)


thr_measure <- seq(min(eval_sign_w$measure),max(eval_sign_w$measure), length.out = 1000)
df_repr <- list()
gr_name <- sapply(comp_name, function(x) strsplit(x, split = '_vs_all')[[1]][1])

for(j in 1:length(thr_measure)){
  # print(j)
  df_repr[[j]] <- data.frame(comp = gr_name, n_pheno_thr = NA, 
                             n_same_sign = NA, precision = NA, recall = NA, 
                             thr_measure = thr_measure[j])
  
  for(i in 1:length(comp_name)){
    
    tmp <- eval_sign_w[eval_sign_w$comp == comp_name[i],]
    df_repr[[j]]$n_pheno_thr[i] <- sum(tmp$measure >= thr_measure[j])
    df_repr[[j]]$n_same_sign[i] <- sum(tmp$same_sign == 1)
    
    if(sum(tmp$measure >= thr_measure[j]) > 0){
      df_repr[[j]]$precision[i] <- sum(tmp$measure >= thr_measure[j] & tmp$same_sign == 1)/sum(tmp$measure >= thr_measure[j])
    }else{
      df_repr[[j]]$precision[i] <- NA
    }
    
    if(sum(tmp$same_sign == 1) > 0){
      df_repr[[j]]$recall[i] <- sum(tmp$measure >= thr_measure[j] & tmp$same_sign == 1)/sum(tmp$same_sign == 1)
    }else{
      df_repr[[j]]$recall[i] <- NA
    }
  }
  
}

df_repr <- do.call(rbind, df_repr)
df_repr$comp <- factor(df_repr$comp, levels = gr_name)

thr <- max(sapply(gr_name, function(x) min(df_repr$thr_measure[df_repr$precision == 1 & df_repr$comp == x], na.rm  = T)))
tmp <- subset(df_repr, thr_measure <= thr)
pl <- ggplot(data =  tmp,aes(x = thr_measure, y = precision, color = comp, group = comp, size = n_pheno_thr))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('F-stat * |beta|')+ylab('Precision\nsame beta sign')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''), size = guide_legend(title = 'n. pheno > thr'))+
  scale_color_d3()

ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_precision.png', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_precision.pdf', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500, compress = F)

pl <- ggplot(data =  tmp, aes(x = thr_measure, y = recall, color = comp, group = comp, size = n_same_sign))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('F-stat * |beta|')+ylab('Recall\nsame beta sign')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''), size = guide_legend(title = 'n. pheno\nsame beta sign'))+
  scale_color_d3()

ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_recall.png', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sriskScores_clusterCases_%s_group_relatedPheno_recall.pdf', outFold, pheno_name), plot = pl, width = 5.5, height = 4.5, dpi = 500, compress = F)

#### plot association for the fixed thr
tmp_eval <- eval_sign_w[eval_sign_w$measure >= measureGoodness_thr, ]
if(nrow(tmp_eval)>0){
  tmp_eval$new_id <- paste0(tmp_eval$comp, '_and_', tmp_eval$pheno_id)
  rs_res_red <- rs_res[rs_res$pheno_id %in% unique(tmp_eval$pheno_id), ]
  rs_res_red$new_id <- paste0(rs_res_red$comp, '_and_', rs_res_red$pheno_id)
  rs_res_red$reliable <- 'not pass threshold'
  rs_res_red$reliable[rs_res_red$new_id %in% tmp_eval$new_id] <- 'pass threshold'
  
  df_red <- rs_res_red
  df_red$new_id <- df_red$Field
  df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
  df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$reliable <- factor(df_red$reliable, levels = c('not pass threshold', 'pass threshold'))
  df_red$type_res <- 'beta'
  pheno_ann_red <- pheno_ann[match(df_red$pheno_type, pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red$comp))
  len_h <- length(unique(df_red$pheno_id))
  # change labels 
  labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  names(labs_new) <- as.character(unique(df_red$comp))
  
  pl_beta <-  ggplot(df_red, aes(x = new_id, y = OR_or_Beta, shape = reliable))+
    geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
          strip.text = element_text(size=8, color = 'white', face = 'bold'))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=pheno_ann_red$color)+
    coord_flip()
  
  pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
  stripr <- which(grepl('strip-t', pl_beta$layout$name))
  fills <- gr_color
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
    pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  ggsave(filename = sprintf('%sriskScores_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.png', outFold, 'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'png')
  ggsave(filename = sprintf('%sriskScore_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.pdf', outFold,'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'pdf')
  
  # original endo
  endo_res_red <- endo_res_tot[endo_res_tot$pheno_id %in% unique(tmp_eval$pheno_id), ]
  endo_res_red$new_id <- paste0(endo_res_red$comp, '_and_', endo_res_red$pheno_id)
  endo_res_red$reliable <- 'not pass threshold'
  endo_res_red$reliable[endo_res_red$new_id %in% tmp_eval$new_id] <- 'pass threshold'
  
  df_red <- endo_res_red
  df_red$new_id <- df_red$Field
  df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
  df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$reliable <- factor(df_red$reliable, levels = c('not pass threshold', 'pass threshold'))
  df_red$type_res <- 'beta'
  id_mod <- df_red$type_pheno != 'CONTINUOUS'
  if(any(id_mod)){
    df_red$CI_low[id_mod] <- log(df_red$CI_low[id_mod])
    df_red$CI_up[id_mod] <- log(df_red$CI_up[id_mod])
  }
  
  pheno_ann_red <- pheno_ann[match(df_red$pheno_type, pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red$comp))
  len_h <- length(unique(df_red$pheno_id))
  # change labels 
  labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  names(labs_new) <- as.character(unique(df_red$comp))
  
  pl_beta <-  ggplot(df_red, aes(x = new_id, y = beta, shape = reliable))+
    geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
          strip.text = element_text(size=8, color = 'white', face = 'bold'))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=pheno_ann_red$color)+
    coord_flip()
  
  pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
  stripr <- which(grepl('strip-t', pl_beta$layout$name))
  fills <- gr_color
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
    pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  ggsave(filename = sprintf('%soriginalEndo_riskScores_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.png', outFold, 'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'png')
  ggsave(filename = sprintf('%soriginalEndo_riskScores_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.pdf', outFold,'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'pdf')
  
}


