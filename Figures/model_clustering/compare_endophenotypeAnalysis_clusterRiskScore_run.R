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
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="compare cluster risk score and endophenotype analysis")
parser$add_argument("--riskScore_analysis_file", type = "character", nargs = '*', help = "")
parser$add_argument("--endopheno_analysis_file", type = "character", help = "")
parser$add_argument("--pval_FDR_pheno", type = "double", default = 0.05, help = "")
parser$add_argument("--pval_pheno_show", type = "double", default = 0.001, help = "")
parser$add_argument("--thr_plot", type = "double", default = 1e-10, help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--meta_analysis", type = "logical",default = F,  help = "")
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
outFold <- args$outFold

#####################################################################################################
# setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/')
# tissue <- 'Liver'
# riskScore_analysis_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis.RData', tissue)
# endopheno_analysis_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', tissue)
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/CAD_HARD_clustering/', tissue)
# pval_FDR_pheno <- 0.05
# pval_pheno_show <- 0.001
# thr_plot <- 1e-50
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# ####################################################################################################

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

# compute overall concordance rate: significant associations for risk scores
thr_pval <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-10, 1e-50,  1e-75, 1e-100, 1e-150)
df_repr <- list()
comp_name <- as.character(unique(df_red$comp))

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
write.table(file = sprintf('%sriskScores_clusterCases_%s_groupCorrelation_relatedPheno_concordance.txt', outFold, pheno_name), x = df_repr, quote = F, sep = '\t', row.names = F, col.names = T)

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




