# plot for publication CMC clustering

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot results group specific differences")
parser$add_argument("--riskScore_ann_file", type = "character", help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--phenoInfo_file", type = "character", help = "")
parser$add_argument("--pheno_plot", type = "character", help = "")
parser$add_argument("--measureGoodness_thr", type = "integer", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
riskScore_ann_file <- args$riskScore_ann_file
color_pheno_file <- args$color_pheno_file
pheno_name <- args$pheno_name
phenoInfo_file <- args$phenoInfo_file
pheno_plot <- args$pheno_plot
measureGoodness_thr <- args$measureGoodness_thr
outFold <- args$outFold

# ###################################################################################################################
# pheno_name <- 'SCZ'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/color_pheno_type_UKBB.txt'
# riskScore_ann_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/SCZ_clustering/matchUKBB_updated_riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis_annotated.txt'
# outFold <- paste0('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/SCZ_clustering/matchUKBB_updated')
# phenoInfo_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/phenotypeDescription_rsSCZ_updated.txt'
# pheno_plot <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/SCZ_clustering/list_phenoid_plot'
# ################################################################################################################

pheno_ann <- read.delim(color_pheno_file, header = T, stringsAsFactors = F)
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4', 'brown', 'brown','chartreuse4', 'brown'), pheno_type = c('ICD9-10_OPCS4', 'Medications', 'Medication',
                                                                                                                                    'Medical_conditions', 'Alcohol', 'Asthma_related_drugs')))
pheno_ann <- rbind(pheno_ann, data.frame(color = c('orange', 'firebrick2'), 
                                         pheno_type = c('dMRI_skeleton', 'Numeric_memory')))
pheno_ann$color[pheno_ann$pheno_type == 'T1_structural_brain_MRI'] <- 'grey40'
pheno_ann$color[pheno_ann$pheno_type == 'Blood_biochemistry'] <- 'blue4'
pheno_ann$color[pheno_ann$pheno_type == 'Blood_count'] <- 'hotpink4'
pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
pheno_ann$color[pheno_ann$pheno_type == 'Trail_making'] <- 'lightpink4'
pheno_ann$color[pheno_ann$pheno_type %in% c('ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')] <- 'grey40'

phenoInfo <- read.delim(phenoInfo_file, header = T, stringsAsFactors = F, sep = '\t')
if(!'pheno_type' %in% colnames(phenoInfo)){
  tmp_name <- sapply(phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
  tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
  phenoInfo$pheno_type <- tmp_name
  phenoInfo$pheno_type[phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
}

rs_res <- read.delim(riskScore_ann_file, h=T, stringsAsFactors = F, sep = '\t')
rs_res <- rs_res[!is.na(rs_res$pvalue), ]

###################################
### plot number of association ####
###################################
rs_res$comp <- factor(rs_res$comp)
rs_res$sign <- 'no'
rs_res$sign[rs_res$pval_corr <= 0.05] <- 'yes'

thr_lims <- seq(1000, 5000, by = 500)
df <- list()
for(i in 1:length(thr_lims)){
  df[[i]] <- rs_res %>% filter(measure >=thr_lims[i]) %>% group_by(comp) %>% count(sign, .drop = FALSE) %>% add_column(CRM_thr = thr_lims[i])
}
df <- bind_rows(df)
df <- df %>% mutate(gr = strsplit(as.character(comp), split = '_vs_all')[[1]])
df$gr <- factor(df$gr)
df$CRM_thr <- factor(df$CRM_thr, levels = thr_lims)

pl_count <-  ggplot(df, aes(x = CRM_thr, y = n, fill = gr))+
  geom_bar(stat = 'identity', width = 0.7,
           position=position_dodge())+
  theme_bw()+ 
  ylab('n. of significant\nphenotypes')+ 
  xlab('cluster-reliable measure = |beta|*F-stat')+ 
  theme(legend.position = 'right', axis.title.y = element_text(size = 9), axis.title.x = element_text(size=9),
        axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))+
  scale_fill_d3()
ggsave(filename = sprintf('%sriskScores_tscore_zscaled_cluster%s_GLM_count.png', outFold, 'Cases'), width = 5, height = 3, plot = pl_count, device = 'png', dpi = 300)
ggsave(filename = sprintf('%sriskScores_tscore_zscaled_cluster%s_GLM_count.pdf', outFold, 'Cases'), width = 5, height = 3, plot = pl_count, device = 'pdf')


#####################################################
### plot selection of association among 1500 thr ####
#####################################################
pheno_to_plot <- read.table(pheno_plot, header = F, stringsAsFactors = F)$V1
# names_red <- unique(rs_res$pheno_id[rs_res$pval_corr <= 0.05 & rs_res$measure>=1500])
rs_tmp <- rs_res[rs_res$pheno_id %in% pheno_to_plot, ]
rs_tmp$new_id <- paste0(rs_tmp$comp, '_and_', rs_tmp$pheno_id)
P <- length(unique(rs_tmp$comp))
gr_color <- pal_d3(palette = 'category20')(P)

df_red <- rs_tmp
df_red$type_m <- 'not reliable'
df_red$type_m[df_red$measure >= measureGoodness_thr] <- 'reliable'
df_red$new_id <- df_red$Field
df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$type_m <- factor(df_red$type_m, levels = c('not reliable', 'reliable'))
df_red$type_res <- 'beta'
df_red$sign <- 'no'
df_red$sign[df_red$pval_corr <= 0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
# divide in 2 part
df_red1 <- subset(df_red, pheno_type %in% c('Blood_count', 'Blood_biochemistry', 'Blood_count_ratio'))
df_red2 <- subset(df_red, !pheno_type %in% c('Blood_count', 'Blood_biochemistry', 'Blood_count_ratio'))
pheno_ann_red1 <- pheno_ann[match(df_red1$pheno_type, pheno_ann$pheno_type), ]
pheno_ann_red2 <- pheno_ann[match(df_red2$pheno_type, pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
len_h <- length(unique(df_red$pheno_id))

# change labels 
labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
names(labs_new) <- as.character(unique(df_red$comp))

pl_beta1 <-  ggplot(df_red1, aes(x = new_id, y = OR_or_Beta, shape = sign, color = type_m))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), size = 0.3, width=.3, position=position_dodge(0.05), show_guide=FALSE, linetype="solid")+
  geom_point(position=position_dodge(0.05), size = 1)+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red1$color), 
        strip.text = element_text(size=8, color = 'white', face = 'bold'))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c('grey', 'black'))+
  coord_flip()

pl_beta1 <- ggplot_gtable(ggplot_build(pl_beta1))
stripr <- which(grepl('strip-t', pl_beta1$layout$name))
fills <- gr_color
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', pl_beta1$grobs[[i]]$grobs[[1]]$childrenOrder))
  pl_beta1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


pl_beta2 <-  ggplot(df_red2, aes(x = new_id, y = OR_or_Beta, shape = sign, color = type_m))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), size = 0.3, width=.3, position=position_dodge(0.05), show_guide=FALSE, linetype="solid")+
  geom_point(position=position_dodge(0.05), size = 1)+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red2$color), 
        strip.text = element_text(size=8, color = 'white', face = 'bold'))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c('grey', 'black'))+
  coord_flip()

pl_beta2 <- ggplot_gtable(ggplot_build(pl_beta2))
stripr <- which(grepl('strip-t', pl_beta2$layout$name))
fills <- gr_color
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', pl_beta2$grobs[[i]]$grobs[[1]]$childrenOrder))
  pl_beta2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tot_pl <- ggarrange(plotlist = list(pl_beta1, pl_beta2), ncol = 1, nrow = 2, align='v', heights = c(1, 0.7))
ggsave(filename = sprintf('%sriskScores_tscore_zscaled_cluster%s_GLM_beta_measureThr%s_selected.png', outFold, 'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.1+2, plot = tot_pl, device = 'png', dpi=300)
ggsave(filename = sprintf('%sriskScore_tscore_zscaled_cluster%s_GLM_beta_measureThr%s_selected.pdf', outFold,'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.1+2, plot = tot_pl, device = 'pdf')



