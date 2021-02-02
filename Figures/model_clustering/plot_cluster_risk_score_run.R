# plot endophenotype analysis from matrix results
# specific for risk score: combaine with goodness of fit

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
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot results group specific differences")
parser$add_argument("--riskScore_res_file", type = "character", help = "")
parser$add_argument("--R2_file", type = "character", help = "")
parser$add_argument("--measureGoodness_thr", type = "double", default = 3000, help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--phenoInfo_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
riskScore_res_file <- args$riskScore_res_file
R2_file <- args$R2_file
measureGoodness_thr <- args$measureGoodness_thr
color_pheno_file <- args$color_pheno_file
pheno_name <- args$pheno_name
phenoInfo_file <- args$phenoInfo_file
outFold <- args$outFold

###################################################################################################################
# tissue <- 'Brain_Caudate_basal_ganglia'
# pheno_name <- 'SCZ'
# measureGoodness_thr <- 300
# color_pheno_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/color_pheno_type_UKBB.txt'
# riskScore_res_file <- paste0('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/clustering_res/',tissue,'/excludeMHC_riskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis.txt')
# outFold <- paste0('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/',tissue,'/SCZ_clustering/excludeMHC_')
# R2_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/predict_UKBB/Brain_Caudate_basal_ganglia/200kb/noGWAS/devgeno0.01_testdevgeno0/tscore_corrThr0.5_relatedPhenotypes_R2_risk_score_phenotype.txt'
# phenoInfo_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/phenotypeDescription_rsSCZ.txt'
#################################################################################################################

pheno_ann <- read.delim(color_pheno_file, header = T, stringsAsFactors = F)
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4', 'brown', 'brown','chartreuse4', 'brown'), pheno_type = c('ICD9-10_OPCS4', 'Medications', 'Medication',
                                                                                                                                    'Medical_conditions', 'Alcohol', 'Asthma_related_drugs')))
pheno_ann <- rbind(pheno_ann, data.frame(color = c('orange', 'firebrick2'), 
                              pheno_type = c('dMRI_skeleton', 'Numeric_memory')))
pheno_ann$color[pheno_ann$pheno_type == 'T1_structural_brain_MRI'] <- 'grey40'
pheno_ann$color[pheno_ann$pheno_type == 'Blood_biochemistry'] <- 'blue4'
pheno_ann$color[pheno_ann$pheno_type == 'Blood_count'] <- 'hotpink4'
pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
pheno_ann$color[pheno_ann$pheno_type %in% c('ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')] <- 'grey40'

phenoInfo <- read.delim(phenoInfo_file, header = T, stringsAsFactors = F, sep = '\t')
if(!'pheno_type' %in% colnames(phenoInfo)){
  tmp_name <- sapply(phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
  tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
  phenoInfo$pheno_type <- tmp_name
  phenoInfo$pheno_type[phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
}

rs_res <- read.delim(riskScore_res_file, h=T, stringsAsFactors = F, sep = '\t')
rs_res <- rs_res[!is.na(rs_res$pvalue), ]
R2_pheno_rs <- read.delim(R2_file, h=T, stringsAsFactors = F, sep = '\t')
common_pheno <- intersect(R2_pheno_rs$pheno_id, unique(rs_res$pheno_id))
R2_pheno_rs <- R2_pheno_rs[match(common_pheno, R2_pheno_rs$pheno_id),]
rs_res_comp <- rs_res[rs_res$pheno_id %in% common_pheno,]

# add measure info
comp_name <- sort(unique(rs_res_comp$comp))
rs_res_withmeas <- list()
for(i in 1:length(comp_name)){
  
  tmp_rs <- rs_res_comp[rs_res_comp$comp %in% comp_name[i], ]
  c_pheno_tmp <- intersect(R2_pheno_rs$pheno_id, tmp_rs$pheno_id)
  tmp_rs <- tmp_rs[match(c_pheno_tmp, tmp_rs$pheno_id),]
  tmp_R2 <- R2_pheno_rs[match(c_pheno_tmp, R2_pheno_rs$pheno_id) , ]
  tmp_rs$R2_risk <- tmp_R2$R2_risk
  tmp_rs$Fstat_risk <- tmp_R2$Fstat_risk
  tmp_rs$measure <- tmp_rs$Fstat*abs(tmp_rs$beta)
  tmp_info <- phenoInfo[match(c_pheno_tmp, phenoInfo$pheno_id),]
  tmp_rs$pheno_type <- tmp_info$pheno_type
  rs_res_withmeas[[i]] <- tmp_rs
  
}
rs_res_withmeas <- do.call(rbind, rs_res_withmeas)
# save updated table
write.table(x = rs_res_withmeas, file = sprintf('%sriskScores_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_metaAnalysis_annotated.txt', outFold),
            col.names = T, row.names = F, sep = '\t', quote = F)
### plot
pheno_sign <- unique(rs_res_withmeas$pheno_id[(rs_res_withmeas$pval_corr <= 0.05 & rs_res_withmeas$measure >= measureGoodness_thr) |(rs_res_withmeas$pvalue <= 0.001 & rs_res_withmeas$measure >= measureGoodness_thr) ])
rs_res_withmeas$new_id <- paste0(rs_res_withmeas$comp, '_and_', rs_res_withmeas$pheno_id)
P <- length(unique(rs_res_withmeas$comp))
gr_color <- pal_d3(palette = 'category20')(P)

df_red <- rs_res_withmeas[rs_res_withmeas$pheno_id %in% pheno_sign, ]
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
pheno_ann_red <- pheno_ann[match(df_red$pheno_type, pheno_ann$pheno_type), ]

len_w <- length(unique(df_red$comp))
len_h <- length(unique(df_red$pheno_id))
# change labels 
labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
names(labs_new) <- as.character(unique(df_red$comp))

pl_beta <-  ggplot(df_red, aes(x = new_id, y = OR_or_Beta, shape = sign, color = type_m))+
  geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
  theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
        strip.text = element_text(size=8, color = 'white', face = 'bold'))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c('grey', 'black'))+
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


# reduced plot: exclude all blood count info
df_red_small <- df_red[df_red$pheno_type != 'Blood_count' | grepl('leukocyte',df_red$Field), ]
if(!identical(df_red, df_red_small)){
  pheno_ann_red <- pheno_ann[match(df_red_small$pheno_type, pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red_small$comp))
  len_h <- length(unique(df_red_small$pheno_id))
  # change labels 
  labs_new <- sapply(as.character(unique(df_red_small$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  names(labs_new) <- as.character(unique(df_red_small$comp))
  
  pl_beta <-  ggplot(df_red_small, aes(x = new_id, y = OR_or_Beta, shape = sign, color = type_m))+
    geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red$color), 
          strip.text = element_text(size=8, color = 'white', face = 'bold'))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=c('grey', 'black'))+
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
  
  ggsave(filename = sprintf('%sreduced_riskScores_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.png', outFold, 'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'png')
  ggsave(filename = sprintf('%sreduced_riskScore_tscore_zscaled_cluster%s_GLM_beta_measureThr%s.pdf', outFold,'Cases', as.character(measureGoodness_thr)), width = len_w+3, height = len_h*0.2+1, plot = pl_beta, device = 'pdf')
}



