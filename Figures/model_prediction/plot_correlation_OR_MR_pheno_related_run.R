# plot correlation with related phenotypes and Mendelian randomization results for CAD

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

parser <- ArgumentParser(description="correlation plots: tissue specific")
parser$add_argument("--corrRes_tissue_file", type = "character", nargs = '*', help = "results for correlation analysis tissue spec")
parser$add_argument("--mrResTscore_tissue_file", type = "character", nargs = '*', help = "results for mendelian randomization")
parser$add_argument("--mrResTscore_rev_tissue_file", type = "character", nargs = '*', help = "results for reverse mendelian randomization")
parser$add_argument("--mrResPath_tissue_file", type = "character", nargs = '*', help = "results for mendelian randomization")
parser$add_argument("--mrResPath_rev_tissue_file", type = "character", nargs = '*', help = "results for reverse mendelian randomization")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--pheno_list_file", type = "character", help = "")
parser$add_argument("--pheno_list_MR_file", type = "character", help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
corrRes_tissue_file <- args$corrRes_tissue_file
tissue_name <- args$tissue_name
pheno_name <- args$pheno_name
mrResTscore_tissue_file <- args$mrResTscore_tissue_file
mrResTscore_rev_tissue_file <- args$mrResTscore_rev_tissue_file
mrResPath_tissue_file <- args$mrResPath_tissue_file
mrResPath_rev_tissue_file <- args$mrResPath_rev_tissue_file
pheno_list_file <- args$pheno_list_file
color_pheno_file <- args$color_pheno_file
color_tissues_file <- args$color_tissues_file
pheno_list_MR_file <- args$pheno_list_MR_file
outFold <- args$outFold

#########################################################################################################################
# tissue_name <- c('DLPC_CMC','Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere','Brain_Cerebellum',
#                  'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus',
#                  'Brain_Nucleus_accumbens_basal_ganglia', 'Cells_EBV-transformed_lymphocytes')
# corrRes_tissue_file <- paste0('Meta_Analysis_SCZ/', tissue_name,'/enrichment_SCZ-UKBB_res/perc0.3_correlation_enrich_SCZ_relatedPheno.RData')
# mrResPath_tissue_file <- paste0('Meta_Analysis_SCZ/', tissue_name,'/enrichment_SCZ-UKBB_res/Mendelian_randomization_tot_path_pvalFDRrel0.05.txt')
# mrResPath_rev_tissue_file <-  paste0('Meta_Analysis_SCZ/', tissue_name,'/enrichment_SCZ-UKBB_res/Mendelian_randomization_reverse_tot_path_pvalFDRpheno0.05.txt')
# mrResTscore_tissue_file <- paste0('Meta_Analysis_SCZ/', tissue_name,'/enrichment_SCZ-UKBB_res/Mendelian_randomization_tscore_pvalFDRrel0.05.txt')
# mrResTscore_rev_tissue_file <-  paste0('Meta_Analysis_SCZ/', tissue_name,'/enrichment_SCZ-UKBB_res/Mendelian_randomization_reverse_tscore_pvalFDRpheno0.05.txt')
# outFold <- 'Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/perc0.3_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/color_pheno_type_UKBB.txt'
# pheno_name <- 'SCZ'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pheno_list_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/keep_pheno_corr_SCZ.txt'
# pheno_list_MR_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/keep_pheno_MRest_SCZ.txt'
# ######################################################################################################################

########################################################################################################################
# tissue_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# corrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_correlation_enrich_CAD_HARD_relatedPheno.RData')
# mrResPath_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_Mendelian_randomization_tot_path_pvalFDRrel0.05.txt')
# mrResPath_rev_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_Mendelian_randomization_reverse_tot_path_pvalFDRpheno0.05.txt')
# mrResTscore_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/Mendelian_randomization_tscore_pvalFDRrel0.05.txt')
# mrResTscore_rev_tissue_file <-  paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/Mendelian_randomization_reverse_tscore_pvalFDRpheno0.05.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/perc0.3_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# pheno_name <- 'CAD_HARD'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pheno_list_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_corr.txt'
# pheno_list_MR_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_MRext.txt'
######################################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissue_name,color_tissues$tissue),]
# color_tissues <- rbind(color_tissues, data.frame(tissues = 'All_tissues', color = 'black', type = 'All_GTEx', nsamples_train = NA))

pheno_keep <- read.table(pheno_list_file, h=F, stringsAsFactors = F, sep = '\t', check.names = F)$V1
pheno_plot_MR <- read.table(pheno_list_MR_file, h=F, stringsAsFactors = F, sep = '\t', check.names = F)$V1

# create matrix correlation/fisher OR
feat_corr <- lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_OR <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est <- lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est_se <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est_pval <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est_pval_FDR <- lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est_low <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mr_est_up <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mrR_est <- lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mrR_est_se <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mrR_est_pval <-lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_mrR_est_pval_FDR <- lapply(1:2, function(x) matrix(ncol = length(tissue_name), nrow = length(pheno_keep)))
feat_corr_tot <- list(NULL, NULL)
feat_corr_pval_tot <- list(NULL, NULL)
feat_mr_est_pval_tot <- list(NULL, NULL)
feat_mr_est_tot <- list(NULL, NULL)
feat_mrR_est_pval_tot <- list(NULL, NULL)
feat_mrR_est_tot <- list(NULL, NULL)

for(i in 1:length(tissue_name)){
  
  print(i)
  
  tmp <- get(load(corrRes_tissue_file[i]))
  # put NA value not significant after correction
  id_c <- match(pheno_keep, tmp$pheno$names_field)
  feat_corr[[1]][,i] <- tmp$tscore$cor_spearman[id_c]
  feat_corr_tot[[1]] <-  cbind(feat_corr_tot[[1]], tmp$tscore$cor_spearman)
  feat_corr_pval_tot[[1]] <-  cbind(feat_corr_pval_tot[[1]], tmp$tscore$cor_pval)
  feat_corr[[1]][tmp$tscore$cor_pval_BHcorr[id_c] > 0.05, i] <- NA
  feat_OR[[1]][,i] <- tmp$tscore$fisher_OR[id_c]
  feat_OR[[1]][tmp$tscore$fisher_pval_BHcorr[id_c] > 0.05, i] <- NA 
  
  feat_corr[[2]][,i] <- tmp$pathScore$cor_spearman[id_c]
  feat_corr_tot[[2]] <- cbind(feat_corr_tot[[2]], tmp$pathScore$cor_spearman)
  feat_corr_pval_tot[[2]] <-  cbind(feat_corr_pval_tot[[2]], tmp$pathScore$cor_pval)
  feat_corr[[2]][tmp$pathScore$cor_pval_BHcorr[id_c] > 0.05, i] <- NA
  feat_OR[[2]][,i] <- tmp$pathScore$fisher_OR[id_c]
  feat_OR[[2]][tmp$pathScore$fisher_pval_BHcorr[id_c] > 0.05, i] <- NA 
  
  tmp_mr <- read.delim(mrResTscore_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  id <- match(pheno_keep, tmp_mr$names_field)
  feat_mr_est[[1]][,i] <- tmp_mr$MREgg_est[id]
  feat_mr_est_tot[[1]] <- cbind(feat_mr_est_tot[[1]], tmp_mr$MREgg_est[match(tmp$pheno$names_field,tmp_mr$names_field)])
  feat_mr_est_pval_tot[[1]] <- cbind(feat_mr_est_pval_tot[[1]], tmp_mr$MREgg_est_pval[match(tmp$pheno$names_field,tmp_mr$names_field)])
  feat_mr_est_se[[1]][,i] <- tmp_mr$MREgg_est_se[id]
  feat_mr_est_pval[[1]][,i] <- tmp_mr$MREgg_est_pval[id]
  feat_mr_est_pval_FDR[[1]][,i] <- tmp_mr$MREgg_est_pval_FDRcorr[id]
  feat_mr_est_low[[1]][,i] <- tmp_mr$MREgg_est_CIl[id]
  feat_mr_est_up[[1]][,i] <- tmp_mr$MREgg_est_CIu[id]
  
  tmp_mr <- read.delim(mrResPath_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  id <- match(pheno_keep, tmp_mr$names_field)
  feat_mr_est_tot[[2]] <- cbind(feat_mr_est_tot[[2]], tmp_mr$MREgg_est[match(tmp$pheno$names_field,tmp_mr$names_field)])
  feat_mr_est_pval_tot[[2]] <- cbind(feat_mr_est_pval_tot[[2]], tmp_mr$MREgg_est_pval[match(tmp$pheno$names_field,tmp_mr$names_field)])
  feat_mr_est[[2]][,i] <- tmp_mr$MREgg_est[id]
  feat_mr_est_se[[2]][,i] <- tmp_mr$MREgg_est_se[id]
  feat_mr_est_pval[[2]][,i] <- tmp_mr$MREgg_est_pval[id]
  feat_mr_est_pval_FDR[[2]][,i] <- tmp_mr$MREgg_est_pval_FDRcorr[id]
  feat_mr_est_low[[2]][,i] <- tmp_mr$MREgg_est_CIl[id]
  feat_mr_est_up[[2]][,i] <- tmp_mr$MREgg_est_CIu[id]
  
  tmp_mrR <- read.delim(mrResTscore_rev_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  feat_mrR_est_tot[[1]] <- cbind(feat_mrR_est_tot[[1]], tmp_mrR$MREgg_est[match(tmp$pheno$names_field,tmp_mrR$names_field)])
  feat_mrR_est_pval_tot[[1]] <- cbind(feat_mrR_est_pval_tot[[1]], tmp_mrR$MREgg_est_pval[match(tmp$pheno$names_field,tmp_mrR$names_field)])
  id <- match(pheno_keep, tmp_mrR$names_field)
  feat_mrR_est[[1]][,i] <- tmp_mrR$MREgg_est[id]
  feat_mrR_est_se[[1]][,i] <- tmp_mrR$MREgg_est_se[id]
  feat_mrR_est_pval[[1]][,i] <- tmp_mrR$MREgg_est_pval[id]
  feat_mrR_est_pval_FDR[[1]][,i] <- tmp_mrR$MREgg_est_pval_FDRcorr[id]
  
  
  tmp_mrR <- read.delim(mrResPath_rev_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  feat_mrR_est_tot[[2]] <- cbind(feat_mrR_est_tot[[2]], tmp_mrR$MREgg_est[match(tmp$pheno$names_field,tmp_mrR$names_field)])
  feat_mrR_est_pval_tot[[2]] <- cbind(feat_mrR_est_pval_tot[[2]], tmp_mrR$MREgg_est_pval[match(tmp$pheno$names_field,tmp_mrR$names_field)])
  id <- match(pheno_keep, tmp_mrR$names_field)
  feat_mrR_est[[2]][,i] <- tmp_mrR$MREgg_est[id]
  feat_mrR_est_se[[2]][,i] <- tmp_mrR$MREgg_est_se[id]
  feat_mrR_est_pval[[2]][,i] <- tmp_mrR$MREgg_est_pval[id]
  feat_mrR_est_pval_FDR[[2]][,i] <- tmp_mrR$MREgg_est_pval_FDRcorr[id]
  
}

pheno_info <- tmp$pheno[id_c, ]
tot_pheno <- tmp$pheno

#### modify to save correct matrix ####
# save matrix
colnames(feat_corr_tot[[1]]) <- colnames(feat_corr_tot[[2]]) <- tissue_name
colnames(feat_corr_pval_tot[[1]]) <- colnames(feat_corr_pval_tot[[2]]) <- tissue_name
colnames(feat_mr_est_tot[[1]]) <- colnames(feat_mr_est_tot[[2]]) <- tissue_name
colnames(feat_mr_est_pval_tot[[1]]) <- colnames(feat_mr_est_pval_tot[[2]]) <- tissue_name
colnames(feat_mrR_est_tot[[1]]) <- colnames(feat_mrR_est_tot[[2]]) <- tissue_name
colnames(feat_mrR_est_pval_tot[[1]]) <- colnames(feat_mrR_est_pval_tot[[2]]) <- tissue_name

type_data <- c('tscore', 'tot_path')

for(i in 1:2){
  
  new <- cbind(feat_corr_tot[[i]], data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_correlation_allTissues.txt', outFold, type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')
  new <- cbind(feat_corr_pval_tot[[i]], data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_correlation_pvalue_allTissues.txt', outFold,  type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')

  new <- cbind(feat_mr_est_tot[[i]],  data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_MRest_allTissues.txt', outFold, type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')
  new <- cbind(feat_mr_est_pval_tot[[i]],  data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_MRest_pvalue_allTissues.txt', outFold, type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')
  
  new <- cbind(feat_mrR_est_tot[[i]],  data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_MRest_rev_allTissues.txt', outFold, type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')
  new <- cbind(feat_mrR_est_pval_tot[[i]],  data.frame(pheno = tot_pheno$names_field, pheno_type = tot_pheno$pheno_type))
  write.table(file = sprintf('%s%s_MRest_rev_pvalue_allTissues.txt', outFold, type_data[i]), x = new, col.names = T, row.names = F, quote = F, sep = '\t')
  
}

pheno_info$name_plot <- pheno_info$names_field
pheno_info$name_plot[grepl('Non-cancer illness code self-reported',pheno_info$name_plot)] <- pheno_info$pheno[grepl('Non-cancer illness code self-reported',pheno_info$name_plot)]
id_l <- which(grepl('\\(left\\)', pheno_info$name_plot))
if(length(id_l)>0){
  for(i in id_l){pheno_info$name_plot[i] <- paste0(strsplit(pheno_info$name_plot[i], split = '\\(left\\)')[[1]][1], 'left', strsplit(pheno_info$name_plot[i], split = '\\(left\\)')[[1]][2])}
}
id_r <- which(grepl('\\(right\\)', pheno_info$name_plot))
if(length(id_r)>0){
  for(i in id_r){pheno_info$name_plot[i] <- paste0(strsplit(pheno_info$name_plot[i], split = '\\(right\\)')[[1]][1], 'right', strsplit(pheno_info$name_plot[i], split = '\\(right\\)')[[1]][2])}
}
id_n <- which(grepl('(normalised for head size)', pheno_info$name_plot))
if(length(id_n)){
  for(i in id_n){pheno_info$name_plot[i] <- paste0(strsplit(pheno_info$name_plot[i], split = '\\(normalised for head size\\)')[[1]][1], 'normalised for head size', strsplit(pheno_info$name_plot[i], split = '\\(normalised for head size\\)')[[1]][2])}
}
# pheno_info$name_plot[grepl('formerly Supplementary Motor Cortex', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('formerly Supplementary Motor Cortex', pheno_info$name_plot)], function(x) paste(strsplit(x, ' \\(formerly Supplementary Motor Cortex)\\ ')[[1]], collapse = ' '))
id_a <-  which(grepl('in group-defined mask', pheno_info$name_plot))
if(length(id_a)>0){
  for(i in id_a){pheno_info$name_plot[i] <- paste0(strsplit(pheno_info$name_plot[i], split = '\\(in group-defined mask\\)')[[1]][1], 'in group-defined mask', strsplit(pheno_info$name_plot[i], split = '\\(in group-defined mask\\)')[[1]][2])}
}
id_b <-  which(grepl('in group-defined amygdala activation mask', pheno_info$name_plot))
if(length(id_b)>0){
  for(i in id_b){pheno_info$name_plot[i] <- paste0(strsplit(pheno_info$name_plot[i], split = '\\(in group-defined amygdala activation mask\\)')[[1]][1], 'in group-defined amygdala activation mask', 
                                    strsplit(pheno_info$name_plot[i], split = '\\(in group-defined amygdala activation mask\\)')[[1]][2])}
}
# pheno_info$name_plot[grepl('Heschl', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Heschl', pheno_info$name_plot)], function(x) paste(strsplit(x, ' \\((includes H1 and H2))\\ ')[[1]], collapse = ' '))
pheno_info$name_plot[grepl('Duration to complete alphanumeric path \\(trail #2\\)', pheno_info$name_plot)] <- 'Duration to complete alphanumeric path (20157)'
pheno_info$name_plot[grepl('Duration to complete numeric path \\(trail #1\\)', pheno_info$name_plot)] <- 'Duration to complete numeric path (20156)'
pheno_info$name_plot[grepl('Daytime dozing / sleeping', pheno_info$name_plot)] <- 'Daytime dozing / sleeping narcolepsy (1220)'
pheno_info$name_plot[pheno_info$name_plot == 'Diagnoses - ICD10 : I10 Essential (primary) hypertension (41270_I10)'] <- 'Diagnoses - ICD10 : I10 Essential primary hypertension (41270_I10)'
pheno_info$name_plot[grepl('Diagnoses - ICD10 : ', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Diagnoses - ', pheno_info$name_plot)], function(x) strsplit(x, split = 'Diagnoses - ')[[1]][2])
pheno_info$name_plot[grepl('Heel bone mineral density', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Heel bone mineral density', pheno_info$name_plot)], function(x) paste0(strsplit(x, split = ' \\(BMD\\)')[[1]], collapse = ''))
pheno_info$name_plot[grepl('Red blood cell', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Red blood cell', pheno_info$name_plot)], function(x) paste0(strsplit(x, split = ' \\(erythrocyte\\)')[[1]], collapse = ''))
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of father : None of the above (group 1) (20107_100)'] <- 'Illnesses of father : None of the above group 1 (20107_100)'
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of father : None of the above (group 2) (20107_101)'] <- 'Illnesses of father : None of the above group 2 (20107_101)'
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of mother : None of the above (group 1) (20110_100)'] <- 'Illnesses of mother : None of the above group 1 (20110_100)'
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of mother : None of the above (group 2) (20110_101)'] <- 'Illnesses of mother : None of the above group 2 (20110_101)'
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of siblings : None of the above (group 1) (20111_100)'] <- 'Illnesses of siblings : None of the above group 1 (20111_100)'
pheno_info$name_plot[pheno_info$name_plot == 'Illnesses of siblings : None of the above (group 2) (20111_101)'] <- 'Illnesses of siblings : None of the above group 2 (20111_101)'
pheno_info$name_plot <- sapply(pheno_info$name_plot, function(x) strsplit(x, split = ' \\(')[[1]][1])
pheno_info$name_plot[is.na(pheno_info$name_plot)] <- ''

pheno_info <- pheno_info[order(pheno_info$pheno_type),]
id <- match(pheno_info$names_field, pheno_keep)
colnames(feat_corr[[1]]) <- colnames(feat_corr[[2]]) <- tissue_name
colnames(feat_OR[[1]]) <- colnames(feat_OR[[2]]) <- tissue_name
colnames(feat_mr_est[[1]]) <- colnames(feat_mr_est_se[[1]]) <- colnames(feat_mr_est_pval[[1]]) <- tissue_name
colnames(feat_mrR_est[[1]]) <- colnames(feat_mrR_est_se[[1]]) <- colnames(feat_mrR_est_pval[[1]]) <- tissue_name
colnames(feat_mr_est_low[[1]]) <- colnames(feat_mr_est_up[[1]])  <- tissue_name
colnames(feat_mr_est[[2]]) <- colnames(feat_mr_est_se[[2]]) <- colnames(feat_mr_est_pval[[2]]) <- tissue_name
colnames(feat_mrR_est[[2]]) <- colnames(feat_mrR_est_se[[2]]) <- colnames(feat_mrR_est_pval[[2]]) <- tissue_name
colnames(feat_mr_est_low[[2]]) <- colnames(feat_mr_est_up[[2]])  <- tissue_name
colnames(feat_mr_est_pval_FDR[[1]]) <- colnames(feat_mr_est_pval_FDR[[2]])  <- tissue_name
colnames(feat_mrR_est_pval_FDR[[1]]) <- colnames(feat_mrR_est_pval_FDR[[2]])  <- tissue_name

feat_corr <- lapply(feat_corr, function(x) x[id, ])
feat_OR <- lapply(feat_OR, function(x) x[id, ])
feat_mr_est <- lapply(feat_mr_est, function(x) x[id, ])
feat_mr_est_se <- lapply(feat_mr_est_se, function(x) x[id, ])
feat_mr_est_pval <- lapply(feat_mr_est_pval, function(x) x[id, ])
feat_mr_est_low <- lapply(feat_mr_est_low, function(x) x[id, ])
feat_mr_est_up <- lapply(feat_mr_est_up, function(x) x[id, ])
feat_mrR_est <- lapply(feat_mrR_est, function(x) x[id, ])
feat_mrR_est_pval <- lapply(feat_mrR_est_pval, function(x) x[id, ])
feat_mrR_est_se <- lapply(feat_mrR_est_se, function(x) x[id, ])
feat_mr_est_pval_FDR <- lapply(feat_mr_est_pval_FDR, function(x) x[id, ])
feat_mrR_est_pval_FDR <- lapply(feat_mrR_est_pval_FDR, function(x) x[id, ])

# plot
############# make plot ################
save_pheatmap_png <- function(x, filename, width=10, height=10, res = 200) {
  png(filename, width = width, height = height, res = res, units = 'in')
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
# used to make heatmaps 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_90",
  ns = asNamespace("pheatmap")
)

plot_heatmap_corr <- function(type_mat, mat_cor, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 0.6, pheno_name, show_rownames = T){
  
  title_sub <- 'Spearman Correlation'
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  coul[50] <- '#ffffff'
  
  tmp_mat <- as.matrix(mat_cor)
  tmp_mat[is.na(tmp_mat)] <- 0
  
  val <- cap_val
  mat_breaks <- seq(-val, val, length.out = 100)
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # Data frame with column annotations.
  mat_row <- data.frame(pheno = pheno_info$pheno_type)
  rownames(mat_row) <- pheno_info$name_plot
  
  mat_colors <- list(pheno = pheno_ann$color)
  names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
  
  rownames(tmp_mat) <- pheno_info$name_plot
  
  mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat), color_tissues$tissue)])
  rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_t <- list(tissue_type = unique(color_tissues$color))
  names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
  
  new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = show_rownames, 
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s\nz-statistic %s', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 12)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_correlation_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_correlation_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}

plot_heatmap_OR <- function(type_mat, mat_OR, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 10, pheno_name, show_rownames = T){
  
  title_sub <- 'Fisher Odds Ratio'
  coul <- colorRampPalette(brewer.pal(9, "Purples"))(100)
  coul[1] <- '#ffffff'
  
  tmp_mat <- as.matrix(mat_OR)
  tmp_mat[is.na(tmp_mat)] <- 0
  
  val <- cap_val
  mat_breaks <- seq(0, val, length.out = 100)
  tmp_mat[tmp_mat>=val] <- val
  
  # Data frame with column annotations.
  mat_row <- data.frame(pheno = pheno_info$pheno_type)
  rownames(mat_row) <- pheno_info$name_plot
  
  mat_colors <- list(pheno = pheno_ann$color)
  names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
  
  rownames(tmp_mat) <- pheno_info$name_plot
  
  mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat), color_tissues$tissue)])
  rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_t <- list(tissue_type = unique(color_tissues$color))
  names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
  
  new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = show_rownames, 
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s\nz-statistic %s', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 12)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_fisherOR_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_fisherOR_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}

plot_heatmap_mr <- function(type_mat, mat_mr, mat_pval, mat_pval_FDR = NULL, pheno_ann, pheno_info, 
                            color_tissues, height_pl = 17, width_pl=11, cap_val = 2, pheno_name, show_rownames = T, title_sub = 'MR-Egger estimate'){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  coul[50] <- '#ffffff'
  
  mat_mr[mat_pval>0.05] <- NA
  labels_sign <- matrix('', nrow = nrow(mat_mr), ncol = ncol(mat_mr))
  if(!is.null(mat_pval_FDR)){
    labels_sign[mat_pval_FDR <= 0.05] <- '*'  
  }
  
  tmp_mat <- as.matrix(mat_mr)
  tmp_mat[is.na(tmp_mat)] <- 0
  
  val <- cap_val
  mat_breaks <- seq(-val, val, length.out = 100)
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # Data frame with column annotations.
  mat_row <- data.frame(pheno = pheno_info$pheno_type)
  rownames(mat_row) <- pheno_info$name_plot
  
  mat_colors <- list(pheno = pheno_ann$color)
  names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
  
  rownames(tmp_mat) <- pheno_info$name_plot
  
  mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat), color_tissues$tissue)])
  rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_t <- list(tissue_type = unique(color_tissues$color))
  names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
  
  new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = show_rownames, 
                    display_numbers = labels_sign, fontsize_number = 12, number_color = 'black',
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s\n%s', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 12)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}

type_data <- c('tscore', 'tot_path')
id_MR <- which(pheno_info$names_field %in% pheno_plot_MR)

for(i in 1:2){
  
  # correlation
  plot_heatmap_corr(type_mat = type_data[i], mat_cor = feat_corr[[i]][,tissue_name], pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                    pheno_info = pheno_info[, ],  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T, 
                    cap_val = min(c(0.6, max(abs(feat_corr[[i]][,tissue_name]), na.rm = T))))

  # Odds ratio
  plot_heatmap_OR(type_mat = type_data[i], mat_OR = feat_OR[[i]][,tissue_name], pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                  pheno_info = pheno_info,  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T, 
                  cap_val = min(c(10, max(abs(feat_OR[[i]][,tissue_name]), na.rm = T))))

  # MR-Egger
  plot_heatmap_mr(type_mat = type_data[i], mat_mr = feat_mr_est[[i]][,tissue_name], mat_pval = feat_mr_est_pval[[i]][,tissue_name],
                  pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], cap_val = 1, 
                  pheno_info = pheno_info[,],  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T)
  # MR-Egger
  plot_heatmap_mr(type_mat = type_data[i], mat_mr = feat_mrR_est[[i]][,tissue_name], mat_pval = feat_mrR_est_pval[[i]][,tissue_name],
                  pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], cap_val = min(c(1, max(abs(feat_mrR_est[[i]][,tissue_name]), na.rm = T))), 
                  pheno_info = pheno_info[,],  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = paste0('reverse_',pheno_name), show_rownames = T, title_sub = 'MR-Egger reverse: estimate')
  
  # MR-Egger (pvalue) and correlation selected phenotypes
  dat <- sign(feat_mr_est[[i]][id_MR,tissue_name])*(-log10(feat_mr_est_pval[[i]][id_MR, tissue_name]))
  id_inf <- dat == Inf & !is.na(dat)
  if(any(id_inf)){
    dat[id_inf] <- -log10(2*pnorm(-abs(feat_mr_est[[i]][id_MR,tissue_name][id_inf]/feat_mr_est_se[[i]][id_MR,tissue_name][id_inf])))*sign(feat_mr_est[[i]][id_MR,tissue_name][id_inf])
  }
  
  plot_heatmap_mr(type_mat = type_data[i], mat_mr = dat, mat_pval = feat_mr_est_pval[[i]][id_MR,tissue_name], mat_pval_FDR = feat_mr_est_pval_FDR[[i]][id_MR,tissue_name], 
                  pheno_ann = color_pheno[match(unique(pheno_info[id_MR,]$pheno_type),color_pheno$pheno_type),], cap_val = min(c(10, max(abs(dat), na.rm = T))), 
                  pheno_info = pheno_info[id_MR,],  color_tissues = color_tissues, height_pl = 9, width_pl = 14, pheno_name = paste0('specPheno',pheno_name), show_rownames = T, 
                  title_sub = 'MR-Egger: signed -log10(pvalue)')
  
  dat <- sign(feat_mrR_est[[i]][id_MR,tissue_name])*(-log10(feat_mrR_est_pval[[i]][id_MR, tissue_name]))
  id_inf <- dat == Inf & !is.na(dat)
  if(any(id_inf)){
    dat[id_inf] <- -log10(2*pnorm(-abs(feat_mrR_est[[i]][id_MR,tissue_name][id_inf]/feat_mrR_est_se[[i]][id_MR,tissue_name][id_inf])))*sign(feat_mrR_est[[i]][id_MR,tissue_name][id_inf])
  }
  
  plot_heatmap_mr(type_mat = type_data[i], mat_mr = dat, mat_pval = feat_mrR_est_pval[[i]][id_MR,tissue_name], mat_pval_FDR = feat_mrR_est_pval_FDR[[i]][id_MR,tissue_name], 
                  pheno_ann = color_pheno[match(unique(pheno_info[id_MR,]$pheno_type),color_pheno$pheno_type),], cap_val = min(c(10, max(abs(dat), na.rm = T))), 
                  pheno_info = pheno_info[id_MR,],  color_tissues = color_tissues, height_pl = 9, width_pl = 14, pheno_name = paste0('reverse_specPheno',pheno_name), show_rownames = T, 
                  title_sub = 'MR-Egger reverse: signed -log10(pvalue)')
  
  plot_heatmap_corr(type_mat = type_data[i], mat_cor = feat_corr[[i]][id_MR,tissue_name], pheno_ann = color_pheno[match(unique(pheno_info[id_MR,]$pheno_type),color_pheno$pheno_type),], 
                    pheno_info = pheno_info[id_MR, ],  color_tissues = color_tissues, height_pl = 9, width_pl = 14, pheno_name = paste0('specPheno',pheno_name), show_rownames = T, 
                    cap_val = min(c(0.6, max(abs(feat_corr[[i]][id_MR,tissue_name]), na.rm = T))))
  
}

# plot correlation path and tscore: 
df_tot <- data.frame(Tscore = as.vector(feat_corr_tot[[1]]), pathScore = as.vector(feat_corr_tot[[2]]), 
                     tissue = unlist(lapply(tissue_name, function(x) rep(x, nrow(feat_corr_tot[[1]])))),
                     stringsAsFactors = F)
df_tot$tissue <- factor(df_tot$tissue, levels = tissue_name)
r2_val <- data.frame(tissue = tissue_name, r = sapply(tissue_name, function(x) cor(df_tot$Tscore[df_tot$tissue == x], df_tot$pathScore[df_tot$tissue == x])))
r2_val$tissue <- factor(r2_val$tissue, levels = tissue_name)
r2_val$r_lab <- paste('R = ', as.character(round(r2_val$r, digits = 3)))

file_name <- sprintf('%scorrelation_%s_relatedPheno_compareTissueSpec_tscore_path', outFold, pheno_name)
pl <- ggplot(data = df_tot, aes(x = Tscore, y = pathScore, color = tissue))+
  geom_point(alpha = 0.6, size = 0.8)+
  facet_wrap(.~tissue, nrow = 2)+
  # geom_text_repel(segment.color = 'grey50', color = 'black', size = 2, min.segment.length = unit(0, 'lines'),
  #                 segment.alpha = 0.6,  force = 15) +
  geom_text(data = r2_val, aes(x = 0, y = ifelse(grepl('CAD', pheno_name), 0.8, 0.3), label = r_lab), 
            color = 'black')+
  xlab('gene T-score')+ylab('Pathway score')+ 
  theme_bw()+ ggtitle('Spearman Correlation Z-statistic')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), 
                                                       legend.text = element_text(size = 6),  legend.key.size = unit(0.3, "cm"), legend.title = element_blank())+
  guides(color=guide_legend(ncol=1))+
  geom_abline(slope = 1, intercept= 0, linetype = 'dashed', alpha = 0.6, color = 'black')+
  # geom_smooth(method=lm, se = FALSE, size = 0.5, linetype = 'dashed', alpha = 0.6, color = 'black')+
  scale_color_manual(values = color_tissues$color)

ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 11, height = 4.5, dpi = 500)
ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 11, height = 4.5, dpi = 500, compress = F)
write.table(r2_val,file = paste0(file_name, '.txt') , col.names = T, row.names = F, sep = '\t', quote = F)

# plot number of phenotypes that have MR estimate sign
df_perc <- data.frame(nMR = colSums(feat_mr_est_pval_tot[[1]] < 0.05, na.rm = T), nCorr = colSums(!is.na(feat_mr_est_pval_tot[[1]])),
                     perc = colSums(feat_mr_est_pval_tot[[1]] < 0.05, na.rm = T)/colSums(!is.na(feat_mr_est_pval_tot[[1]])), 
                     type = rep('T-score', length(tissue_name)), tissue = tissue_name, stringsAsFactors = F)
df_perc <- rbind(df_perc, data.frame(nMR = colSums(feat_mr_est_pval_tot[[2]] < 0.05, na.rm = T), nCorr = colSums(!is.na(feat_mr_est_pval_tot[[2]])), 
                                     perc = colSums(feat_mr_est_pval_tot[[2]] < 0.05, na.rm = T)/colSums(!is.na(feat_mr_est_pval_tot[[2]])), 
                                     type = rep('Pathway score', length(tissue_name)), tissue = tissue_name, stringsAsFactors = F))
rownames(df_perc) <- NULL
df_perc$type <- factor(df_perc$type, levels = c('T-score', 'Pathway score'))
df_perc$tissue <- factor(df_perc$tissue, levels = tissue_name)

file_name <- sprintf('%sMRestCount_%s_relatedPheno_compareTissueSpec_tscore_path', outFold, pheno_name)
pl <- ggplot(data = df_perc, aes(x = tissue, y = nMR, fill = type))+
  geom_bar(alpha = 0.6, width = 0.8, color = 'black', stat = 'identity', position="dodge")+
  ylab('Number of significant MR estimates')+ 
  theme_bw()+ theme(axis.text.y = element_text(color = color_tissues$color), axis.title.y = element_blank(), legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
                     legend.title = element_blank())+
  scale_fill_manual(values = c('lightgrey', 'grey20'))+
  coord_flip()

ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 6, height = 4.5, dpi = 500)
ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 6, height = 4.5, dpi = 500, compress = F)

write.table(df_perc,file = paste0(file_name, '.txt') , col.names = T, row.names = F, sep = '\t', quote = F)

### reverse
df_perc <- data.frame(nMR = colSums(feat_mrR_est_pval_tot[[1]] < 0.05, na.rm = T), nCorr = colSums(!is.na(feat_mrR_est_pval_tot[[1]])),
                      perc = colSums(feat_mrR_est_pval_tot[[1]] < 0.05, na.rm = T)/colSums(!is.na(feat_mrR_est_pval_tot[[1]])), 
                      type = rep('T-score', length(tissue_name)), tissue = tissue_name, stringsAsFactors = F)
df_perc <- rbind(df_perc, data.frame(nMR = colSums(feat_mrR_est_pval_tot[[2]] < 0.05, na.rm = T), nCorr = colSums(!is.na(feat_mrR_est_pval_tot[[2]])), 
                                     perc = colSums(feat_mrR_est_pval_tot[[2]] < 0.05, na.rm = T)/colSums(!is.na(feat_mrR_est_pval_tot[[2]])), 
                                     type = rep('Pathway score', length(tissue_name)), tissue = tissue_name, stringsAsFactors = F))
rownames(df_perc) <- NULL
df_perc$type <- factor(df_perc$type, levels = c('T-score', 'Pathway score'))
df_perc$tissue <- factor(df_perc$tissue, levels = tissue_name)

file_name <- sprintf('%sMRestCount_reverse_%s_relatedPheno_compareTissueSpec_tscore_path', outFold, pheno_name)
pl <- ggplot(data = df_perc, aes(x = tissue, y = nMR, fill = type))+
  geom_bar(alpha = 0.6, width = 0.8, color = 'black', stat = 'identity', position="dodge")+
  ylab('Number of significant MR estimates')+ 
  theme_bw()+ theme(axis.text.y = element_text(color = color_tissues$color), axis.title.y = element_blank(), legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
                    legend.title = element_blank())+
  scale_fill_manual(values = c('lightgrey', 'grey20'))+
  coord_flip()

ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 6, height = 4.5, dpi = 500)
ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 6, height = 4.5, dpi = 500, compress = F)

write.table(df_perc,file = paste0(file_name, '.txt') , col.names = T, row.names = F, sep = '\t', quote = F)


# df_perc = data.frame(nMR = rowSums(feat_mr_est_pval[[1]] < 0.05, na.rm = T), nCorr = rowSums(!is.na(feat_mr_est_pval[[1]])),
#                      perc = rowSums(feat_mr_est_pval[[1]] < 0.05, na.rm = T)/rowSums(!is.na(feat_mr_est_pval[[1]])), 
#                      type = rep('tscore', length(pheno_keep)), pheno = pheno_keep, stringsAsFactors = F)
# df_perc <- rbind(df_perc, data.frame(nMR = rowSums(feat_mr_est_pval[[2]] < 0.05, na.rm = T), nCorr = rowSums(!is.na(feat_mr_est_pval[[2]])), 
#                                      perc = rowSums(feat_mr_est_pval[[2]] < 0.05, na.rm = T)/rowSums(!is.na(feat_mr_est_pval[[2]])), 
#                                      type = rep('tot_path', length(pheno_keep)), pheno = pheno_keep, stringsAsFactors = F))
# rownames(df_perc) <- NULL

#### plot specific MR results with pvalue and SE ###
id_MR <- which(pheno_info$names_field %in% pheno_plot_MR)
id_keep <- 2
df_mr <- data.frame(MR_est = as.vector(feat_mr_est[[id_keep]][id_MR, ]), MR_est_CIl = as.vector(feat_mr_est_low[[id_keep]][id_MR, ]), 
                    MR_est_CIu = as.vector(feat_mr_est_up[[id_keep]][id_MR, ]), MR_est_pval = as.vector(feat_mr_est_pval[[id_keep]][id_MR, ]), MR_est_se = as.vector(feat_mr_est_se[[id_keep]][id_MR, ]), 
                    pheno_name = rep(pheno_info$name_plot[id_MR], ncol(feat_mr_est[[id_keep]])), tissue = unlist(lapply(colnames(feat_mr_est[[id_keep]]), function(x) rep(x, length(id_MR)))), stringsAsFactors = F)
df_mr$sign <- 'no'
df_mr$sign[df_mr$MR_est_pval <= 0.05] <- 'yes'
df_mr$sign <- factor(df_mr$sign, levels = c('no', 'yes'))
df_mr$pheno_name <- factor(df_mr$pheno_name, levels = pheno_info$name_plot[id_MR])
df_mr$tissue <- factor(df_mr$tissue, levels =color_tissues$tissues)
df_mr$log_pval <- -log10(df_mr$MR_est_pval)
id_inf <- df_mr$log_pval == Inf & !is.na(df_mr$log_pval)
if(any(id_inf)){
  df_mr$log_pval[id_inf] <- -log10(2*pnorm(-abs(df_mr$MR_est[id_inf]/df_mr$MR_est_se[id_inf])))
}


if(grepl('CAD', pheno_name)){
  myPalette <- colorRampPalette(c('grey20', rev(brewer.pal(11, "Spectral"))))
  pl <-  ggplot(subset(df_mr, tissue %in% c('Adipose_Subcutaneous', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')),
              aes(x = pheno_name, y = MR_est, shape = sign, color = log_pval))+
  # pl <-  ggplot(df_mr,aes(x = pheno_name, y = MR_est, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=MR_est_CIl, ymax=MR_est_CIu), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    facet_wrap(.~tissue, nrow = 1)+
    ylab('MR-Egger estimate') + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    theme(legend.position = 'bottom', legend.title = element_text(size = 7), legend.text = element_text(size = 7), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 8, angle = 0, hjust = 1), axis.text.y = element_text(size = 7.5), 
          strip.text = element_text(size=6.7))+
    scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(df_mr$log_pval, na.rm = T)))+
    scale_shape_manual(values=c(0, 15))+
    guides(shape=FALSE, color = guide_colourbar(barwidth=10,label.position="bottom", barheight = 1))+labs(color = '-log10(pvalue)')+coord_flip()
  # scale_color_manual(values=color_tissues$color[color_tissues$tissue %in% c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')])

  pl <- ggplot_gtable(ggplot_build(pl))
  stripr <- which(grepl('strip-t', pl$layout$name))
  fills <- color_tissues$color[color_tissues$tissue %in% c('Adipose_Subcutaneous',  'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')]
  # fills <- color_tissues$color
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', pl$grobs[[i]]$grobs[[1]]$childrenOrder))
    pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$alpha <- 0.7
    k <- k+1
  }

  ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.png", outFold, type_data[id_keep], pheno_name), width = 11, height = 4, plot = pl, device = 'png')
  ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.pdf", outFold, type_data[id_keep], pheno_name), width = 11, height = 4, plot = pl, device = 'pdf')
}


# myPalette <- colorRampPalette(c('grey20', 'red'))
# if(grepl('SCZ', pheno_name)){
#   
#   df_mr$MR_est_CIl[df_mr$MR_est_CIl <= -2000] <- -2000
#   df_mr$MR_est_CIu[df_mr$MR_est_CIu >= 2000] <- 2000
#   pl <-  ggplot(subset(df_mr, tissue %in% c('DLPC_CMC', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
#                                             'Brain_Cerebellum', 'Cells_EBV-transformed_lymphocytes')),
#                 aes(x = pheno_name, y = MR_est, shape = sign, color = log_pval))+
#     # pl <-  ggplot(df_mr,aes(x = pheno_name, y = MR_est, shape = sign))+
#     geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=MR_est_CIl, ymax=MR_est_CIu), width=.2, position=position_dodge(0.5))+
#     theme_bw()+ 
#     facet_wrap(.~tissue, nrow = 1)+
#     ylab('MR-Egger estimate') + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
#     theme(legend.position = 'bottom', legend.title = element_text(size = 7), legend.text = element_text(size = 7), 
#           plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
#           axis.text.x = element_text(size = 8, angle = 45, hjust = 1), axis.text.y = element_text(size = 7.5), 
#           strip.text = element_text(size=6.7))+
#     scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(df_mr$log_pval, na.rm = T)))+
#     scale_shape_manual(values=c(0, 15))+
#     scale_y_continuous(limits = c(-2000, 2000), labels = function(x) format(x, scientific = TRUE))+
#     guides(shape=FALSE, color = guide_colourbar(barwidth=10,label.position="bottom", barheight = 1))+labs(color = '-log10(pvalue)')+coord_flip()
#   # scale_color_manual(values=color_tissues$color[color_tissues$tissue %in% c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')])
#   
#   pl <- ggplot_gtable(ggplot_build(pl))
#   stripr <- which(grepl('strip-t', pl$layout$name))
#   fills <- color_tissues$color[color_tissues$tissue %in% c('DLPC_CMC', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
#                                                            'Brain_Cerebellum', 'Cells_EBV-transformed_lymphocytes')]
#   # fills <- color_tissues$color
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', pl$grobs[[i]]$grobs[[1]]$childrenOrder))
#     pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#     pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$alpha <- 0.7
#     k <- k+1
#   }
#   
#   ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.png", outFold, type_data[id_keep], pheno_name), width = 11, height = 5, plot = pl, device = 'png')
#   ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.pdf", outFold, type_data[id_keep], pheno_name), width = 11, height = 5, plot = pl, device = 'pdf')
#   
#   
# }



