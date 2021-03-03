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
suppressPackageStartupMessages(library(corrplot))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="correlation and MR plots: tissue specific")
parser$add_argument("--corrRes_file", type = "character", help = "")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--mrRes_Egg_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_pval_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_pval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_intpval_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_intpval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_pval_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_pval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_FDRpval_file", type = "character", help = "")
parser$add_argument("--mrRes_IVW_FDRpval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_intFDRpval_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_intFDRpval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_FDRpval_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_FDRpval_rev_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_int_file", type = "character", help = "")
parser$add_argument("--mrRes_Egg_int_rev_file", type = "character", help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_list_MR_file", type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--data_type", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")


args <- parser$parse_args()
corrRes_file <- args$corrRes_file
tissue_name <- args$tissue_name
pheno_name <- args$pheno_name
mrRes_Egg_file <- args$mrRes_Egg_file
mrRes_Egg_rev_file <- args$mrRes_Egg_rev_file
mrRes_Egg_pval_file <- args$mrRes_Egg_pval_file
mrRes_Egg_pval_rev_file <- args$mrRes_Egg_pval_rev_file
mrRes_Egg_intpval_file <- args$mrRes_Egg_intpval_file
mrRes_Egg_intpval_rev_file <- args$mrRes_Egg_intpval_rev_file
mrRes_IVW_file <- args$mrRes_IVW_file
mrRes_IVW_pval_file <- args$mrRes_IVW_pval_file
mrRes_IVW_rev_file <- args$mrRes_IVW_rev_file
mrRes_IVW_pval_rev_file <- args$mrRes_IVW_pval_rev_file

mrRes_IVW_FDRpval_file <- args$mrRes_IVW_FDRpval_file
mrRes_IVW_FDRpval_rev_file <- args$mrRes_IVW_FDRpval_rev_file
mrRes_Egg_intFDRpval_file <- args$mrRes_Egg_intFDRpval_file
mrRes_Egg_intFDRpval_rev_file <- args$mrRes_Egg_intFDRpval_rev_file
mrRes_Egg_FDRpval_file <- args$mrRes_Egg_FDRpval_file
mrRes_Egg_FDRpval_rev_file <- args$mrRes_Egg_FDRpval_rev_file
mrRes_Egg_int_file <- args$mrRes_Egg_int_file
mrRes_Egg_int_rev_file <- args$mrRes_Egg_int_rev_file

color_pheno_file <- args$color_pheno_file
color_tissues_file <- args$color_tissues_file
pheno_list_MR_file <- args$pheno_list_MR_file
data_type <- args$data_type
outFold <- args$outFold

#########################################################################################################################
# tissue_name <- c('DLPC_CMC','Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere','Brain_Cerebellum',
#                  'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus',
#                  'Brain_Nucleus_accumbens_basal_ganglia', 'Cells_EBV-transformed_lymphocytes')
# corrRes_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt')
# mrRes_Egg_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_Egger_estimates_allTissues.txt')
# mrRes_Egg_rev_file <-  paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_allTissues.txt')
# mrRes_Egg_pval_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_Egger_estimates_pvalue_allTissues.txt')
# mrRes_Egg_pval_rev_file <-  paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_pvalue_allTissues.txt')
# mrRes_Egg_intpval_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_pvalue_allTissues.txt')
# mrRes_Egg_intpval_rev_file <-  paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_pvalue_allTissues.txt')
# 
# mrRes_IVW_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_IVW_estimates_allTissues.txt')
# mrRes_IVW_pval_file <-  paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_IVW_estimates_pvalue_allTissues.txt')
# mrRes_IVW_rev_file <- paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_allTissues.txt')
# mrRes_IVW_pval_rev_file <-  paste0('Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_pvalue_allTissues.txt')
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# outFold <- 'Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/dist200000b_'
# pheno_name <- 'SCZ'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/color_pheno_type_UKBB.txt'
# pheno_list_MR_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/keep_pheno_MRest_SCZ_new.txt'
# data_type <- 'tscore'
# #####################################################################################################################

########################################################################################################################
# tissue_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# corrRes_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt')
# mrRes_Egg_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_estimates_allTissues.txt')
# mrRes_Egg_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_allTissues.txt')
# mrRes_Egg_pval_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_estimates_pvalue_allTissues.txt')
# mrRes_Egg_FDRpval_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_estimates_FDRpvalue_allTissues.txt')
# mrRes_Egg_FDRpval_rev_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_FDRpvalue_allTissues.txt')
# mrRes_Egg_pval_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_pvalue_allTissues.txt')
# mrRes_Egg_intpval_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_pvalue_allTissues.txt')
# mrRes_Egg_intpval_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_pvalue_allTissues.txt')
# mrRes_Egg_int_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_allTissues.txt')
# mrRes_Egg_int_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_allTissues.txt')
# mrRes_Egg_intFDRpval_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_FDRpvalue_allTissues.txt')
# mrRes_Egg_intFDRpval_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_FDRpvalue_allTissues.txt')
# mrRes_IVW_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_IVW_estimates_allTissues.txt')
# mrRes_IVW_pval_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_IVW_estimates_pvalue_allTissues.txt')
# mrRes_IVW_FDRpval_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_IVW_estimates_FDRpvalue_allTissues.txt')
# mrRes_IVW_rev_file <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_allTissues.txt')
# mrRes_IVW_pval_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_pvalue_allTissues.txt')
# mrRes_IVW_FDRpval_rev_file <-  paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB//enrichment_CADHARD_res/CADCardioG_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_FDRpvalue_allTissues.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/CADCardioG_dist200000b_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# data_type <- 'tscore'
# pheno_name <- 'CAD'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pheno_list_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_corr.txt'
# pheno_list_MR_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_MRext.txt'
# ######################################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissue_name,color_tissues$tissue),]

pheno_plot_MR <- read.table(pheno_list_MR_file, h=F, stringsAsFactors = F, sep = '\t', check.names = F)$V1

# load matrices
corr_mat <- read.delim(corrRes_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)

mrRes_Egg <- read.delim(mrRes_Egg_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_rev <- read.delim(mrRes_Egg_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_pval <- read.delim(mrRes_Egg_pval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_pval_rev <- read.delim(mrRes_Egg_pval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_FDRpval <- read.delim(mrRes_Egg_FDRpval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_FDRpval_rev <- read.delim(mrRes_Egg_FDRpval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)

mrRes_Egg_int <- read.delim(mrRes_Egg_int_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_int_rev <- read.delim(mrRes_Egg_int_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_int_pval <- read.delim(mrRes_Egg_intpval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_int_pval_rev <- read.delim(mrRes_Egg_intpval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_int_FDRpval <- read.delim(mrRes_Egg_intFDRpval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_Egg_int_FDRpval_rev <- read.delim(mrRes_Egg_intFDRpval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)

mrRes_IVW <- read.delim(mrRes_IVW_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_IVW_rev <- read.delim(mrRes_IVW_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_IVW_pval <- read.delim(mrRes_IVW_pval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_IVW_pval_rev <- read.delim(mrRes_IVW_pval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_IVW_FDRpval <- read.delim(mrRes_IVW_FDRpval_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)
mrRes_IVW_FDRpval_rev <- read.delim(mrRes_IVW_FDRpval_rev_file, h=T, stringsAsFactors = F, sep = '\t', check.names = F)

pheno_ann <- mrRes_Egg[, 1:3]

# build signed pvalue matrices
mrRes_Egg_sign <- sign(as.matrix(mrRes_Egg[, -c(1:3)]))*-log10(as.matrix(mrRes_Egg_pval[, -c(1:3)]))
# mrRes_Egg_sign[mrRes_Egg_pval[, -c(1:3)] > 0.05 & !is.na(mrRes_Egg_pval[, -c(1:3)])] <- NA
mrRes_Egg_sign <- as.data.frame(mrRes_Egg_sign)

mrRes_Egg_int_sign <- sign(as.matrix(mrRes_Egg_int[, -c(1:3)]))*-log10(as.matrix(mrRes_Egg_int_pval[, -c(1:3)]))
# mrRes_Egg_sign[mrRes_Egg_pval[, -c(1:3)] > 0.05 & !is.na(mrRes_Egg_pval[, -c(1:3)])] <- NA
mrRes_Egg_int_sign <- as.data.frame(mrRes_Egg_int_sign)

mrRes_Egg_rev_sign <- sign(as.matrix(mrRes_Egg_rev[, -c(1:3)]))*-log10(as.matrix(mrRes_Egg_pval_rev[, -c(1:3)]))
# mrRes_Egg_rev_sign[mrRes_Egg_pval_rev[, -c(1:3)] > 0.05 & !is.na(mrRes_Egg_pval_rev[, -c(1:3)])] <- NA
mrRes_Egg_rev_sign <- as.data.frame(mrRes_Egg_rev_sign)

mrRes_Egg_int_rev_sign <- sign(as.matrix(mrRes_Egg_int_rev[, -c(1:3)]))*-log10(as.matrix(mrRes_Egg_int_pval_rev[, -c(1:3)]))
# mrRes_Egg_sign[mrRes_Egg_pval[, -c(1:3)] > 0.05 & !is.na(mrRes_Egg_pval[, -c(1:3)])] <- NA
mrRes_Egg_int_rev_sign <- as.data.frame(mrRes_Egg_int_rev_sign)

mrRes_IVW_sign <- sign(as.matrix(mrRes_IVW[, -c(1:3)]))*-log10(as.matrix(mrRes_IVW_pval[, -c(1:3)]))
# mrRes_IVW_sign[mrRes_IVW_pval[, -c(1:3)] > 0.05 & !is.na(mrRes_IVW_pval[, -c(1:3)])] <- NA
mrRes_IVW_sign <- as.data.frame(mrRes_IVW_sign)

mrRes_IVW_rev_sign <- sign(as.matrix(mrRes_IVW_rev[, -c(1:3)]))*-log10(as.matrix(mrRes_IVW_pval_rev[, -c(1:3)]))
# mrRes_IVW_rev_sign[mrRes_IVW_pval_rev[, -c(1:3)] > 0.05 & !is.na(mrRes_IVW_pval_rev[, -c(1:3)])] <- NA
mrRes_IVW_rev_sign <- as.data.frame(mrRes_IVW_rev_sign)

## filter based on chosen phenotypes
id_keep <- which(pheno_ann$names_field %in% pheno_plot_MR)
corr_mat <- corr_mat[corr_mat$names_field %in% pheno_plot_MR, ]
mrRes_Egg_sign <- mrRes_Egg_sign[id_keep, ]
mrRes_Egg_rev_sign <- mrRes_Egg_rev_sign[id_keep, ]
mrRes_Egg_FDRpval <- mrRes_Egg_FDRpval[id_keep, -c(1:3)]
mrRes_Egg_FDRpval_rev <- mrRes_Egg_FDRpval_rev[id_keep, -c(1:3)]
mrRes_Egg_int_pval <- mrRes_Egg_int_pval[id_keep, -c(1:3)]
mrRes_Egg_int_pval_rev <- mrRes_Egg_int_pval_rev[id_keep, -c(1:3)]
mrRes_Egg_int_FDRpval <- mrRes_Egg_int_FDRpval[id_keep, -c(1:3)]
mrRes_Egg_int_FDRpval_rev <- mrRes_Egg_int_FDRpval_rev[id_keep, -c(1:3)]
mrRes_Egg_int_sign <- mrRes_Egg_int_sign[id_keep, ]
mrRes_Egg_int_rev_sign <- mrRes_Egg_int_rev_sign[id_keep, ]

mrRes_IVW_sign <- mrRes_IVW_sign[id_keep,]
mrRes_IVW_rev_sign <- mrRes_IVW_rev_sign[id_keep,]
mrRes_IVW_FDRpval <- mrRes_IVW_FDRpval[id_keep, -c(1:3)]
mrRes_IVW_FDRpval_rev <- mrRes_IVW_FDRpval_rev[id_keep, -c(1:3)]

pheno_info <- pheno_ann[id_keep, ]

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

# reorder
mrRes_Egg_rev_sign <- mrRes_Egg_rev_sign[order(pheno_info$pheno_type), ]
mrRes_Egg_sign <- mrRes_Egg_sign[order(pheno_info$pheno_type), ]
mrRes_Egg_FDRpval_rev <- mrRes_Egg_FDRpval_rev[order(pheno_info$pheno_type), ]
mrRes_Egg_FDRpval <- mrRes_Egg_FDRpval[order(pheno_info$pheno_type), ]

mrRes_Egg_int_rev_sign <- mrRes_Egg_int_rev_sign[order(pheno_info$pheno_type), ]
mrRes_Egg_int_sign <- mrRes_Egg_int_sign[order(pheno_info$pheno_type), ]
mrRes_Egg_int_pval_rev <- mrRes_Egg_int_pval_rev[order(pheno_info$pheno_type), ]
mrRes_Egg_int_pval <- mrRes_Egg_int_pval[order(pheno_info$pheno_type), ]
mrRes_Egg_int_FDRpval_rev <- mrRes_Egg_int_FDRpval_rev[order(pheno_info$pheno_type), ]
mrRes_Egg_int_FDRpval <- mrRes_Egg_int_FDRpval[order(pheno_info$pheno_type), ]

mrRes_IVW_rev_sign <- mrRes_IVW_rev_sign[order(pheno_info$pheno_type), ]
mrRes_IVW_sign <- mrRes_IVW_sign[order(pheno_info$pheno_type), ]
mrRes_IVW_FDRpval <- mrRes_IVW_FDRpval[order(pheno_info$pheno_type), ]
mrRes_IVW_FDRpval_rev <- mrRes_IVW_FDRpval_rev[order(pheno_info$pheno_type), ]

pheno_info <- pheno_info[order(pheno_info$pheno_type),]

#### create matrix that put Egger estimate if interept p-value <= 0.05 ###
# mat_sign_rev <- mrRes_IVW_rev_sign
# id <- mrRes_Egg_int_pval_rev <= 0.05 & !is.na(mrRes_Egg_int_pval_rev) 
# mat_sign_rev[id] <- mrRes_Egg_rev_sign[id]
# pvalue_dummy_mat_rev <- matrix(NA, nrow = nrow(mat_sign_rev), ncol = ncol(mat_sign_rev))
# pvalue_dummy_mat_rev[!is.na(mat_sign_rev)] <- ' '
# pvalue_dummy_mat_rev[id] <- '•'
# rownames(mat_sign_rev) <- pheno_info$name_plot 
# mat_sign_rev <- as.matrix(mat_sign_rev)
# 
# mat_sign <- mrRes_IVW_sign
# id <- mrRes_Egg_int_pval <= 0.05 & !is.na(mrRes_Egg_int_pval) 
# mat_sign[id] <- mrRes_Egg_sign[id]
# pvalue_dummy_mat <- matrix(NA, nrow = nrow(mat_sign), ncol = ncol(mat_sign))
# pvalue_dummy_mat[!is.na(mat_sign)] <- ' '
# pvalue_dummy_mat[id] <- '•'
# rownames(mat_sign) <- pheno_info$name_plot 
# mat_sign <- as.matrix(mat_sign)

corr_mat_new <- corr_mat[match(pheno_info$pheno, corr_mat$pheno), tissue_name]
# png('provat.png', width = 1500, height = 1500)
# corrplot(mat_sign, is.corr = FALSE,  method = "color", na.label = "square",
#          na.label.col = "grey80",
#          p.mat = pvalue_dummy_mat, sig.level = 0.05)
# dev.off()


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

plot_heatmap_mr <- function(type_mat, mat_mr, mat_corr, pheno_ann, pheno_info, 
                            color_tissues, height_pl = 17, width_pl=11, cap_val = 2, pheno_name, show_rownames = T, title_sub = 'MR-Egger estimate'){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  coul[50] <- '#ffffff'
  
  # mat_mr[mat_pval>0.1] <- 0
  labels_sign <- matrix('', nrow = nrow(mat_mr), ncol = ncol(mat_mr))
  if(!is.null(mat_corr)){
    labels_sign[mat_corr <= 0.05] <- '*'
  }
  # labels_sign <- mat_type
  # labels_sign[is.na(labels_sign)] <- ''
  
  tmp_mat <- as.matrix(mat_mr)
  tmp_mat[!is.na(tmp_mat) & abs(tmp_mat) < -log10(0.05)] <- 0
  
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
                    main = sprintf('%s\n%s', title_sub, type_mat), na_col = "grey90",
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 12)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}

# correlation
plot_heatmap_corr(type_mat = data_type, mat_cor = corr_mat_new[,tissue_name], 
                  pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                  pheno_info = pheno_info, color_tissues = color_tissues, 
                  height_pl = 9, width_pl = 14, pheno_name = pheno_name, show_rownames = T, 
                  cap_val = min(c(0.6, max(abs(corr_mat_new[,tissue_name]), na.rm = T))))

plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_IVW_sign, mat_corr = mrRes_IVW_FDRpval, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_IVW_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('IVW_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-IVW: signed -log10(pvalue)')

plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_IVW_rev_sign, mat_corr = mrRes_IVW_FDRpval_rev, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_IVW_rev_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('reverse_IVW_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-IVW reverse: signed -log10(pvalue)')

plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_Egg_sign, mat_corr = mrRes_Egg_FDRpval, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_Egg_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('Egger_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-Egger: signed -log10(pvalue)')

plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_Egg_rev_sign, mat_corr = mrRes_Egg_FDRpval_rev, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_Egg_rev_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('reverse_Egger_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-Egger reverse: signed -log10(pvalue)')

# intercept MR egger
plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_Egg_int_sign, mat_corr = mrRes_Egg_int_FDRpval, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_Egg_int_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('EggerInt_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-Egger intercept: signed -log10(pvalue)')

plot_heatmap_mr(type_mat = data_type, mat_mr = mrRes_Egg_int_rev_sign, mat_corr = mrRes_Egg_int_FDRpval_rev, 
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                cap_val = min(c(10, max(abs(mrRes_Egg_int_rev_sign), na.rm = T))), 
                pheno_info = pheno_info,  color_tissues = color_tissues, 
                height_pl = 9, width_pl = 14, pheno_name = paste0('reverse_EggerInt_specPheno',pheno_name), show_rownames = T, 
                title_sub = 'MR-Egger reverse intercept: signed -log10(pvalue)')



