library(ggplot2)
library(ggpubr)
library(qvalue)
library(corrplot)
library(RColorBrewer)
library(ggrepel)
library(argparse)
library(Matrix)
library(gdata)
library(VennDiagram)
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot pheno association")
parser$add_argument("--fold", type = "character", help = "fold output")
parser$add_argument("--functR", type = "character", help = "script with functions")
parser$add_argument("--fold_tissue", type = "character",default = 'NA', help = "fold input case example enrichment")
parser$add_argument("--fold_geno_input", type = "character", default = 'NA',help = "fold genotype input info (for case pathway)")
parser$add_argument("--keep_path_file", type = "character", default = 'NA', help = "pathways to keep")
parser$add_argument("--train_fold", type = "character", help = "train model fold")
parser$add_argument("--train_fold_original", type = "character", default = 'NA', help = "train model fold")
parser$add_argument("--color_file", type = "character", help = "file with tissues color code")
parser$add_argument("--gwas_known_file", type = "character", default = 'NA', help = "previusly associated genes")
parser$add_argument("--pheno", type = "character", help = "name phenotype")
parser$add_argument("--type_dat", type = "character", help = "name to append to plot file")
parser$add_argument("--pval_FDR", type = "double", default = 0.05, help = "pvalue threshold")


args <- parser$parse_args()
fold <- args$fold
fold_tissue <- args$fold_tissue
train_fold <- args$train_fold
color_file <- args$color_file
pheno <- args$pheno
type_dat <- args$type_dat
pval_FDR <- args$pval_FDR
fold_geno_input <- args$fold_geno_input
genes_known_file <- args$genes_known_file
keep_path_file <- args$keep_path_file
train_fold_original <- args$train_fold_original
functR <- args$functR

# fold <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/OUTPUT_all/'
# color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pval_FDR <- 0.05
# train_fold <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/'
# pheno <- 'SCZ'
# type_dat <- 'SCZ-PGC'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Figures/model_prediction/plot_prediction_functions.R'
# keep_path_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/keep_path_SCZ_plot.csv'
# fold_tissue <- c('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/',
#                  '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Brain_Cerebellar_Hemisphere/',
#                  '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Brain_Hypothalamus/')
# fold_geno_input <- c('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/Genotype_VariantsInfo_maf001_info06_CMC-PGCgwas-SCZ-PGCall_',
#                      '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/INPUT_DATA_GTEx/Genotype_VariantsInfo_maf001_info06_GTEx-PGCgwas-SCZ-PGCall_')
# train_fold_original <- c('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/DLPC_CMC/','/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/')
# #gwas_known_file = '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB_CADSchukert/GWAS_CAD_HARD_loci_annotatedGenes.txt'
# #priler_loci_file = '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB_CADSchukert/newloci_Priler_tscore_CAD_HARD.txt'

source(functR)

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)

# load results
tscore <- read.delim(sprintf('%stscore_pval_%s_covCorr.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathR <- read.delim(sprintf('%spath_Reactome_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathGO <- read.delim(sprintf('%spath_GO_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathWiki <- read.delim(sprintf('%scustomPath_WikiPath2019Human_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathCMC <- read.delim(sprintf('%scustomPath_CMC_GeneSets_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')

tissues <- unique(tscore$tissue)

# gene location
tscore$start_position <- NA
tscore$chrom <- NA
tscore$TSS_start <- NA

train_fold <- paste0(train_fold, tissues, '/')

for(i in 1:length(train_fold)){
  
  tmp <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[i]), h=T,stringsAsFactors = F)
  tmp <- tmp[match(tscore$ensembl_gene_id[tscore$tissue == tissues[i]], tmp$ensembl_gene_id),]
  tscore$start_position[tscore$tissue == tissues[i]] <- tmp$start_position
  tscore$chrom[tscore$tissue == tissues[i]] <- tmp$chrom
  tscore$TSS_start[tscore$tissue == tissues[i]] <- tmp$TSS_start 
  
}

###################################################################################################################################################

if(length(tissues)>2){
  
  ### correlation tissues ###
  tscore_cor <- create_cor(tissues_name = tissues, res = tscore, id_z = 7)
  pathR_cor <- create_cor(tissues_name = tissues, res = pathR, id_z = 12)
  pathGO_cor <- create_cor(tissues_name = tissues, res = pathGO, id_z = 14)
  pathwiki_cor <- create_cor(tissues_name = tissues, res = pathWiki, id_z = 12)
  
  pl_corr(tscore_cor, type_mat = 'tscore', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold, width_pl = 12, height_pl = 9)
  pl_corr(pathR_cor, type_mat = 'path_Reactome', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold,  width_pl = 12, height_pl = 9)
  pl_corr(pathGO_cor, type_mat = 'path_GO', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold, width_pl = 12, height_pl = 9)
  pl_corr(pathwiki_cor, type_mat = 'path_Wiki2019Human', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold, width_pl = 12, height_pl = 9)
  
}

### number of associated elements ###
tscore_nsgin <- creat_dfnsign(tissues_name = tissues, res = tscore, id_pval_corr = 10, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 10)
pathR_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathR, id_pval_corr = 15, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 15)
pathGO_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathGO, id_pval_corr = 17, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 17)
pathwiki_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathWiki, id_pval_corr = 15, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 15)

write.table(file = sprintf('%snsignificant_tscore.txt', fold), x = tscore_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(file = sprintf('%snsignificant_pathR.txt', fold), x = pathR_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(file = sprintf('%snsignificant_pathGO.txt', fold), x = pathGO_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(file = sprintf('%snsignificant_pathWiki.txt', fold), x = pathwiki_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')

pl_number_function(df = tscore_nsgin$plot, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathR_nsgin$plot, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathGO_nsgin$plot, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathwiki_nsgin$plot, type_mat = 'path_Wiki2019Human', outFold = fold, type_dat = type_dat)

tscore_red <- tscore
HLA_reg <- c(26000000, 34000000)
tscore_red <- tscore_red[!(tscore_red$chrom %in% 'chr6' & tscore_red$start_position <=HLA_reg[2] & tscore_red$start_position >= HLA_reg[1]) , ]
# # recompute pvalue correction:
# tmp <- list()
# for(i in 1:length(tissues)){
#   tmp[[i]] <- tscore_red[tscore_red$tissue == tissues[i], ]
#   tmp[[i]][, 10] <- p.adjust(tmp[[i]][,8], method = 'BH')
# }
# tscore_red <- do.call(rbind, tmp)
tscore_nsgin <- creat_dfnsign(tissues_name = tissues, res = tscore_red, id_pval_corr = 10, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 10)
write.table(file = sprintf('%snsignificant_tscore_noMHC.txt', fold), x = tscore_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')

### number of associated gene per tissue number ###
tscore_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = tscore, id_pval_corr = 10, pval_FDR = pval_FDR)
pathR_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = pathR, id_pval_corr = 15, pval_FDR = pval_FDR)
pathGO_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = pathGO, id_pval_corr = 17, pval_FDR = pval_FDR)
pathwiki_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = pathWiki, id_pval_corr = 15, pval_FDR = pval_FDR)

pl_numberSpec_function(df = tscore_nsgin_tissue, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_numberSpec_function(df = pathR_nsgin_tissue, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_numberSpec_function(df = pathGO_nsgin_tissue, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)
pl_numberSpec_function(df = pathwiki_nsgin_tissue, type_mat = 'path_Wiki2019Human', outFold = fold, type_dat = type_dat)

### manhattan plot ###
if(grepl('SCZ', pheno)){n_sign=10}
# if(grepl('CAD', pheno)){n_sign=4}
# if(grepl('MDD', pheno)){n_sign=20}
# if(grepl('T1D', pheno)){n_sign=10}

tscore_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = n_sign, gene = T)
pl_manhattan_function(data_input = tscore_df, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_manhattan_forpubl_function(data_input = tscore_df, type_mat = 'tscore', outFold = fold, type_dat = type_dat)

# remove MHC6 and plot again
tscore_red_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore_red, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = 2, gene = T)
pl_manhattan_function(data_input = tscore_red_df, type_mat = 'tscore', outFold = fold, type_dat = paste(type_dat, '(no MHC region)'))
pl_manhattan_forpubl_function(data_input = tscore_red_df, type_mat = 'tscore', outFold = fold, type_dat = paste(type_dat, '(no MHC region)'))

# Venn diagram for significnat genes #
if(gwas_known_file != 'NA'){
  
  venn_plot(gwas_known_file = gwas_known_file, tscore = tscore, pval_FDR = pval_FDR, type_dat = type_dat, type_mat = 'tscore')  

  # manhattan plot for new associaiton
  new_loci <- read.table(priler_loci_file, h=T, stringsAsFactors = F, sep = '\t')
  new_loci_ann <- read.table(priler_loci_ann_file, h=T, stringsAsFactors = F, sep = '\t')
  new_loci_ann = new_loci_ann[!new_loci_ann$best_GWAS_sign,]
  tscore_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, gene = T)
  tscore_df$df$external_gene_name <- sapply(tscore_df$df$name, 
                                            function(x) strsplit(x, split = '\n')[[1]][1])
  # modify annotation to plot only interested genes, only 1 gene per locus
  new_loci$id <- paste0(new_loci$external_gene_name, '\n', new_loci$tissue)
  tscore_df$df$sign_name[!tscore_df$df$name %in% new_loci$id] <- 'no'
  for(i in 1:nrow(new_loci_ann)){
    tmp <- new_loci_ann$external_gene_name[i]
    tmp <- strsplit(tmp, split = '[,]')[[1]]
    id_keep <- which.max(abs(tscore_df$df$zstat[tscore_df$df$external_gene_name %in% tmp]))
    names_excl <- tscore_df$df$name[tscore_df$df$external_gene_name %in% tmp][-id_keep]
    tscore_df$df$sign_name[tscore_df$df$name %in% names_excl] <- 'no'
  }
  tscore_df$df$name[tscore_df$df$sign_name == 'no'] <- '' 
  pl_manhattan_function(data_input = tscore_df, type_mat = 'tscore', outFold = paste0(fold, 'newloci_'), type_dat = type_dat)
  
}

###################################################################################################################################################
# specific to phenotype
pathR$type <- 'Reactome'
pathGO$type <- 'GO'
pathWiki$type <- 'Wiki2019Human'
pathCMC$type <- 'CMC'
common_h <- intersect(colnames(pathCMC), intersect(colnames(pathWiki), intersect(colnames(pathR), colnames(pathGO))))
tot_path <- rbind(pathR[, match(common_h, colnames(pathR))], pathGO[, match(common_h, colnames(pathGO))], 
                  pathWiki[, match(common_h, colnames(pathWiki))], pathCMC[match(common_h, colnames(pathCMC))])

# best_path <- best_res_fun(tot_path, tissues, id_pval = 13, n_top = 50, min_ngenes = 3, min_cov = 0)
tmp_path <- read.csv(keep_path_file, h=T, stringsAsFactors = F, sep = ',')
tmp_path <- apply(tmp_path, 1, function(x) paste0(x, collapse = '_'))
tot_path_id <- apply(tot_path[, c('path', 'tissue', 'type')], 1, function(x) paste0(x, collapse = '_'))
best_path <- tot_path[match(tmp_path, tot_path_id), ]
# save table
file_name <- sprintf('%s/table_%s_best_%s.txt', fold,  'path', type_dat)
write.table(x = best_path, file = file_name, quote = F, sep = '\t', col.names = T, row.names = F)
  
best_path$logpval <-  -log10(best_path[, 13]) 
best_path$zstat <- best_path[, 12]
plot_best_path(best_res = best_path, color_tissues = color_tissues,
               title_plot = pheno, type_mat = 'path', outFold = fold, type_dat = type_dat, 
                 tissues = unique(best_path$tissue), height_plot = 7, width_plot = 11, id_pval = 13)

### example pathway enrichment
id_pval <- 1
tissue <- 'DLPC_CMC'
pathways <- c('De novos:SCZ LoF', 'Regulation of Complement cascade')

for(pathway in pathways){
  
  print(pathway)
  color_tmp <- color_tissues$color[color_tissues$tissues == tissue]
  
  # load
  if(pathway == 'De novos:SCZ LoF'){
    info_res <- get(load(sprintf('%spval_Dx_pheno_covCorr_customPath_CMC_GeneSets.RData', fold_tissue[1])))
    id <- which(info_res$pathScore[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore[[id_pval]][[id]]  
  }else{
    info_res <- get(load(sprintf('%spval_Dx_pheno_covCorr.RData', fold_tissue[1])))
    id <- which(info_res$pathScore_reactome[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore_reactome[[id_pval]][[id]]  
  }
  
  gene_res <- info_res$tscore[[id_pval]]
  gene_info <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)
  resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
  train_fold_tissue <- train_fold[grepl(tissue, train_fold)]
  
  plot_showcase(gene_res = gene_res, gene_info = gene_info, genes_path = genes_path, 
                tissue = tissue, pathway = pathway, color_tmp = color_tmp, id_pval_path = 13, 
                pheno = pheno, fold = fold, resBeta = resBeta, train_fold_tissue = train_fold_tissue, fold_geno_input_tmp = fold_geno_input[1], 
                train_fold_original_tmp = train_fold_original[1], name_gwas_pval = 'PVAL')
}

###
tissue <- 'Brain_Cerebellar_Hemisphere'
pathways <- c('MAPK Signaling Pathway WP382','potassium ion transport')

for(pathway in pathways){
  
  print(pathway)
  color_tmp <- color_tissues$color[color_tissues$tissues == tissue]
  
  # load
  if(pathway == 'MAPK Signaling Pathway WP382'){
    info_res <- get(load(sprintf('%spval_Dx_pheno_covCorr_customPath_WikiPath2019Human.RData', fold_tissue[2])))
    id <- which(info_res$pathScore[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore[[id_pval]][[id]]
    id_pval_path <- 13
  }else{
    info_res <- get(load(sprintf('%spval_Dx_pheno_covCorr.RData', fold_tissue[2])))
    id <- which(info_res$pathScore_GO[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore_GO[[id_pval]][[id]]
    id_pval_path <- 15
  }
  
  # load
  gene_res <- info_res$tscore[[id_pval]]
  gene_info <- read.table(sprintf('%sresPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)
  resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
  train_fold_tissue <- train_fold[grepl(tissue, train_fold)]
  
  plot_showcase(gene_res = gene_res, gene_info = gene_info, genes_path = genes_path, 
                tissue = tissue, pathway = pathway, color_tmp = color_tmp, id_pval_path = id_pval_path, 
                pheno = pheno, fold = fold, resBeta = resBeta, train_fold_tissue = train_fold_tissue, fold_geno_input_tmp = fold_geno_input[2], 
                train_fold_original_tmp = train_fold_original[2], name_gwas_pval = 'PGC_PVAL')
  
}

###
tissue <- 'Brain_Hypothalamus'
pathway <- c('Axon guidance')

print(pathway)
color_tmp <- color_tissues$color[color_tissues$tissues == tissue]

# load
info_res <- get(load(sprintf('%spval_Dx_pheno_covCorr.RData', fold_tissue[3])))
id <- which(info_res$pathScore_reactome[[id_pval]]$path == pathway)
genes_path <- info_res$info_pathScore_reactome[[id_pval]][[id]]
gene_res <- info_res$tscore[[id_pval]]
gene_info <- read.table(sprintf('%sresPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)
resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
train_fold_tissue <- train_fold[grepl(tissue, train_fold)]

plot_showcase(gene_res = gene_res, gene_info = gene_info, genes_path = genes_path, 
              tissue = tissue, pathway = pathway, color_tmp = color_tmp, id_pval_path = 13, 
              pheno = pheno, fold = fold, resBeta = resBeta, train_fold_tissue = train_fold_tissue, fold_geno_input_tmp = fold_geno_input[2], 
              train_fold_original_tmp = train_fold_original[2], name_gwas_pval = 'PGC_PVAL')






