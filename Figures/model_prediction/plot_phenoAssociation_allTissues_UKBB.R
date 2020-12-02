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
functR <- args$functR

# fold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/'
# fold_tissue <- c('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Artery_Aorta/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/',
# '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/')
# color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pval_FDR <- 0.05
# train_fold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/'
# pheno <- 'CAD_HARD'
# type_dat <- 'CAD_HARD-UKBB'
# fold_geno_input <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/GTEX_v6/Genotyping_data/Genotype_VariantsInfo_GTEx-PGCgwas-CADgwas-CADall-UKBB_'
# gwas_known_file = '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB_CADSchukert/GWAS_CAD_HARD_loci_annotatedGenes.txt'
# priler_loci_file = '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB_CADSchukert/newloci_Priler_tscore_CAD_HARD.txt'
# keep_path_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/keep_path_CAD_plot.csv'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Figures/model_prediction/plot_prediction_functions.R'

source(functR)

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)

# load results
tscore <- read.delim(sprintf('%stscore_pval_%s_covCorr.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathR <- read.delim(sprintf('%spath_Reactome_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathGO <- read.delim(sprintf('%spath_GO_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
if(file.exists(sprintf('%spath_GO_pval_%s_covCorr_filt.txt', fold, pheno))){
  
}

train_fold_original <- train_fold

tissues <- unique(tscore$tissue)
if(grepl('CAD', pheno)){
  train_fold <- paste0(train_fold, tissues, '/200kb/CAD_GWAS_bin5e-2/')
}else{
  if(grepl('T1D', pheno)){
    train_fold <- paste0(train_fold, tissues, '/200kb/noGWAS/')
  }else{
    train_fold <- paste0(train_fold, tissues, '/')
  }
}

# gene location
tscore$start_position <- NA
tscore$chrom <- NA
tscore$TSS_start <- NA

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
  
  pl_corr(tscore_cor, type_mat = 'tscore', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold, width_pl = 10, height_pl = 7)
  pl_corr(pathR_cor, type_mat = 'path_Reactome', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold,  width_pl = 10, height_pl = 7)
  pl_corr(pathGO_cor, type_mat = 'path_GO', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold, width_pl = 10, height_pl = 7)
  
}
### number of associated elements ###
tscore_nsgin <- creat_dfnsign(tissues_name = tissues, res = tscore, id_pval_corr = 10, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 10)
pathR_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathR, id_pval_corr = 15, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 15)
pathGO_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathGO, id_pval_corr = 17, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 17)

write.table(file = sprintf('%snsignificant_tscore.txt', fold), x = tscore_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(file = sprintf('%snsignificant_pathR.txt', fold), x = pathR_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(file = sprintf('%snsignificant_pathGO.txt', fold), x = pathGO_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')

pl_number_function(df = tscore_nsgin$plot, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathR_nsgin$plot, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathGO_nsgin$plot, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)

tscore_red <- tscore
HLA_reg <- c(26000000, 34000000)
tscore_red <- tscore_red[!(tscore_red$chrom %in% 'chr6' & tscore_red$start_position <=HLA_reg[2] & tscore_red$start_position >= HLA_reg[1]) , ]
tscore_nsgin <- creat_dfnsign(tissues_name = tissues, res = tscore_red, id_pval_corr = 10, pval_FDR = pval_FDR, df_color = color_tissues, id_pval_corr_tot = 12)
write.table(file = sprintf('%snsignificant_tscore_noMHC.txt', fold), x = tscore_nsgin$table, quote = F, col.names = T, row.names = F, sep = '\t')

### number of associated gene per tissue number ###
tscore_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = tscore, id_pval_corr = 10, pval_FDR = pval_FDR)
pathR_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = pathR, id_pval_corr = 15, pval_FDR = pval_FDR)
pathGO_nsgin_tissue <- creat_dfnsign_tissueSpec(tissues_name = tissues, res = pathGO, id_pval_corr = 17, pval_FDR = pval_FDR)

pl_numberSpec_function(df = tscore_nsgin_tissue, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_numberSpec_function(df = pathR_nsgin_tissue, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_numberSpec_function(df = pathGO_nsgin_tissue, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)

### manhattan plot ###
if(grepl('SCZ', pheno)){n_sign=10}
if(grepl('CAD', pheno)){n_sign=4}
if(grepl('MDD', pheno)){n_sign=20}
if(grepl('T1D', pheno)){n_sign=10}

tscore_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = n_sign, gene = T)
# pathR_df <- create_df_manhattan_plot(tissues_name = tissues, res = pathR, id_pval = 13, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 1, n_sign = 2)
# pathGO_df <- create_df_manhattan_plot(tissues_name = tissues, res = pathGO, id_pval = 15, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = 2)

pl_manhattan_function(data_input = tscore_df, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
# pl_manhattan_function(data_input = pathR_df, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
# pl_manhattan_function(data_input = pathGO_df, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)

# if T1D remove MHC6 and plot again
if(grepl('T1D', pheno) | grepl('SCZ', pheno)){
  
  tscore_red <- tscore
  HLA_reg <- c(26000000, 34000000)
  tscore_red <- tscore_red[!(tscore_red$chrom %in% 'chr6' & tscore_red$start_position <=HLA_reg[2] & tscore_red$start_position >= HLA_reg[1]) , ]
  tscore_red_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore_red, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = 10, gene = T)
  pl_manhattan_function(data_input = tscore_red_df, type_mat = 'tscore', outFold = fold, type_dat = paste(type_dat, '(no MHC region)'))
  
}

### pathway plot ngenes info ###
# pathR_df <- create_df_manhattan_plot_path(tissues_name = tissues,  res = pathR, id_pval = 13, thr_genes = 0.1, pval_thr = 10^-4, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 1)
# pathGO_df <- create_df_manhattan_plot_path(tissues_name = tissues, res = pathGO, id_pval = 15, pval_thr = 10^-4, pval_FDR = pval_FDR,df_color = color_tissues, id_name = 2)
# 
# pl_manhattan_function_path(data_input = pathR_df, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat, pval_thr = 10^-4)
# pl_manhattan_function_path(data_input = pathGO_df, type_mat = 'path_GO', outFold = fold, type_dat = type_dat, pval_thr = 10^-4)

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
common_h <- intersect(colnames(pathR), colnames(pathGO))
tot_path <- rbind(pathR[, match(common_h, colnames(pathR))], pathGO[, match(common_h, colnames(pathGO))])

if(grepl('CAD', pheno)){
  
  # best_path <- best_res_fun(tot_path, tissues, id_pval = 13, n_top = 3, min_ngenes = 3, min_cov = 0.1)
  tmp_path <- read.csv(keep_path_file, h=T, stringsAsFactors = F)
  tmp_path <- apply(tmp_path, 1, function(x) paste0(x, collapse = '_'))
  tot_path_id <- apply(tot_path[, c('path', 'tissue', 'type')], 1, function(x) paste0(x, collapse = '_'))
  best_path <- tot_path[match(tmp_path, tot_path_id), ]
  # save table
  file_name <- sprintf('%s/table_%s_best_%s.txt', fold,  'path', type_dat)
  write.table(x = best_path, file = file_name, quote = F, sep = '\t', col.names = T, row.names = F)
  
  best_path$logpval <-  -log10(best_path[, 13]) 
  best_path$zstat <- best_path[, 12]
  plot_best_path(best_res = best_path, color_tissues = color_tissues, title_plot = pheno, type_mat = 'path', outFold = fold, type_dat = type_dat, 
                 tissues = tissues, height_plot = 7, width_plot = 10, id_pval = 13)
 
}

if(grepl('MDD', pheno)){
  
  best_pathR <- best_res_fun(pathR, tissues, id_pval = 13, n_top = 10)
  best_pathGO <- best_res_fun(pathGO, tissues, id_pval = 15,  n_top = 10)
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = type_dat,
                 tissues = c('DLPC_CMC', 'Whole_Blood'), height_plot = 5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno),  type_mat = 'path_GO', outFold = fold, type_dat = type_dat,
                 tissues = c('DLPC_CMC', 'Whole_Blood'), height_plot = 5)
  
}

if(grepl('T1D', pheno)){
  
  best_pathR <- best_res_fun(pathR, tissues, id_pval = 13, n_top = 15)
  best_pathGO <- best_res_fun(pathGO, tissues, id_pval = 15,  n_top = 15)
  
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = type_dat,
                 tissues = c('Adipose_Subcutaneous'), height_plot = 3.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno),  type_mat = 'path_GO', outFold = fold, type_dat = type_dat,
                 tissues = c('Adipose_Subcutaneous'), height_plot = 3.5, id_name = 2)
  
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = paste0(type_dat, 'v2'),
                 tissues = c('Pancreas'), height_plot = 3.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno), type_mat = 'path_GO',outFold = fold, type_dat = paste0(type_dat, 'v2'),
                 tissues = c('Pancreas'), height_plot = 3.5,  id_name = 2)
  
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = paste0(type_dat, 'v3'),
                 tissues = c('Whole_Blood'), height_plot = 3.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno), type_mat = 'path_GO',outFold = fold, type_dat = paste0(type_dat, 'v3'),
                 tissues = c('Whole_Blood'), height_plot = 3.5, id_name = 2)
  
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = paste0(type_dat, 'v4'),
                 tissues = c('Liver'), height_plot = 3.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno), type_mat = 'path_GO',outFold = fold, type_dat = paste0(type_dat, 'v4'),
                 tissues = c('Liver'), height_plot = 3.5, id_name = 2)
  
}


### example pathway enrichment
if(grepl('CAD', pheno)){
  
  if(pheno == 'CAD_HARD'){id_pval <- 1}
  if(pheno == 'CAD_SOFT'){id_pval <- 2}
  
  tissue <- 'Artery_Aorta'
  pathways <- c('Death Receptor Signalling', 'Cardiac conduction', 'Collagen formation')
  
  for(pathway in pathways){
    print(pathway)
    new_path <- paste0(strsplit(pathway, split = '[ ]')[[1]], collapse = '_')
    color_tmp <- color_tissues$color[color_tissues$tissues == tissue]
    
    # load
    info_res <- get(load(sprintf('%spval_CAD_pheno_covCorr.RData', fold_tissue[1])))
    id <- which(info_res$pathScore_reactome[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore_reactome[[id_pval]][[id]]
    gene_res <- info_res$tscore[[id_pval]]
    gene_info <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)
    
    id <- sapply(gene_res$external_gene_name, function(x) which(gene_info$external_gene_name == x))
    if(any(sapply(id, length)>1)){
      rm_id <- names(which(sapply(id, length)>1))
      gene_res <- gene_res[! gene_res$external_gene_name %in% rm_id, ]
      id <- id[-which(sapply(id, length)>1)]
      id <- unlist(id)
    }
    gene_info <- gene_info[id, ] 
    identical(gene_info$external_gene_name, gene_res$external_gene_name)
    
    id <- which(colnames(gene_info) == 'train_dev')-1
    gene_res <- cbind(gene_res, gene_info[,1:id])
    
    gene_tmp <- lapply(paste0('chr', 1:22), function(x) gene_res[gene_res$chrom==x,])
    start_pos <- sapply(gene_tmp, function(x) min(x$start_position))
    end_pos <- sapply(gene_tmp, function(x) max(x$end_position))
    df_add <- data.frame(start = start_pos, end = end_pos)
    df_add$start_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))+df_add$start
    df_add$add_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))
    new_pos <- mapply(function(x, y) x + y$start_position, x = df_add$add_plot, y = gene_tmp, SIMPLIFY = T)
    new_pos <- unlist(new_pos)
    
    gene_tmp <- do.call(rbind,gene_tmp)
    gene_tmp$new_pos <- new_pos
    
    new_df <- data.frame(pos = new_pos, val = -log10(gene_tmp[,8]), name = gene_tmp$external_gene_name, path = 0, stringsAsFactors = F, 
                         chr = sapply(gene_tmp$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))
    new_df$chr <- as.numeric(new_df$chr)
    new_df$path[new_df$name %in% genes_path$tscore$external_gene_name] <- 1 
    new_df$path[!(new_df$name %in% genes_path$tscore$external_gene_name) & (new_df$chr %% 2 ==0)] <- 2
    new_df$path <- factor(new_df$path, levels = c(1,0,2))
    new_df$type <- 0
    new_df$type[new_df$name%in% genes_path$tscore$external_gene_name] <- 1 
    new_df$type <- factor(new_df$type, levels = c(1,0))
    
    int_val = -log10(genes_path$path[,13])
    pl_manh <- ggplot(new_df, aes(x = pos, y = val, color = path, alpha = type)) + 
      geom_abline(intercept = int_val, slope = 0, color=color_tmp, linetype="dashed")+
      geom_point() + ggtitle('Tscore association with CAD')+
      scale_color_manual(values = c(color_tmp, "grey90", 'grey70')) +
      ylim(0, max(c(int_val, min(c(max(new_df$val), 8)))))+
      geom_text(x=new_df$pos[nrow(new_df)-1000], y=int_val+0.5, label=genes_path$path$path, size = 5, color = color_tmp)+
      geom_text_repel(
        data = subset(new_df, path == 1), aes(label = new_df$name[new_df$path==1]), size = 3.9, color = 'black', box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+
      xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15),
            axis.text.x=element_text(size = 11, angle = 45),
            axis.text.y=element_text(size = 11), legend.position = 'none')+
      scale_size_discrete(range = c(1.2, 0.3, 0.3))+
      scale_alpha_discrete(range = c(1, 0.5))+
      scale_x_continuous(breaks = df_add$start_plot, labels = c(1:22))
    ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, dpi=500, device = 'png')
    ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, device = 'pdf')
    
    ### plot SNPs gwas info for genes in the pathway ###
    df_genePath <- gene_tmp[gene_tmp$external_gene_name %in%  new_df$name[new_df$path == 1],]
    resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
    # load only certain chr
    beta_res <- NULL
    id_chr <- unique(df_genePath$chrom)
    for(i in id_chr){
      
      print(i)
      ind_chr <- as.numeric(strsplit(i, split = 'chr')[[1]][2])
      gene_pos <- read.table(sprintf('%shg19_ENSEMBL_TSS_%s_matched.txt', train_fold[grepl(tissue, train_fold)], i), header = T, stringsAsFactors = F)
      snp_pos <- read.table(sprintf('%s%s.txt', fold_geno_input, i), header = T, stringsAsFactors = F)
      id <- which(gene_pos$ensembl_gene_id %in% df_genePath$ensembl_gene_id)
      tmp <- resBeta[[ind_chr]][,id]
      if(length(id)==1){
        beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F, val = tmp[tmp!=0], id = which(tmp!=0), 
                                                     gene = gene_pos$external_gene_name[id]), snp_pos[which(tmp!=0), ]))  
      }else{
        for(j in id){
          tmp <- resBeta[[ind_chr]][,j]
          beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F,val = tmp[tmp!=0], id = which(tmp!=0), 
                                                       gene = gene_pos$external_gene_name[j]), snp_pos[which(tmp!=0), ]))    
        }
      }
      
    }
    
    # find double of beta_res and merge
    dup_snp <- names(which(table(beta_res$ID)>1))
    if(length(dup_snp)>1){
      for(i in 1:length(dup_snp)){
        id <- which(beta_res$ID == dup_snp[i])
        beta_res <- rbind(beta_res, beta_res[id[1], ])
        beta_res$val[nrow(beta_res)] <- mean(beta_res$val[id])
        beta_res$gene[nrow(beta_res)] <- paste(beta_res$gene[id], sep = '', collapse = '_')
        beta_res <- beta_res[-id, ]
      }
    }
    
    # plot only gwas results: SNPs that influece have not so much relevance
    df_start_end <- data.frame(chr = 1:22, start = 0, end = 0)
    for(i in 1:22){
      print(i)
      snp_pos <- read.table(sprintf('%shg19_SNPs_chr%i_matched.txt', train_fold_original, i), header = T, stringsAsFactors = F)
      df_start_end[i,2:3] <- c(snp_pos$position[1], snp_pos$position[nrow(snp_pos)]) 
    }
    
    beta_res <- lapply(1:22, function(x) beta_res[beta_res$CHR == x, ])
    df_start_end$add_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))
    new_pos <- mapply(function(x, y) x + y$POS, x = df_start_end$add_plot, y = beta_res, SIMPLIFY = T)
    new_pos <- unlist(new_pos)
    beta_res <- do.call(rbind, beta_res)
    beta_res$new_pos <- new_pos
    beta_res$transf_pvalue <- -log10(beta_res$CAD_p_dgc)
    df_start_end$start_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))+df_start_end$start
    
    pl_manh_snps <- ggplot(beta_res, aes(x = new_pos, y = transf_pvalue, color = val)) + 
      geom_point(size = 1) + ggtitle('GWAS pvalues for CAD')+
      xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 11),
            axis.text.x=element_text(size = 11, angle = 45),
            axis.text.y=element_text(size = 11), legend.position = 'bottom')+
      scale_x_continuous(breaks = df_start_end$start_plot, labels = c(1:22), limits = c(df_start_end$start_plot[1], sum(df_start_end$end)))+
      scale_color_gradient2(midpoint=0, low="blue", mid="grey",
                            high="red", space ="Lab")+
      labs(color = "reg. coefficient")
    
    ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, dpi=500, device = 'png')
    ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, device = 'pdf')
  
  }
  
  
  tissue <- 'Liver'
  pathway <- c('glutathione peroxidase activity')
    
  print(pathway)
    new_path <- paste0(strsplit(pathway, split = '[ ]')[[1]], collapse = '_')
    color_tmp <- color_tissues$color[color_tissues$tissues == tissue]
    
    # load
    info_res <- get(load(sprintf('%spval_CAD_pheno_covCorr.RData', fold_tissue[2])))
    id <- which(info_res$pathScore_GO[[id_pval]]$path == pathway)
    genes_path <- info_res$info_pathScore_GO[[id_pval]][[id]]
    gene_res <- info_res$tscore[[id_pval]]
    gene_info <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)
    
    id <- sapply(gene_res$external_gene_name, function(x) which(gene_info$external_gene_name == x))
    if(any(sapply(id, length)>1)){
      rm_id <- names(which(sapply(id, length)>1))
      gene_res <- gene_res[!gene_res$external_gene_name %in% rm_id, ]
      id <- id[-which(sapply(id, length)>1)]
      id <- unlist(id)
    }
    gene_info <- gene_info[id, ] 
    identical(gene_info$external_gene_name, gene_res$external_gene_name)
    
    id <- which(colnames(gene_info) == 'train_dev')-1
    gene_res <- cbind(gene_res, gene_info[,1:id])
    
    gene_tmp <- lapply(paste0('chr', 1:22), function(x) gene_res[gene_res$chrom==x,])
    start_pos <- sapply(gene_tmp, function(x) min(x$start_position))
    end_pos <- sapply(gene_tmp, function(x) max(x$end_position))
    df_add <- data.frame(start = start_pos, end = end_pos)
    df_add$start_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))+df_add$start
    df_add$add_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))
    new_pos <- mapply(function(x, y) x + y$start_position, x = df_add$add_plot, y = gene_tmp, SIMPLIFY = T)
    new_pos <- unlist(new_pos)
    
    gene_tmp <- do.call(rbind,gene_tmp)
    gene_tmp$new_pos <- new_pos
    
    new_df <- data.frame(pos = new_pos, val = -log10(gene_tmp[,8]), name = gene_tmp$external_gene_name, path = 0, stringsAsFactors = F, 
                         chr = sapply(gene_tmp$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))
    new_df$chr <- as.numeric(new_df$chr)
    new_df$path[new_df$name %in% genes_path$tscore$external_gene_name] <- 1 
    new_df$path[!(new_df$name %in% genes_path$tscore$external_gene_name) & (new_df$chr %% 2 ==0)] <- 2
    new_df$path <- factor(new_df$path, levels = c(1,0,2))
    new_df$type <- 0
    new_df$type[new_df$name%in% genes_path$tscore$external_gene_name] <- 1 
    new_df$type <- factor(new_df$type, levels = c(1,0))
    
    int_val = -log10(genes_path$path[,13])
    pl_manh <- ggplot(new_df, aes(x = pos, y = val, color = path, alpha = type)) + 
      geom_abline(intercept = int_val, slope = 0, color=color_tmp, linetype="dashed")+
      geom_point() + ggtitle('Tscore association with CAD')+
      scale_color_manual(values = c(color_tmp, "grey90", 'grey70')) +
      ylim(0, max(c(int_val, min(c(max(new_df$val), 8)))))+
      geom_text(x=new_df$pos[nrow(new_df)-1000], y=int_val+0.5, label=genes_path$path$path, size = 5, color = color_tmp)+
      geom_text_repel(
        data = subset(new_df, path == 1), aes(label = new_df$name[new_df$path==1]), size = 3.9, color = 'black', box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+
      xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15),
            axis.text.x=element_text(size = 11, angle = 45),
            axis.text.y=element_text(size = 11), legend.position = 'none')+
      scale_size_discrete(range = c(1.2, 0.3, 0.3))+
      scale_alpha_discrete(range = c(1, 0.5))+
      scale_x_continuous(breaks = df_add$start_plot, labels = c(1:22))
    ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, dpi=500, device = 'png')
    ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, device = 'pdf')
    
    ### plot SNPs gwas info for genes in the pathway ###
    df_genePath <- gene_tmp[gene_tmp$external_gene_name %in%  new_df$name[new_df$path == 1],]
    resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
    # load only certain chr
    beta_res <- NULL
    id_chr <- unique(df_genePath$chrom)
    for(i in id_chr){
      
      print(i)
      ind_chr <- as.numeric(strsplit(i, split = 'chr')[[1]][2])
      gene_pos <- read.table(sprintf('%shg19_ENSEMBL_TSS_%s_matched.txt', train_fold[grepl(tissue, train_fold)], i), header = T, stringsAsFactors = F)
      snp_pos <- read.table(sprintf('%s%s.txt', fold_geno_input, i), header = T, stringsAsFactors = F)
      id <- which(gene_pos$ensembl_gene_id %in% df_genePath$ensembl_gene_id)
      tmp <- resBeta[[ind_chr]][,id]
      if(length(id)==1){
        beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F, val = tmp[tmp!=0], id = which(tmp!=0), 
                                                     gene = gene_pos$external_gene_name[id]), snp_pos[which(tmp!=0), ]))  
      }else{
        for(j in id){
          tmp <- resBeta[[ind_chr]][,j]
          beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F,val = tmp[tmp!=0], id = which(tmp!=0), 
                                                       gene = gene_pos$external_gene_name[j]), snp_pos[which(tmp!=0), ]))    
        }
      }
      
    }
    
    # find double of beta_res and merge
    dup_snp <- names(which(table(beta_res$ID)>1))
    if(length(dup_snp)>1){
      for(i in 1:length(dup_snp)){
        id <- which(beta_res$ID == dup_snp[i])
        beta_res <- rbind(beta_res, beta_res[id[1], ])
        beta_res$val[nrow(beta_res)] <- mean(beta_res$val[id])
        beta_res$gene[nrow(beta_res)] <- paste(beta_res$gene[id], sep = '', collapse = '_')
        beta_res <- beta_res[-id, ]
      }
    }
    
    # plot only gwas results: SNPs that influece have not so much relevance
    df_start_end <- data.frame(chr = 1:22, start = 0, end = 0)
    for(i in 1:22){
      print(i)
      snp_pos <- read.table(sprintf('%shg19_SNPs_chr%i_matched.txt', train_fold_original, i), header = T, stringsAsFactors = F)
      df_start_end[i,2:3] <- c(snp_pos$position[1], snp_pos$position[nrow(snp_pos)]) 
    }
    
    beta_res <- lapply(1:22, function(x) beta_res[beta_res$CHR == x, ])
    df_start_end$add_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))
    new_pos <- mapply(function(x, y) x + y$POS, x = df_start_end$add_plot, y = beta_res, SIMPLIFY = T)
    new_pos <- unlist(new_pos)
    beta_res <- do.call(rbind, beta_res)
    beta_res$new_pos <- new_pos
    beta_res$transf_pvalue <- -log10(beta_res$CAD_p_dgc)
    df_start_end$start_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))+df_start_end$start
    
    pl_manh_snps <- ggplot(beta_res, aes(x = new_pos, y = transf_pvalue, color = val)) + 
      geom_point(size = 1) + ggtitle('GWAS pvalues for CAD')+
      xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 11),
            axis.text.x=element_text(size = 11, angle = 45),
            axis.text.y=element_text(size = 11), legend.position = 'bottom')+
      scale_x_continuous(breaks = df_start_end$start_plot, labels = c(1:22), limits = c(df_start_end$start_plot[1], sum(df_start_end$end)))+
      scale_color_gradient2(midpoint=0, low="blue", mid="grey",
                            high="red", space ="Lab")+
      labs(color = "reg. coefficient")
    
    ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, dpi=500, device = 'png')
    ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, device = 'pdf')
    
  
}



