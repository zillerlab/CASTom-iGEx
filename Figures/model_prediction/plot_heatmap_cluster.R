# plot specific cluster (CAD)

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Liver'
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissues_name)
#####################################################################################################################
source(functR)
pheat_pl_gr <- function(mat, type_mat, height_pl = 10, width_pl = 7, color_df, outFile, cap = NA){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  
  tmp_mat <- as.matrix(mat[, !colnames(mat) %in% c('id', 'tissue')])
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  rownames(tmp_mat) <- mat$id
  
  mat_row <- data.frame(tissue = mat$tissue)
  rownames(mat_row) <- rownames(tmp_mat)
  
  mat_colors_gr <- list(tissue = color_df$color)
  names(mat_colors_gr$tissue) <- unique(mat$tissue)
  
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = T, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=T, cluster_cols=F, border_color = NA,
                    annotation_colors = mat_colors_gr,
                    annotation_row = mat_row, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 10, fontsize = 10, 
                    main =  sprintf("%s", type_mat),
                    cellwidth = 15, 
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  save_pheatmap_png(hm_pl, paste0(outFile, '.png'), height =height_pl , width =width_pl)
  save_pheatmap_pdf(hm_pl, paste0(outFile, '.pdf'), height =height_pl , width =width_pl)
  
}

pheat_pl <- function(mat, cl, type_mat, height_pl = 10, width_pl = 7, outFile, cap = NA, res_pl = 200){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # match mat and cl
  cl <- cl[match(colnames(mat), cl$id),]
  
  tmp_mat <- as.matrix(mat)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  tmp_mat <- tmp_mat[,id]
  
  mat_col <- data.frame(group = paste0('cl', cl$gr[id]))
  rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_gr <- list(group = rep(c('#444444', '#C1C1C1'),P)[1:P])
  names(mat_colors_gr$group) <- paste0('cl', 1:P)
  
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = F, show_rownames = T, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=T, cluster_cols=F, border_color = NA, annotation_colors = mat_colors_gr,
                    annotation_col = mat_col, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 10, fontsize = 12, 
                    main =  sprintf("%s", type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  save_pheatmap_png(hm_pl, paste0(outFile, '.png'), height =height_pl , width =width_pl, res = res_pl)
  save_pheatmap_pdf(hm_pl, paste0(outFile, '.pdf'), height =height_pl , width =width_pl)
  
}




color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))


#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)

df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
id_feat <- res_tscore$test_diff_gr$id[res_tscore$test_diff_gr$pval_corr< 0.05]
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% id_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:35]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(pathGO_original[, colnames(pathGO_original) %in% keep_feat])
rownames(red_mat) <- res_pathGO$res_pval$path[match(rownames(red_mat), res_pathGO$res_pval[, 1])]
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'path_GO', height_pl = 7, width_pl = 10,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Adipose_Visceral_Omentum'
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 11, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(pathGO_original[, colnames(pathGO_original) %in% keep_feat])
rownames(red_mat) <- res_pathGO$res_pval$path[match(rownames(red_mat), res_pathGO$res_pval[, 1])]
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'path_GO', height_pl = 7, width_pl = 10,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Adipose_Subcutaneous'
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)

df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Adrenal_Gland'
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)

df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(pathGO_original[, colnames(pathGO_original) %in% keep_feat])
rownames(red_mat) <- res_pathGO$res_pval$path[match(rownames(red_mat), res_pathGO$res_pval[, 1])]
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'path_GO', height_pl = 7, width_pl = 10,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

#######################################################################################################
library(corrplot)

# plot NMI clustering:
NMI_mat <- read.table('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric_nmiCL_singleTissue.txt', h=T, stringsAsFactors=F)
NMI_mat <- as.matrix(NMI_mat)
rownames(NMI_mat) <- colnames(NMI_mat)
color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
tissues_name <- colnames(NMI_mat)
# correlation
diag(NMI_mat) <- NA
col <- colorRampPalette(brewer.pal(9, 'PuBuGn'))(100)
ord <- corrMatOrder(NMI_mat, order="hclust", hclust.method = 'ward.D')
newcolours <- color_tissues$color[match(tissues_name, color_tissues$tissues)][ord]
title_pl <- sprintf('NMI tscore z-scaled')
outFold <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/'

pdf(file = sprintf('%s/NMI_%s_%s.pdf', outFold, 'tscore', 'zscaled'), width = 9, height = 6, compress = F, pointsize = 12)
corrplot(NMI_mat, type="upper", order = 'hclust', hclust.method = 'ward.D',
         tl.col = newcolours, tl.cex=1.2,
         col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
         addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.9, mar = c(0,0,1,0))
dev.off()
png(file = sprintf('%s/NMI_%s_%s.png', outFold, 'tscore', 'zscaled'), width = 9, height = 6, res = 300, units = 'in')
corrplot(NMI_mat, type="upper", order = 'hclust', hclust.method = 'ward.D',
         tl.col = newcolours, tl.cex=1.2,
         col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
         addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.9, mar = c(0,0,1,0))
dev.off()

########################################
################ Asthma ################
########################################


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Whole_Blood'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/', tissues_name)
#####################################################################################################################

source(functR)

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))


#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Lung'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/Asthma_pheno/Asthma_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 1)

red_mat <- t(pathGO_original[, colnames(pathGO_original) %in% keep_feat])
rownames(red_mat) <- res_pathGO$res_pval$path[match(rownames(red_mat), res_pathGO$res_pval[, 1])]
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'path_GO', height_pl = 7, width_pl = 10,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)


#######################################
################ T1D ##################
#######################################

#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Adipose_Subcutaneous'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissues_name)
#####################################################################################################################

source(functR)

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 9, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 13, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Liver'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissues_name)
#####################################################################################################################


color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 11, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'))


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 13, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)


#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Pancreas'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 11, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 13, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)

#####################################################################################################################
color_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
tissues_name <- 'Whole_Blood'
setwd('/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/')
clustFile_tscore <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathR <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_Reactome_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
clustFile_pathGO <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/path_GO_zscaled_clusterCases_PGmethod_HKmetric.RData', tissues_name) 
functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
outFold <- sprintf('OUTPUT_GTEx/predict_UKBB/%s/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/T1D_clustering/', tissues_name)
#####################################################################################################################

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

res_tscore <- get(load(clustFile_tscore))
res_pathR <- get(load(clustFile_pathR))
res_pathGO <- get(load(clustFile_pathGO))

#### heatmap: scaled version (sd=1, mean=0)
tscore_original <- sapply(1:ncol(res_tscore$input_data), function(x) res_tscore$input_data[, x]/res_tscore$res_pval[res_tscore$res_pval[,2] == colnames(res_tscore$input_data)[x], 7])
colnames(tscore_original) <- colnames(res_tscore$input_data)

pathR_original <- sapply(1:ncol(res_pathR$input_data), function(x) res_pathR$input_data[, x]/res_pathR$res_pval[res_pathR$res_pval[,1] == colnames(res_pathR$input_data)[x], 12])
colnames(pathR_original) <- colnames(res_pathR$input_data)

pathGO_original <- sapply(1:ncol(res_pathGO$input_data), function(x) res_pathGO$input_data[, x]/res_pathGO$res_pval[res_pathGO$res_pval[,1] == colnames(res_pathGO$input_data)[x], 14])
colnames(pathGO_original) <- colnames(res_pathGO$input_data)


df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(tscore_original))

df_gr_mean[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(tscore_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(tscore_original)

gr_tscore_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
identical(res_tscore$test_diff_gr$id, gr_tscore_original$mean$id)
mat <- gr_tscore_original$mean[res_tscore$test_diff_gr$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_tscoreOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, width_pl = 7, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)

red_mat <- t(tscore_original[, colnames(tscore_original) %in% keep_feat])
pheat_pl(mat = red_mat, cl = res_tscore$cl_best, type_mat = 'tscore', height_pl = 7, width_pl = 6,  outFile = paste0(file_name, '_heatmap_gr'), 
         cap = 3.5, res_pl = 300)

### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathR_original))

df_gr_mean[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathR_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathR_original)

gr_pathR_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathR_original), pval = apply(pathR_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathR_original$mean[test_diff$pval_corr < 0.05,]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathROriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_Reactome', height_pl = 7, width_pl = 11, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)


### plot pathways
df_gr_mean <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))
df_gr_sd <- matrix(ncol = length(unique(res_tscore$cl_best$gr)), nrow = ncol(pathGO_original))

df_gr_mean[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) mean(x[res_tscore$cl_best$gr == y]))))
df_gr_sd[,] <- t(apply(pathGO_original, 2, function(x) 
  sapply(sort(unique(res_tscore$cl_best$gr)), function(y) sd(x[res_tscore$cl_best$gr == y]))))
df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(res_tscore$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(pathGO_original)

gr_pathGO_original <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
# test differences
test_diff <- data.frame(id = colnames(pathGO_original), pval = apply(pathGO_original, 2, function(x) kruskal.test(x = x, g = factor(res_tscore$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
mat <- gr_pathGO_original$mean[test_diff$pval_corr < 0.05,]
mat <- mat[!duplicated(mat$id),]
keep_feat <- c()
for(i in 1:length(unique(res_tscore$cl_best$gr))){
  keep_feat <- c(keep_feat, mat$id[order(abs(mat[, i]), decreasing = T)[1:20]])
}
keep_feat <- unique(keep_feat)
mat <- mat[mat$id %in% keep_feat, ]
mat$tissue <- tissues_name

mat$id <- res_pathGO$res_pval$path[match(mat$id, res_pathGO$res_pval[, 1])]

file_name <- sprintf('%stscore_zscaled_clusterCases_PGmethod_HKmetric_pathGOOriginal', outFold)
pheat_pl_gr(mat, type_mat = 'path_GO', height_pl = 7, width_pl = 13, color_df = color_tissues, outFile = paste0(file_name, '_mean_heatmap_gr'), cap = 2)



