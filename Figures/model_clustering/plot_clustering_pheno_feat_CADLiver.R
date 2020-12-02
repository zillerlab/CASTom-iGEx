# plot specific pathway, filter lists

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggradar))
suppressPackageStartupMessages(library(scales))
options(bitmapType = 'cairo', device = 'png')

## specific for Liver ##
tissue_name <- 'Liver'
clustFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissue_name)
res_cl <- get(load(clustFile))
color_tissue <- read.table('/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt', h=T, stringsAsFactors=F)

path_info <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_path_ReactomeOriginal_tscoreClusterCases_infoGenes.txt',  tissue_name), h=T, stringsAsFactors = F, sep = '\t')
path_go_info <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_path_GOOriginal_tscoreClusterCases_infoGenes.txt',tissue_name), h=T, stringsAsFactors = F, sep = '\t')
path_feat <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_path_ReactomeOriginal_tscoreClusterCases_featAssociation.txt',  tissue_name), h=T, stringsAsFactors = F, sep = '\t')
path_go_feat <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_path_GOOriginal_tscoreClusterCases_featAssociation.txt', tissue_name), h=T, stringsAsFactors = F, sep = '\t')
gene_feat <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_tscoreOriginal_tscoreClusterCases_featAssociation.txt', tissue_name), h=T, stringsAsFactors = F, sep = '\t')
gene_info <- read.delim(sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_tscoreOriginal_tscoreClusterCases_infoGenes.txt', tissue_name), h=T, stringsAsFactors = F, sep = '\t')

tissues_comp <- unique(path_info$tissue)
tmp <- sapply(tissues_comp, function(x) strsplit(x, split = '_')[[1]])
tissues_red <- sapply(tmp, function(y) paste0(sapply(y, function(x) substr(x, start = 1, stop = 1)), collapse = ''))
df_tissue <- data.frame(name = tissues_comp, short = tissues_red)
color_tissue <- color_tissue[color_tissue$tissue %in% tissues_comp, ]

# keep_path <- c('Golgi Associated Vesicle Biogenesis', 'Vesicle-mediated transport', 'clathrin-coated vesicle',
#                'Biosynthesis of DHA-derived SPMs', 'Free fatty acid receptors', 'Elastic fibre formation', 'glycine binding',
#                'G alpha (q) signalling events', 'Ethanol oxidation', 'Metabolism of carbohydrates',
#                'xenobiotic catabolic process', 'Transcriptional regulation of white adipocyte differentiation', 'antioxidant activity', 'DAP12 interactions', 'Endosomal/Vacuolar pathway', 'Interferon Signaling',
#                'Cytokine Signaling in Immune system', 'PD-1 signaling', 'Endogenous sterols', 'Steroid hormones',
#                'Cytochrome P450 - arranged by substrate type', 'Metabolism of lipids', 'ABC-family proteins mediated transport','Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)',
#                'heme binding', 'Synthesis of very long-chain fatty acyl-CoAs', 'Lysosphingolipid and LPA receptors')

keep_path <- c('Vesicle-mediated transport', 'clathrin-coated vesicle',
               'Free fatty acid receptors', 'Elastic fibre formation', 'glycine binding',
               'Ethanol oxidation', 'Metabolism of carbohydrates',
               'xenobiotic catabolic process', 'Transcriptional regulation of white adipocyte differentiation', 'antioxidant activity', 
               'DAP12 interactions', 'Endosomal/Vacuolar pathway', 'Interferon Signaling',
               'Cytokine Signaling in Immune system', 'PD-1 signaling', 'Endogenous sterols', 'Steroid hormones',
               'Cytochrome P450 - arranged by substrate type', 'Metabolism of lipids', 'ABC-family proteins mediated transport',
               'Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)',
               'heme binding', 'Lysosphingolipid and LPA receptors', 'Death Receptor Signalling', 'Cardiac conduction')


path_info_tot <- rbind(path_info, path_go_info[, match(colnames(path_info), colnames(path_go_info))])
path_feat_tot <- rbind(path_feat, path_go_feat[, match(colnames(path_feat), colnames(path_go_feat))])
path_info_tot$impr <- F
tmp <- list()
comp <- sort(unique(path_feat_tot$comp))
for(i in 1:length(keep_path)){
  print(i)
  tmp[[i]] <- path_info_tot[path_info_tot$path %in% keep_path[i], ]
  gene_names <- lapply(tmp[[i]]$genes_id, function(x) strsplit(x, split = '[,]')[[1]])
  for(j in 1:length(gene_names)){
    gene_new <- gene_feat[gene_feat$feat %in% gene_names[[j]] & gene_feat$tissue == tmp[[i]]$tissue[j], ]
    tmp[[i]]$impr[j] <- any(sapply(comp, function(x) all(gene_new$pval[gene_new$comp == x] > path_feat_tot$pval[path_feat_tot$tissue == tmp[[i]]$tissue[j] & path_feat_tot$feat == keep_path[i] & path_feat_tot$comp == x])))
  }
}
path_info_tot <- do.call(rbind, tmp)

path_info_tot$database <- 'Reactome'
path_info_tot$database[path_info_tot$path %in% path_go_info$path] <- 'Gene Ontology'
path_info_tot$new_id <- paste0(path_info_tot$path, ' (', df_tissue$short[match(path_info_tot$tissue, df_tissue$name)], ')')

path_feat_tot$new_id <- paste0(path_feat_tot$feat, ' (', df_tissue$short[match(path_feat_tot$tissue, df_tissue$name)], ')')
path_feat_tot$database <- 'Reactome'
path_feat_tot$database[path_feat_tot$feat %in% path_go_info$path] <- 'Gene Ontology'
path_feat_tot <- path_feat_tot[path_feat_tot$new_id %in% path_info_tot$new_id, ]
keep_new <- unique(path_feat_tot$new_id[path_feat_tot$pval <= 0.001 | (path_feat_tot$feat %in% c('Death Receptor Signalling', 'Cardiac conduction') & path_feat_tot$pval <= 0.05)])
path_feat_tot <- path_feat_tot[path_feat_tot$new_id %in% keep_new, ]
path_info_tot <- path_info_tot[path_info_tot$new_id %in% keep_new, ]

# load pathway scores:
pathR_input <- list()
pathGO_input <- list()
for(i in 1:length(tissues_comp)){
  print(i)
  pathR_featRelFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_ReactomeOriginal_tscoreClusterCases_featAssociation.RData', tissues_comp[i])
  tmp <- get(load(pathR_featRelFile))
  pathR_input[[i]] <- tmp$scaleData
  colnames(pathR_input[[i]]) <- paste0(colnames(pathR_input[[i]]), ' (', df_tissue$short[match(tissues_comp[i], df_tissue$name)], ')')

  pathGO_featRelFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GOOriginal_tscoreClusterCases_featAssociation.RData', tissues_comp[i])
  tmp <- get(load(pathGO_featRelFile))
  pathGO_input[[i]] <- tmp$scaleData
  colnames(pathGO_input[[i]]) <- tmp$res_pval$path[match(colnames(pathGO_input[[i]]), tmp$res_pval$path_id)]
  colnames(pathGO_input[[i]]) <- paste0(colnames(pathGO_input[[i]]), ' (', df_tissue$short[match(tissues_comp[i], df_tissue$name)], ')')
}
pathR_input <- do.call(cbind, pathR_input)
pathGO_input <- do.call(cbind, pathGO_input)
path_input_tot <- cbind(pathR_input, pathGO_input)
# old_path_input_tot <- path_input_tot
# path_input_tot <- old_path_input_tot

# merge_path <- c('DAP12 interactions', 'Endosomal/Vacuolar pathway', 'Cytokine Signaling in Immune system', 'PD-1 signaling')
merge_path <- keep_path
path_info_tot$new_id_v2 <- path_info_tot$new_id
path_feat_tot$new_id_v2 <- path_feat_tot$new_id

for(i in 1:length(merge_path)){
  
  print(i)
  
  tmp <- path_info_tot[path_info_tot$path %in% merge_path[i],]
  tmp_feat <- path_feat_tot[path_feat_tot$feat %in% merge_path[i],]
  tmp_feat_mat <- sapply(unique(tmp_feat$tissue), function(x) tmp_feat$estimates[tmp_feat$tissue == x])
  tmp_feat_pval <- sapply(unique(tmp_feat$tissue), function(x) tmp_feat$pval[tmp_feat$tissue == x])
  tmp_feat_mat[tmp_feat_pval > 0.05] <- 0
  
  if(ncol(tmp_feat_mat) >1){
    
    keep_tissue <- list()
    collapse_tissue <- list()
    tmp_dist <- as.matrix(dist(t(sign(tmp_feat_mat)), method = 'manhattan'))
    tissue_list <- lapply(1:nrow(tmp_dist), function(x) names(which(tmp_dist[x,] == 0)))
    tissue_list <- sapply(tissue_list, function(x) paste0(x, collapse = ','))
    tissue_list <- unname(tissue_list[!duplicated(tissue_list)])
    
    for(j in 1:length(tissue_list)){
      tmp_tissue <- strsplit(tissue_list[j], split = '[,]')[[1]]
      
      if(length(tmp_tissue) > 1){
        onlyt <- tmp_feat[tmp_feat$tissue %in% tmp_tissue, ]
        keep_tissue[[j]] <- onlyt$tissue[which.min(onlyt$pval)]
        collapse_tissue[[j]] <- setdiff(tmp_tissue, keep_tissue[[j]])
      }else{
        keep_tissue[[j]] <- tmp_tissue
        collapse_tissue[[j]] <- NA
      }
      
    }
  
    keep_tissue_tot <- unlist(keep_tissue)
    path_info_tot <- path_info_tot[!path_info_tot$new_id %in% tmp$new_id[!tmp$tissue %in% keep_tissue_tot],]
    path_feat_tot <- path_feat_tot[!path_feat_tot$new_id %in% tmp$new_id[!tmp$tissue %in% keep_tissue_tot],]
  
    for(j in 1:length(collapse_tissue)){
      if(!any(is.na(collapse_tissue[[j]]))){
        path_info_tot$new_id_v2[path_info_tot$new_id %in% tmp$new_id[tmp$tissue %in% keep_tissue[[j]]]] <- paste0(tmp$new_id[tmp$tissue %in% keep_tissue[[j]]], 
                                                                                                             ' -same in ', paste0(tissues_red[names(tissues_red) %in% collapse_tissue[[j]]], collapse = ' '), '-')
        
        path_feat_tot$new_id_v2[path_feat_tot$new_id %in% tmp$new_id[tmp$tissue %in% keep_tissue[[j]]]] <- paste0(tmp$new_id[tmp$tissue %in% keep_tissue[[j]]], 
                                                                                                           ' -same in ', paste0(tissues_red[names(tissues_red) %in% collapse_tissue[[j]]], collapse = ' '), '-')
      }
    }
 
  }
}

path_input_tot <- path_input_tot[, match(path_info_tot$new_id, colnames(path_input_tot))]
colnames(path_input_tot) <- path_info_tot$new_id_v2
path_info_tot$new_id <- path_info_tot$new_id_v2
path_feat_tot$new_id <- path_feat_tot$new_id_v2

### plot genes ###
gene_info <- gene_info[!is.na(gene_info$relevant_cluster_pruned), ]
gene_info <- gene_info[gene_info$relevant_cluster_pruned, ]
# gene_info <- gene_info[gene_info$relevant_cluster, ]
gene_info$new_id <- paste0(gene_info$external_gene_name, ' (', df_tissue$short[match(gene_info$tissue, df_tissue$name)], ')')
gene_feat$new_id <- paste0(gene_feat$feat, ' (', df_tissue$short[match(gene_feat$tissue, df_tissue$name)], ')')
gene_feat <- gene_feat[gene_feat$new_id %in% gene_info$new_id, ]
# futher filter based on association
keep_genes <- unique(gene_feat$new_id[gene_feat$pval <= 10^-6])
gene_info <- gene_info[gene_info$new_id %in% keep_genes, ]
gene_feat <- gene_feat[gene_feat$new_id %in% gene_info$new_id, ]

# load tscore:
tscore_input <- list()
for(i in 1:length(tissues_comp)){
  print(i)
  tscore_featRelFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscoreOriginal_tscoreClusterCases_featAssociation.RData', tissues_comp[i])
  tmp <- get(load(tscore_featRelFile))
  tscore_input[[i]] <- tmp$scaleData
  colnames(tscore_input[[i]]) <- paste0(colnames(tscore_input[[i]]), ' (', df_tissue$short[match(tissues_comp[i], df_tissue$name)], ')')
}
tscore_input <- do.call(cbind, tscore_input)
tscore_input <- tscore_input[, match(gene_info$new_id, colnames(tscore_input))]

merge_genes <- unique(gene_info$external_gene_name[duplicated(gene_info$external_gene_name)])
gene_info$new_id_v2 <- gene_info$new_id
gene_feat$new_id_v2 <- gene_feat$new_id

for(i in 1:length(merge_genes)){
  print(i)
  tmp <- gene_info[gene_info$external_gene_name %in% merge_genes[i],]
  tmp_feat <- gene_feat[gene_feat$feat %in% merge_genes[i],]
  tmp_feat_mat <- sapply(unique(tmp_feat$tissue), function(x) tmp_feat$estimates[tmp_feat$tissue == x])
  
  if(all(apply(tmp_feat_mat, 1, function(x) length(unique(sign(x)))) == 1)){
    
    keep_tissue <- tmp_feat$tissue[which.min(tmp_feat$pval)]
    excl_tissue <- setdiff(tmp$tissue, keep_tissue)
    gene_info <- gene_info[!gene_info$new_id %in% tmp$new_id[!tmp$tissue %in% keep_tissue],]
    gene_info$new_id_v2[gene_info$new_id %in% tmp$new_id[tmp$tissue %in% keep_tissue]] <- paste0(tmp$new_id[tmp$tissue %in% keep_tissue], 
                                                                                                         ' -same in ', paste0(tissues_red[names(tissues_red) %in% excl_tissue], collapse = ' '), '-')
    gene_feat <- gene_feat[!gene_feat$new_id %in% tmp$new_id[!tmp$tissue %in% keep_tissue],]
    gene_feat$new_id_v2[gene_feat$new_id %in% tmp$new_id[tmp$tissue %in% keep_tissue]] <- paste0(tmp$new_id[tmp$tissue %in% keep_tissue], 
                                                                                                         ' -same in ', paste0(tissues_red[names(tissues_red) %in% excl_tissue], collapse = ' '), '-')
  }
}

tscore_input <- tscore_input[, match(gene_info$new_id, colnames(tscore_input))]
colnames(tscore_input) <- gene_info$new_id_v2
gene_info$new_id <- gene_info$new_id_v2
gene_feat$new_id <- gene_feat$new_id_v2


#### function ####
pheat_pl_tot <- function(pheno_name, mat_score, info_feat, test_feat, pval_thr_est = 0.05, 
                         cl, height_pl = 10, width_pl = 7, outFile, cap = NA, cap_est = 0, res_pl = 200, tissue = NULL, color_tissue){
  
  if(!is.null(tissue)){
    info_feat <- info_feat[info_feat$tissue %in% tissue, ]
    mat_score <- mat_score[, match(info_feat$new_id, colnames(mat_score))]
    test_feat <- test_feat[info_feat$new_id %in% test_feat$new_id, ]
  }
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # cap mat if necessary
  tmp_mat <- as.matrix(mat_score)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # sort cl, match mat according cl
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  cl <- cl[id, ]
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  tmp_mat <- t(tmp_mat)
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  zstat_col_fun = colorRamp2(c(min(info_feat$Zstat, na.rm = T), 0, max(info_feat$Zstat, na.rm = T)), 
                             c("blue","#F0F0F0", "red"))
  if(length(unique(info_feat$ngenes_tscore)) == 1){
    ngenes_col_fun = colorRamp2(c(0,max(info_feat$ngenes_tscore)), c("white", "#035F1D"))
    perc_col_fun = colorRamp2(c(0, max(info_feat$ngenes_tscore/info_feat$ngenes_path)), c("white", "#316879"))
  }else{
    ngenes_col_fun = colorRamp2(c(min(info_feat$ngenes_tscore),max(info_feat$ngenes_tscore)), c("white", "#035F1D"))
    perc_col_fun = colorRamp2(c(min(info_feat$ngenes_tscore/info_feat$ngenes_path), max(info_feat$ngenes_tscore/info_feat$ngenes_path)), c("white", "#316879"))
  }
  database_color <- c('#260096', '#49bcd7')
  names(database_color) <- c('Gene Ontology', 'Reactome')
  tissue_color <- color_tissue$color
  names(tissue_color) <- color_tissue$tissues
  row_ha <- rowAnnotation(n_genes = info_feat$ngenes_tscore, perc = info_feat$ngenes_tscore/info_feat$ngenes_path, zstat = info_feat$Zstat,
                          database = info_feat$database, tissue = info_feat$tissue, 
                          annotation_label = list(tissue = 'tissue', database = 'database', n_genes = 'n. genes', perc = '% tot genes', zstat = sprintf('z-statistic %s', pheno_name)), 
                          col = list(tissue = tissue_color, database = database_color, n_genes = ngenes_col_fun, perc = perc_col_fun, zstat = zstat_col_fun))
  
  if(!is.na(cap_est)){
    test_feat$estimates[test_feat$estimates < -cap_est] <- -cap_est
    test_feat$estimates[test_feat$estimates > cap_est] <- cap_est
  }
  estimate_col_fun = colorRamp2(c(min(test_feat$estimates), 0, max(test_feat$estimates)), c("#00677B", "#F0F0F0", "#BF443B"))
  lgd_est = Legend(title = "wilcoxon estimates", col = estimate_col_fun, 
                   at = round(c(seq(min(test_feat$estimates), 0, length.out = 4), 
                                seq(0, max(test_feat$estimates), length.out = 4)[-1]), digits = 2),
                   labels = as.character(round(c(seq(min(test_feat$estimates), 0, length.out = 4), 
                                                 seq(0, max(test_feat$estimates), length.out = 4)[-1]), digits = 2)))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = sprintf("FDR pvalue < %s", as.character(pval_thr_est)))
  
  keep_path <- rownames(tmp_mat)
  # add pvalue info for each group
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat[test_feat$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_path, tmp_gr$new_id)]
    is_sign <- tmp_gr$pval_corr[match(keep_path, tmp_gr$new_id)] < pval_thr_est
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  row_ha_gr <- rowAnnotation(gr = anno_simple(df_est, col = estimate_col_fun, pch = df_pch, border = T, pt_gp = gpar(fontface = df_font)), 
                             annotation_label = 'wilcoxon\nestimates', annotation_name_side = 'bottom', annotation_name_rot = 0, simple_anno_size_adjust = T)
  font_path <- rep('plain', nrow(tmp_mat))
  font_path[info_feat$impr] <- 'bold'
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled scores", col = coul, cluster_rows = F, cluster_columns = FALSE,  show_row_dend = F, show_column_names = F, 
                   top_annotation = column_ha, column_split = cl$gr, column_title = NULL, 
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10, col = 'black', fontface = font_path),
                   left_annotation = row_ha, right_annotation = row_ha_gr, 
                   border = TRUE, use_raster = T)
  tot_pl_p <- hm_pl
  side_par <- round(max(sapply(rownames(tmp_mat), nchar))) + 30
  ht_list <- tot_pl_p
  
  png(file=paste0(outFile, '.png'), res = res_pl, units = 'in', width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
  pdf(file=paste0(outFile, '.pdf'), width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
}

pheat_pl_tscore <- function(pheno_name, mat_tscore, info_feat_tscore, test_feat_tscore, pval_thr_est = 0.05, 
                         cl, height_pl = 10, width_pl = 7, outFile, cap = NA, res_pl = 200, cap_est = 0, tissue = NULL, color_tissue){
  
  if(!is.null(tissue)){
    info_feat_tscore <- info_feat_tscore[info_feat_tscore$tissue %in% tissue, ]
    mat_tscore <- mat_tscore[, match(info_feat_tscore$new_id, colnames(mat_tscore))]
    test_feat_tscore <- test_feat_tscore[info_feat_tscore$new_id %in% test_feat_tscore$new_id, ]
  }
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # cap mat if necessary
  tmp_mat <- as.matrix(mat_tscore)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # sort cl, match mat according cl
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  cl <- cl[id, ]
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  
  # order gene according location
  info_feat_tscore <- info_feat_tscore[order(info_feat_tscore$start_position), ]
  info_feat_tscore <- info_feat_tscore[order(as.numeric(sapply(info_feat_tscore$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))), ]
  keep_gene <- info_feat_tscore$new_id
  tmp_mat <- tmp_mat[, match(keep_gene, colnames(tmp_mat)), drop = F]
  tmp_mat <- t(tmp_mat)
  chr_fact <- factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom))
  test_feat_tscore <- test_feat_tscore[test_feat_tscore$new_id %in% keep_gene, ,  drop = F]
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  mat_colors_chr <- list(chrom = rep(c('#7C7C7C', '#C1C1C1'),length(unique(info_feat_tscore$chrom)))[1:length(unique(info_feat_tscore$chrom))])
  names(mat_colors_chr$chrom) <- unique(info_feat_tscore$chrom)
  
  zstat_col_fun = colorRamp2(c(min(c(info_feat_tscore$Zstat), na.rm = T), 0, max(c(info_feat_tscore$Zstat), na.rm = T)), 
                             c("blue","#F0F0F0", "red"))
  # add pvalue info for each group
  if(!is.na(cap_est)){
    test_feat_tscore$estimates[test_feat_tscore$estimates < -cap_est] <- -cap_est
    test_feat_tscore$estimates[test_feat_tscore$estimates > cap_est] <- cap_est
  }
  estimate_col_fun = colorRamp2(c(min(test_feat_tscore$estimates), 0, max(test_feat_tscore$estimates)), c("#00677B", "#F0F0F0", "#BF443B"))
  lgd_est = Legend(title = "wilcoxon estimates", col = estimate_col_fun, 
                   at = round(c(seq(min(test_feat_tscore$estimates), 0, length.out = 4), 
                                seq(0, max(test_feat_tscore$estimates), length.out = 4)[-1]), digits = 2),
                   labels = as.character(round(c(seq(min(test_feat_tscore$estimates), 0, length.out = 4), 
                                                 seq(0, max(test_feat_tscore$estimates), length.out = 4)[-1]), digits = 2)))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = sprintf("FDR pvalue < %s", as.character(pval_thr_est)))
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  
  tissue_color <- color_tissue$color
  names(tissue_color) <- color_tissue$tissues
  row_ha <- rowAnnotation(chrom = anno_block(gp = gpar(fill = mat_colors_chr$chrom),
                                             labels = factor(names(mat_colors_chr$chrom), levels = names(mat_colors_chr$chrom)),
                                             labels_gp = gpar(col = "black", fontsize = 10), 
                                             labels_rot = 0), 
                          zstat = info_feat_tscore$Zstat,
                          tissue = info_feat_tscore$tissue, 
                          col = list(zstat = zstat_col_fun, tissue = tissue_color), 
                          annotation_label = list(tissue = 'tissue', zstat = sprintf('z-statistic %s', pheno_name)), 
                          annotation_name_gp = gpar(col = 'white'))
  
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat_tscore[test_feat_tscore$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_gene, tmp_gr$new_id)]
    is_sign <- tmp_gr$pval_corr[match(keep_gene, tmp_gr$new_id)] < pval_thr_est
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  row_ha_gr <- rowAnnotation(gr = anno_simple(df_est, col = estimate_col_fun, pch = df_pch, border = T, pt_gp = gpar(fontface = df_font)), 
                             annotation_label = '', annotation_name_side = 'top', annotation_name_rot = 0, simple_anno_size_adjust = T)
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled\nT-scores", col = coul, cluster_rows = FALSE, cluster_columns = FALSE,  show_column_names = F, 
                   top_annotation = column_ha, column_split = cl$gr, column_title = NULL,
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10),
                   left_annotation = row_ha,  row_split  = factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom)), row_title = NULL, row_gap = unit(0, "mm"),
                   right_annotation = row_ha_gr, 
                   border = TRUE, use_raster = T, show_heatmap_legend = F)
  tot_pl <- hm_pl
  ht_list <- tot_pl
  
  png(file=paste0(outFile, '.png'), res = res_pl, units = 'in', width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T)
  dev.off()
  
  pdf(file=paste0(outFile, '.pdf'), width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T)
  dev.off()
 
}

##################

pheat_pl_tot(pheno_name = 'CAD', mat_score = path_input_tot, info_feat = path_info_tot, test_feat = path_feat_tot, cl = res_cl$cl_best, pval_thr_est = 0.05, 
             outFile = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_pathOriginal_tscoreClusterCases',  tissue_name), 
             res_pl = 300, height_pl = 12, width_pl = 19, cap = 3, cap_est = 0.6, color_tissue = color_tissue)


pheat_pl_tscore(pheno_name = 'CAD', mat_tscore = tscore_input, info_feat_tscore = gene_info, test_feat_tscore = gene_feat , cl = res_cl$cl_best,  pval_thr_est = 0.05, 
                outFile = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_tscoreOriginal_tscoreClusterCases',  tissue_name), 
                res_pl = 150, height_pl = 20, width_pl = 13, cap = 3, cap_est = 1 , color_tissue = color_tissue)


##################################################################
# plot specific treatment response
tissue_name <- 'Liver'
treatmentResponsePairwiseFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_TreatResponse_pairwise.txt', tissue_name)
res_treat <- fread(treatmentResponsePairwiseFile, h=T, stringsAsFactors = F, data.table = F)
# consider only medicine specific for CAD
p1 <- res_treat[res_treat$pheno_Field %in% c('LDL direct','Lymphocyte count', 'C-reactive protein') & 
                  res_treat$treat_meaning %in% c('Cholesterol lowering medication'),]

p2 <- res_treat[res_treat$pheno_Field %in% c('Platelet count', 'Platelet crit', 'Platelet distribution width') & 
                  res_treat$treat_meaning %in% c('Paracetamol', 'Aspirin', 'Vitamin D'),]

p3 <- res_treat[res_treat$pheno_Field %in% c('Haemoglobin concentration', 'Haematocrit percentage', 'Mean corpuscular haemoglobin concentration') & 
                  res_treat$treat_meaning %in% c('Iron'),]

p4 <- res_treat[res_treat$pheno_Field %in% c('C-reactive protein') & 
                  res_treat$treat_meaning %in% c('Glucosamine'),]

p5 <- res_treat[res_treat$pheno_Field %in% c('Glycated haemoglobin (HbA1c)') & 
                  res_treat$treat_meaning %in% c('Folic acid or Folate (Vit B9)'),]

res_treat <- rbind(p1, p2, p3, p4, p5)

res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
df_red <- res_treat

df_red$new_id <- df_red$pheno_Field
df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
df_red$type_res <- 'beta'
df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Paracetamol', 'Aspirin', 'Vitamin D', 'Iron','Cholesterol lowering medication', 'Glucosamine', 'Folic acid or Folate (Vit B9)'))
df_red$sign <- 'no'
df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
df_red$gr1 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
df_red$gr2 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])

tot_gr <- sort(unique(c(as.character(df_red$gr1),  as.character(df_red$gr2))))

df_tot_gr <- list()
for(i in 1:length(tot_gr)){
  
  tmp <- df_red[df_red$gr1 == tot_gr[i] | df_red$gr2 == tot_gr[i], ]
  feat_name <- unique(tmp$comb_name)
  df_tot_gr[[i]] <- data.frame(comb_name =c(),gr =c(), sign =c(),pheno_class =c(),new_id =c(),treat_meaning  =c(),
                               gr_ORorBeta =c(),gr_CI_low =c(), gr_CI_up =c())
  for(j in 1:length(feat_name)){
    df_new <- data.frame(comb_name = feat_name[j], gr = tot_gr[i], sign = any(tmp$sign[tmp$comb_name == feat_name[j]] == 'yes'), 
                         pheno_class = tmp$pheno_class[tmp$comb_name == feat_name[j]][1], 
                         new_id = tmp$new_id[tmp$comb_name == feat_name[j]][1], 
                         treat_meaning = tmp$treat_meaning[tmp$comb_name == feat_name[j]][1], 
                         pheno_type = tmp$pheno_type[tmp$comb_name == feat_name[j]][1])
    name_gr <- ifelse(any(tot_gr[i] %in% tmp$gr1[tmp$comb_name == feat_name[j]]), 'gr1', 'gr2')
    df_new$gr_ORorBeta <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_ORorBeta')][1]
    df_new$gr_CI_low <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_low')][1]
    df_new$gr_CI_up <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_up')][1]
    df_tot_gr[[i]] <- rbind(df_tot_gr[[i]], df_new)
  }
}
df_tot_gr <- do.call(rbind, df_tot_gr)
df_tot_gr$gr <- factor(df_tot_gr$gr, levels = unique(df_tot_gr$gr))

len_w <- 7
len_h <- 7 + length(unique(df_tot_gr$gr))*0.1

gr_color <- pal_d3(palette = 'category20')(length(unique(df_tot_gr$gr)))

pl_beta_p1 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Cholesterol lowering medication')), 
                      aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=gr_color)+guides(shape=FALSE)+
  coord_flip()

pl_beta_p2 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Paracetamol', 'Aspirin', 'Vitamin D')), 
                      aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
  theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=gr_color)+guides(shape=FALSE)+
  coord_flip()

pl_beta_p3 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Iron')), 
                      aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
  theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=gr_color)+guides(shape=FALSE)+
  coord_flip()

pl_beta_p4 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Glucosamine')), 
                      aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=gr_color)+guides(shape=FALSE)+
  coord_flip()

pl_beta_p5 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Folic acid or Folate (Vit B9)')), 
                      aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=gr_color)+guides(shape=FALSE)+
  coord_flip()

outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/',tissue_name)

tot_pl <- ggarrange(plotlist = list(pl_beta_p1, pl_beta_p2, pl_beta_p3), align = 'h', nrow = 1, common.legend = T, widths = c(1, 2.2, 1.4))
ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll_subset1.png', outFold, 'Cases', 'CAD'), width = 15, height = 3.5, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll_subset1.pdf', outFold, 'Cases', 'CAD'), width = 15, height = 3.5, plot = tot_pl, device = 'pdf')

tot_pl <- ggarrange(plotlist = list(pl_beta_p4, pl_beta_p5), align = 'v', nrow = 1, widths = c(1, 1.1))
ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll_subset2.png', outFold, 'Cases', 'CAD'), width = 7, height = 1.7, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll_subset2.pdf', outFold, 'Cases', 'CAD'), width = 7, height = 1.7, plot = tot_pl, device = 'pdf')

####################################################
### spider plot selected pathways and phenotypes ###
####################################################

## pathways ##
keep_path_sp <- c('Vesicle-mediated transport (L)', 'Free fatty acid receptors (AS)',
                  'Ethanol oxidation (HAA) -same in AA HLV-', 
                  'Interferon Signaling (HAA) -same in AVO AS AA HLV CT WB L AG AC-', 'Endogenous sterols (CS) -same in HAA L-',
                  'Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) (HAA) -same in HLV-', 
                  'heme binding (AA) -same in WB-')

path_input_red <- path_input_tot[, keep_path_sp]
gr_id <- sort(unique(res_cl$cl_best$gr))
df_mean <- t(sapply(gr_id, function(x) colMeans(path_input_red[res_cl$cl_best$gr == x,])))
df_mean <- cbind(data.frame(group = paste0('gr_',gr_id)), df_mean)
# df_mean <- data.frame(val = as.vector(df_mean), group = as.vector(sapply(gr_id, function(x) rep(paste0('gr_',x), ncol(path_input_red)))))
# df_mean$path <- rep(colnames(path_input_red), length(gr_id))
# df_mean$new_id <- df_mean$path
colnames(df_mean)[colnames(df_mean) == 'Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) (HAA) -same in HLV-'] <- 'Regulation of IGF transport and uptake by IGFBPs (HAA) -same in HLV-'
new_name <- unname(sapply(colnames(df_mean[, -1]), function(x) paste(strsplit(x, split = ' [(]')[[1]], collapse = '\n(')))
colnames(df_mean)[-1] <- new_name
df_mean[, -1] <- apply(df_mean[,-1], 2, rescale)

gr_color <- pal_d3(palette = 'category20')(length(unique(gr_id)))
pl <- ggradar(df_mean,  
              #grid.min = min(df_mean[, -1]), grid.max = max(df_mean[, -1]), grid.mid = 0, 
              #values.radar = round(c(min(df_mean[, -1]), 0, max(df_mean[, -1])), digits = 2),
              grid.min = 0, grid.max = 1, grid.mid = 0.5, 
              values.radar = c('0%', '50%', '100%'),
              group.colours = gr_color, 
        grid.label.size = 5,
        axis.label.size = 3.5, 
        group.point.size = 2,
        group.line.width = 1,
        legend.text.size= 10, 
        legend.position = 'top', 
        plot.extent.x.sf = 2, 
        plot.extent.y.sf = 1.2)
# pl + theme(plot.margin = margin(2, 4, 4, 4, "cm"))

ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotPathway_tscoreClusterCases.png',  tissue_name), device = 'png', width = 10, height = 10, dpi = 320)
ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotPathway_tscoreClusterCases.pdf',  tissue_name), device = 'pdf', width = 10, height = 10)


## phenotypes ## 
keep_pheno_sp <- c('LDL direct', 'Apolipoprotein A', 'Lymphocyte count', 'Haemoglobin concentration', 'Hyperlipidemia', 'Peripheral_vascular_disease', 'Age_stroke')
tmp <- get(load(sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', tissue_name)))
tmp_nom <- get(load(sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/nominalAnalysis_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', tissue_name)))

phenoDat_tot <- cbind(tmp$phenoDat[, match(tmp$phenoInfo$pheno_id[tmp$phenoInfo$Field %in% keep_pheno_sp],colnames(tmp$phenoDat))], tmp_nom$phenoDat[, colnames(tmp_nom$phenoDat) %in% keep_pheno_sp])
colnames(phenoDat_tot)[1:sum(tmp$phenoInfo$Field %in% keep_pheno_sp)] <- tmp$phenoInfo$Field[tmp$phenoInfo$Field %in% keep_pheno_sp]
# phenoDat_tot$Age_stroke <- rescale(phenoDat_tot$Age_stroke)
# phenoDat_tot[, 'Apolipoprotein A'] <- rescale(phenoDat_tot[, 'Apolipoprotein A'])
# phenoDat_tot[, 'LDL direct'] <- rescale(phenoDat_tot[, 'LDL direct'])
# phenoDat_tot[, 'Haemoglobin concentration'] <- rescale(phenoDat_tot[, 'Haemoglobin concentration'])
# phenoDat_tot[, 'Lymphocyte count'] <- rescale(phenoDat_tot[, 'Lymphocyte count'])

# compute mean across groups
df_mean <- t(sapply(gr_id, function(x) colMeans(phenoDat_tot[res_cl$cl_best$gr == x,], na.rm = T)))
df_mean <- cbind(data.frame(group = paste0('gr_',gr_id)), df_mean)
df_mean[, -1] <- apply(df_mean[,-1], 2, rescale)

pl <- ggradar(df_mean,  grid.min = 0, grid.max = 1, grid.mid = 0.5, 
              values.radar = c('0%', '50%', '100%'),
              group.colours = gr_color, 
              grid.label.size = 5,
              axis.label.size = 3.5, 
              group.point.size = 2,
              group.line.width = 1,
              legend.text.size= 10, 
              legend.position = 'top', 
              plot.extent.x.sf = 2, 
              plot.extent.y.sf = 1.2)
# pl + theme(plot.margin = margin(2, 4, 4, 4, "cm"))

ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotPheno_tscoreClusterCases.png',  tissue_name), device = 'png', width = 10, height = 10, dpi = 320)
ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotPheno_tscoreClusterCases.pdf',  tissue_name), device = 'pdf', width = 10, height = 10)

## treatment response ##
keep_treat_pheno_sp <- data.frame(pheno = c('LDL direct', 'Platelet crit', 'Platelet crit', 'Haemoglobin concentration', 'Glycated haemoglobin (HbA1c)', 'C-reactive protein'), 
                                  treat = c('Cholesterol lowering medication' , 'Paracetamol', 'Vitamin D', 'Iron', 'Folic acid or Folate (Vit B9)', 'Glucosamine'))
  
covDat <- fread('INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All_phenoAssoc_withMedication.txt', h=T, stringsAsFactors = F, data.table = F)
covDat <- covDat[match(res_cl$samples_id, covDat$Individual_ID), ]
treat_pheno <- colnames(covDat)[!colnames(covDat) %in% c('Individual_ID', 'Dx', paste0('PC',1:10), 'Age', 'Gender')]

phenoDat <- fread('INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeMatrix_CADHARD_All_phenoAssoc_withMedication.txt', h=T, stringsAsFactors = F, data.table = F)
phenoDat <- phenoDat[match(res_cl$samples_id, phenoDat$Individual_ID), ]
phenoInfo <- fread('INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeDescription_withMedication.txt', h=T, stringsAsFactors = F, data.table = F)
# consider only certain class of phenotypes
phenoInfo <- phenoInfo[phenoInfo$pheno_type %in% c('Arterial_stiffness', 'Blood_biochemistry', 'Blood_count', 'Blood_pressure', 'Body_size_measures', 'Hand_grip_strength', 'Impedance_measures', 
                                                   'Residential_air_pollution', 'Spirometry'),]
# exclude Nucleated red blood cell
phenoInfo <- phenoInfo[!grepl('Nucleated red blood cell',phenoInfo$Field), ]
phenoInfo <- phenoInfo[!grepl('Traffic intensity on the nearest road',phenoInfo$Field), ]
phenoInfo_treat <- fread('INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeDescription_covariateMatrix_withMedication.txt', h=T, stringsAsFactors = F, data.table = F)
phenoInfo_treat <- rbind(phenoInfo_treat,
                         data.frame(pheno_id = c('6153_6177_1', '6153_6177_2', '6153_6177_3'), FieldID = c('6153_6177', '6153_6177', '6153_6177'),
                                    Field = rep('Medication for cholesterol, blood pressure or diabetes', 3), Path = rep(NA, 3), Strata  = rep(NA, 3),
                                    Sexed = rep('Unisex', 3), Coding = rep(NA, 3), Coding_meaning = c('Cholesterol lowering medication', 'Blood pressure medication', 'Insulin'),
                                    original_type = rep('CAT_MULTIPLE', 3), transformed_type = rep('CAT_MUL_BINARY_VAR',3), nsamples = rep(NA,3), nsamples_T= rep(NA,3),
                                    nsamples_F= rep(NA,3), pheno_type = rep('Medication', 3)))

phenoInfo_treat <- phenoInfo_treat[!phenoInfo_treat$FieldID %in% c('6177', '6153'),]

phenoInfo_treat <- phenoInfo_treat[phenoInfo_treat$Coding_meaning %in% keep_treat_pheno_sp$treat, ]
phenoInfo <- phenoInfo[phenoInfo$Field %in% keep_treat_pheno_sp$pheno, ]

df_diff <- matrix(nrow = length(gr_id), ncol = nrow(keep_treat_pheno_sp))
for(i in 1:nrow(keep_treat_pheno_sp)){
  pv <- phenoDat[, colnames(phenoDat) == phenoInfo$pheno_id[phenoInfo$Field == keep_treat_pheno_sp$pheno[i]]]
  cv <- covDat[, colnames(covDat) == paste0('p', phenoInfo_treat$pheno_id[phenoInfo_treat$Coding_meaning == keep_treat_pheno_sp$treat[i]])]
  pv_gr <- lapply(gr_id, function(x) pv[res_cl$cl_best$gr == x])
  cv_gr <- lapply(gr_id, function(x) cv[res_cl$cl_best$gr == x])
  df_diff[, i] <- mapply(function(x, y) mean(x[y == 1], na.rm = T) - mean(x[y == 0], na.rm = T), x = pv_gr, y = cv_gr)
}
colnames(df_diff) <- paste0(keep_treat_pheno_sp$treat , '\n', keep_treat_pheno_sp$pheno)
df_diff <- cbind(data.frame(group = paste0('gr_',gr_id)), df_diff)
df_diff[, -1] <- apply(df_diff[,-1], 2, rescale)

pl <- ggradar(df_diff,  grid.min = 0, grid.max = 1, grid.mid = 0.5, 
              values.radar = c('0%', '50%', '100%'),
              group.colours = gr_color, 
              grid.label.size = 5,
              axis.label.size = 3.5, 
              group.point.size = 2,
              group.line.width = 1,
              legend.text.size= 10, 
              legend.position = 'top', 
              plot.extent.x.sf = 2, 
              plot.extent.y.sf = 1.2)
# pl + theme(plot.margin = margin(2, 4, 4, 4, "cm"))

ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotTreat_tscoreClusterCases.png',  tissue_name), device = 'png', width = 10, height = 10, dpi = 320)
ggsave(plot = pl, filename = sprintf('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/cl%s_spiderPlotTreat_tscoreClusterCases.pdf',  tissue_name), device = 'pdf', width = 10, height = 10)

