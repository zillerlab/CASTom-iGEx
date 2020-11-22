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

parser <- ArgumentParser(description="correlation plots: tissue specific CAD")
parser$add_argument("--corrRes_tissue_file", type = "character", nargs = '*', help = "results for correlation analysis tissue spec")
parser$add_argument("--mrRes_tissue_file", type = "character", nargs = '*', help = "results for mendelian randomization")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--pheno_list_file", type = "character", help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
corrRes_tissue_file <- args$corrRes_tissue_file
tissue_name <- args$tissue_name
pheno_name <- args$pheno_name
mrRes_tissue_file <- args$mrRes_tissue_file
pheno_list_file <- args$pheno_list_file
color_pheno_file <- args$color_pheno_file
color_tissues_file <- args$color_tissues_file
outFold <- args$outFold

#########################################################################################################################
# tissue_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# corrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_correlation_enrich_CAD_HARD_relatedPheno.RData')
# mrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_Mendelian_randomization_tot_path_pvalFDRrel0.05.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/perc0.3_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# pheno_name <- 'CAD_HARD'
# type_data <- 'tot_path'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# pheno_list_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_corr.txt'
# pheno_list_MR_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/keep_pheno_MR.txt'
# ######################################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissue_name,color_tissues$tissue),]
# color_tissues <- rbind(color_tissues, data.frame(tissues = 'All_tissues', color = 'black', type = 'All_GTEx', nsamples_train = NA))

pheno_keep <- read.table(pheno_list_file, h=F, stringsAsFactors = F, sep = '\t', check.names = F, quote = "")$V1
pheno_plot_MR <- read.table(pheno_list_MR_file, h=F, stringsAsFactors = F, sep = '\t', check.names = F, quote = "")$V1

# create matrix correlation/fisher OR
feat_corr <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_OR <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_mr_est <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_mr_est_se <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_mr_est_pval <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_mr_est_low <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))
feat_mr_est_up <- matrix(ncol = length(tissue_name), nrow = length(pheno_keep))

id_keep <- ifelse(type_data == 'tscore', 1, 2)

for(i in 1:length(tissue_name)){
  
  print(i)
  
  tmp_mr <- read.delim(mrRes_tissue_file[i], h=T, stringsAsFactors = F, sep = '\t')
  id <- match(pheno_keep, tmp_mr$names_field)
  feat_mr_est[,i] <- tmp_mr$MREgg_est[id]
  feat_mr_est_se[,i] <- tmp_mr$MREgg_est_se[id]
  feat_mr_est_pval[,i] <- tmp_mr$MREgg_est_pval[id]
  feat_mr_est_low[,i] <- tmp_mr$MREgg_est_CIl[id]
  feat_mr_est_up[,i] <- tmp_mr$MREgg_est_CIu[id]
  # feat_mr_est[tmp_mr$MREgg_est_pval[id] > 0.05, i] <- NA 
  # feat_mr_se_est[tmp_mr$MREgg_est_pval[id] > 0.05, i] <- NA 
  
  tmp <- get(load(corrRes_tissue_file[i]))
  # put NA value not significant after correction
  id <- match(pheno_keep, tmp$pheno$names_field)
  feat_corr[,i] <- tmp[[id_keep]]$cor_spearman[id]
  feat_corr[tmp[[id_keep]]$cor_pval_BHcorr[id] > 0.05, i] <- NA
  feat_OR[,i] <- tmp[[id_keep]]$fisher_OR[id]
  feat_OR[tmp[[id_keep]]$fisher_pval_BHcorr[id] > 0.05, i] <- NA 
  
}

pheno_info <- tmp$pheno[id, ]
pheno_info$name_plot <- pheno_info$names_field
pheno_info$name_plot[pheno_info$name_plot == 'Diagnoses - ICD10 : I10 Essential (primary) hypertension (41270_I10)'] <- 'Diagnoses - ICD10 : I10 Essential primary hypertension (41270_I10)'
pheno_info$name_plot[grepl('Diagnoses - ICD10 : ', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Diagnoses - ', pheno_info$name_plot)], function(x) strsplit(x, split = 'Diagnoses - ')[[1]][2])
pheno_info$name_plot[pheno_info$name_plot == 'Hand grip strength (left) (46)'] <- 'Hand grip strength left (46)'
pheno_info$name_plot[pheno_info$name_plot == 'Hand grip strength (right) (47)'] <- 'Hand grip strength right (47)'
pheno_info$name_plot[grepl('Heel bone mineral density', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Heel bone mineral density', pheno_info$name_plot)], function(x) paste0(strsplit(x, split = ' \\(BMD\\) ')[[1]], collapse = ''))
pheno_info$name_plot[grepl('Red blood cell', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Red blood cell', pheno_info$name_plot)], function(x) paste0(strsplit(x, split = ' \\(erythrocyte\\) ')[[1]], collapse = ''))
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
colnames(feat_corr) <- c(tissue_name)
colnames(feat_OR) <-  c(tissue_name)
colnames(feat_mr_est) <- colnames(feat_mr_est_se) <- colnames(feat_mr_est_pval) <- c(tissue_name)

feat_corr <- feat_corr[id, ]
feat_OR <- feat_OR[id, ]
feat_mr_est <- feat_mr_est[id, ]
feat_mr_est_se <- feat_mr_est_se[id, ]
feat_mr_est_pval <- feat_mr_est_pval[id, ]
feat_mr_est_low <- feat_mr_est_low[id, ]
feat_mr_est_up <- feat_mr_est_up[id, ]

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

plot_heatmap_mr <- function(type_mat, mat_mr, mat_pval, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 2, pheno_name, show_rownames = T){
  
  title_sub <- 'MR-Egger estimate'
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  coul[50] <- '#ffffff'
  
  mat_mr[mat_pval>0.05] <- NA
  
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
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s\nz-statistic %s', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 12)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_MRestimates_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}


# correlation
plot_heatmap_corr(type_mat = type_data, mat_cor = feat_corr[,tissue_name], pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                   pheno_info = pheno_info,  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T)

# Odds ratio
plot_heatmap_OR(type_mat = type_data, mat_OR = feat_OR[,tissue_name], pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], 
                  pheno_info = pheno_info,  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T)

# MR-Egger
plot_heatmap_mr(type_mat = type_data, mat_mr = feat_mr_est[,tissue_name], mat_pval = feat_mr_est_pval[,tissue_name],
                pheno_ann = color_pheno[match(unique(pheno_info$pheno_type),color_pheno$pheno_type),], cap_val = 1, 
                pheno_info = pheno_info,  color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = T)


#### plot specific MR results with pvalue and SE ###
id_MR <- which(pheno_info$names_field %in% pheno_plot_MR)
df_mr <- data.frame(MR_est = as.vector(feat_mr_est[id_MR, ]), MR_est_CIl = as.vector(feat_mr_est_low[id_MR, ]), 
                    MR_est_CIu = as.vector(feat_mr_est_up[id_MR, ]), MR_est_pval = as.vector(feat_mr_est_pval[id_MR, ]), 
                    pheno_name = rep(pheno_info$name_plot[id_MR], ncol(feat_mr_est)), tissue = unlist(lapply(colnames(feat_mr_est), function(x) rep(x, length(id_MR)))), stringsAsFactors = F)
df_mr$sign <- 'no'
df_mr$sign[df_mr$MR_est_pval <= 0.05] <- 'yes'
df_mr$sign <- factor(df_mr$sign, levels = c('no', 'yes'))
df_mr$pheno_name <- factor(df_mr$pheno_name, levels = pheno_info$name_plot[id_MR])
df_mr$tissue <- factor(df_mr$tissue, levels =color_tissues$tissues)

pl <-  ggplot(subset(df_mr, tissue %in% c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')), 
              aes(x = pheno_name, y = MR_est, shape = sign))+
  geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=MR_est_CIl, ymax=MR_est_CIu), width=.2, position=position_dodge(0.5))+
  theme_bw()+ 
  facet_wrap(.~tissue, nrow = 1)+
  ylab('MR-Egger estimate') + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 7), 
        plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 1), axis.text.y = element_text(size = 7.5), 
        strip.text = element_text(size=6.7))+
  scale_shape_manual(values=c(0, 15))+guides(shape=FALSE)+coord_flip()
  # scale_color_manual(values=color_tissues$color[color_tissues$tissue %in% c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')])

pl <- ggplot_gtable(ggplot_build(pl))
stripr <- which(grepl('strip-t', pl$layout$name))
fills <- color_tissues$color[color_tissues$tissue %in% c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Liver' ,'Whole_Blood')]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', pl$grobs[[i]]$grobs[[1]]$childrenOrder))
  pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  pl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$alpha <- 0.7
  k <- k+1
}

ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.png", outFold, type_data, pheno_name), width = 9, height = 3.5, plot = pl, device = 'png')
ggsave(filename = sprintf("%s%s_MRestimates_%s_relatedPheno_specPheno.pdf", outFold, type_data, pheno_name), width = 9, height = 3.5, plot = pl, device = 'pdf')




