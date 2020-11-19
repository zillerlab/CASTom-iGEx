# plot correlation with related phenotypes

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

parser <- ArgumentParser(description="correlation plots: general and tissue specific")
parser$add_argument("--corrRes_tot_file", type = "character", help = "results for correlation analysis")
parser$add_argument("--corrRes_tissue_file", type = "character", nargs = '*', help = "results for correlation analysis tissue spec")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--diff_thr_plot", type = "double", default = 0.09,  help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
corrRes_tot_file <- args$corrRes_tot_file
corrRes_tissue_file <- args$corrRes_tissue_file
tissue_name <- args$tissue_name
pheno_name <- args$pheno_name
color_pheno_file <- args$color_pheno_file
color_tissues_file <- args$color_tissues_file
diff_thr_plot <- args$diff_thr_plot
outFold <- args$outFold

#########################################################################################################################
# tissue_name <- c('Liver', 'Artery_Aorta', 'Whole_Blood', 'Artery_Coronary', 'Colon_Sigmoid')
# corrRes_tot_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/allPath_correlation_enrich_CAD_HARD_relatedPheno.RData'
# corrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/allPath_correlation_enrich_CAD_HARD_relatedPheno.RData')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/allPath_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# pheno_name <- 'CAD_HARD'
# diff_thr_plot <- 0.09
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
##########################################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissue_name,color_tissues$tissue),]
color_tissues <- rbind(color_tissues, data.frame(tissues = 'All_tissues', color = 'black', type = 'All_GTEx', nsamples_train = NA))

cor_tot <- get(load(corrRes_tot_file))
# create matrix correlation/fisher OR
tscore_corr <- matrix(ncol = length(tissue_name)+1, nrow = nrow(cor_tot$pheno))
tscore_OR <- matrix(ncol = length(tissue_name)+1, nrow = nrow(cor_tot$pheno))
path_corr <- matrix(ncol = length(tissue_name)+1, nrow = nrow(cor_tot$pheno))
path_OR <- matrix(ncol = length(tissue_name)+1, nrow = nrow(cor_tot$pheno))

# put NA value not significant after correction
tscore_corr[,1] <- cor_tot$tscore$cor_spearman
tscore_corr[cor_tot$tscore$cor_pval_BHcorr > 0.05, 1] <- NA 
tscore_OR[,1] <- cor_tot$tscore$fisher_OR
tscore_OR[cor_tot$tscore$fisher_pval_BHcorr > 0.05, 1] <- NA 

path_corr[,1] <- cor_tot$pathScore$cor_spearman
path_corr[cor_tot$pathScore$cor_pval_BHcorr > 0.05, 1] <- NA 
path_OR[,1] <- cor_tot$pathScore$fisher_OR
path_OR[cor_tot$pathScore$fisher_pval_BHcorr > 0.05, 1] <- NA 

# save pheno
pheno_info <- data.frame(pheno_name = cor_tot$pheno$names_field, pheno_type = cor_tot$pheno$pheno_type)
pheno_info$name_plot <- pheno_info$pheno_name
pheno_info$name_plot[pheno_info$name_plot == 'Diagnoses - ICD10 : I10 Essential (primary) hypertension (41270_I10)'] <- 'Diagnoses - ICD10 : I10 Essential primary hypertension (41270_I10)'
pheno_info$name_plot[grepl('Diagnoses - ICD10 : ', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Diagnoses - ', pheno_info$name_plot)], function(x) strsplit(x, split = 'Diagnoses - ')[[1]][2])
pheno_info$name_plot[grepl('Types of transport used', pheno_info$name_plot)] <- sapply(pheno_info$name_plot[grepl('Types of transport used', pheno_info$name_plot)], function(x) paste0(strsplit(x, split = ' \\(excluding work\\) ')[[1]], collapse = ''))
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
 

for(i in 1:length(tissue_name)){
  print(i)
  tmp <- get(load(corrRes_tissue_file[i]))
  # put NA value not significant after correction
  tscore_corr[,i+1] <- tmp$tscore$cor_spearman
  tscore_corr[tmp$tscore$cor_pval_BHcorr > 0.05, i+1] <- NA 
  tscore_OR[,i+1] <- tmp$tscore$fisher_OR
  tscore_OR[tmp$tscore$fisher_pval_BHcorr > 0.05, i+1] <- NA 
  
  path_corr[,i+1] <- tmp$pathScore$cor_spearman
  path_corr[tmp$pathScore$cor_pval_BHcorr > 0.05, i+1] <- NA 
  path_OR[,i+1] <- tmp$pathScore$fisher_OR
  path_OR[tmp$pathScore$fisher_pval_BHcorr > 0.05, i+1] <- NA 
  
}
colnames(path_OR) <- colnames(path_corr) <- c('All_tissues', tissue_name)
colnames(tscore_OR) <- colnames(tscore_corr) <- c('All_tissues', tissue_name)
tscore_corr <- cbind(pheno_info, tscore_corr)
tscore_OR <- cbind(pheno_info, tscore_OR)
path_OR <- cbind(pheno_info, path_OR)
path_corr <- cbind(pheno_info, path_corr)

# save results
write.table(file = sprintf("%s%s_correlation_%s_relatedPheno_combined.txt", outFold, 'tscore', pheno_name), x = tscore_corr, col.names = T, sep = '\t', row.names = F, quote = F)
write.table(file = sprintf("%s%s_enrichmentOR_%s_relatedPheno_combined.txt", outFold, 'tscore', pheno_name), x = tscore_OR, col.names = T, sep = '\t', row.names = F, quote = F)
write.table(file = sprintf("%s%s_correlation_%s_relatedPheno_combined.txt", outFold, 'path', pheno_name), x = path_corr, col.names = T, sep = '\t', row.names = F, quote = F)
write.table(file = sprintf("%s%s_enrichmentOR_%s_relatedPheno_combined.txt", outFold, 'path', pheno_name), x = path_OR, col.names = T, sep = '\t', row.names = F, quote = F)

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

plot_heatmap_split <- function(type_mat, mat_cor, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 0.6, pheno_name, show_rownames = T){
  
  title_sub <- 'Spearman Correlation'
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  
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


# pvalue
plot_heatmap_split(type_mat = 'tscore', mat_cor = tscore_corr[,c('All_tissues',tissue_name)], pheno_ann = color_pheno, pheno_info = tscore_corr[,!colnames(tscore_corr) %in% c('All_tissues',tissue_name)], 
                   color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = F)

plot_heatmap_split(type_mat = 'pathScore', mat_cor = path_corr[,c('All_tissues',tissue_name)], pheno_ann = color_pheno, pheno_info = path_corr[,!colnames(tscore_corr) %in% c('All_tissues',tissue_name)], 
                   color_tissues = color_tissues, height_pl = 13, width_pl = 14, pheno_name = pheno_name, show_rownames = F)


