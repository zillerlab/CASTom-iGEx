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
# tissue_name <- c('Liver', 'Artery_Aorta')
# corrRes_tot_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/randomChoice_correlation_CAD_HARD_relatedPheno.RData'
# corrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/randomChoice_correlation_CAD_HARD_relatedPheno.RData')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/randomChoice_'
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt' 
# pheno_name <- 'CAD_HARD'
# diff_thr_plot <- 0.09
##########################################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissue_name,color_tissues$tissue),]

cor_tot <- get(load(corrRes_tot_file))

df_tot <- data.frame(tscore_pval = -log10(cor_tot$tscore$cor_pval), tscore_corr = cor_tot$tscore$cor_spearman, 
                     path_pval = -log10(cor_tot$pathScore$cor_pval), path_corr = cor_tot$pathScore$cor_spearman, 
                     pheno_name = cor_tot$pheno$names_field, pheno_type = cor_tot$pheno$pheno_type)

color_pheno <- color_pheno[color_pheno$pheno_type %in% df_tot$pheno_type, ]

df_tot$pheno_type <- factor(df_tot$pheno_type, levels = color_pheno$pheno_type)
df_tot$name_plot <- df_tot$pheno_name
df_tot$name_plot[df_tot$name_plot == 'Diagnoses - ICD10 : I10 Essential (primary) hypertension (41270_I10)'] <- 'Diagnoses - ICD10 : I10 Essential primary hypertension (41270_I10)'
df_tot$name_plot[abs(df_tot$tscore_corr) - abs(df_tot$path_corr) >= -diff_thr_plot & abs(df_tot$tscore_corr) - abs(df_tot$path_corr) <= diff_thr_plot] <- ''
df_tot$name_plot[grepl('Diagnoses - ICD10 : ', df_tot$name_plot)] <- sapply(df_tot$name_plot[grepl('Diagnoses - ', df_tot$name_plot)], function(x) strsplit(x, split = 'Diagnoses - ')[[1]][2])
df_tot$name_plot <- sapply(df_tot$name_plot, function(x) strsplit(x, split = ' \\(')[[1]][1])
df_tot$name_plot[is.na(df_tot$name_plot)] <- ''

# plot tscore vs pathway correlation:
file_name <- sprintf('%scorrelation_%s_relatedPheno_compare_tscore_path', outFold, pheno_name)
pl <- ggplot(data = df_tot, aes(x = tscore_corr, y = path_corr, color = pheno_type, label = name_plot))+
  geom_point(alpha = 0.8, size = 1)+
  geom_text_repel(segment.color = 'grey50', color = 'black', size = 2, min.segment.length = unit(0, 'lines'),
                  segment.alpha = 0.6,  force = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  xlab('gene T-score')+ylab('Pathway score')+ 
  theme_bw()+ ggtitle('Correlation Z-statistic')+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                                                       legend.text = element_text(size = 6),  legend.key.size = unit(0.3, "cm"), legend.title = element_blank())+
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = color_pheno$color)

ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 6, height = 4.5, dpi = 500)
ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 6, height = 4.5, dpi = 500, compress = F)

# heatmap for tissue specific
res_t <- list()
for(i in 1:length(tissue_name)){
  
  tmp <- get(load(corrRes_tissue_file[i]))
  
  df_tmp <- data.frame(tscore_corr = tmp$tscore$cor_spearman, 
                       path_corr = tmp$pathScore$cor_spearman, 
                       pheno_name = tmp$pheno$names_field, pheno_type = tmp$pheno$pheno_type)
  df_tmp$keep <- F
  df_tmp$keep[abs(abs(df_tmp$path_corr) - abs(df_tmp$tscore_corr)) >= diff_thr_plot & (abs(df_tmp$tscore_corr)>=0.3 | abs(df_tmp$path_corr)>=0.3)] <- T
  df_tmp$tissue <- tissue_name[i]
  res_t[[i]] <- df_tmp

}

# create dataframe for heatmap plot
keep_tot <- rowSums(sapply(res_t, function(x) x$keep))>0
tscore_mat <- sapply(1:length(tissue_name), function(x) res_t[[x]]$tscore_corr[keep_tot])
colnames(tscore_mat) <- tissue_name
path_mat <- sapply(1:length(tissue_name), function(x) res_t[[x]]$path_corr[keep_tot])
colnames(path_mat) <- tissue_name
pheno_ann_t <- res_t[[1]][keep_tot, c('pheno_name', 'pheno_type')]
pheno_ann_t$name_plot <- pheno_ann_t$pheno_name
pheno_ann_t$name_plot[pheno_ann_t$name_plot == 'Diagnoses - ICD10 : I10 Essential (primary) hypertension (41270_I10)'] <- 'Diagnoses - ICD10 : I10 Essential primary hypertension (41270_I10)'
pheno_ann_t$name_plot[grepl('Diagnoses - ICD10 : ', pheno_ann_t$name_plot)] <- sapply(pheno_ann_t$name_plot[grepl('Diagnoses - ICD10 : ', pheno_ann_t$name_plot)], function(x) strsplit(x, split = 'Diagnoses - ICD10 : ')[[1]][2])

# remove ICD10 Circulatory and family history
id <- pheno_ann_t$pheno_type %in%  c('ICD10_Circulatory_system', 'Family_history')
pheno_ann_t <- pheno_ann_t[!id,]
tscore_mat <- tscore_mat[!id,]
path_mat <- path_mat[!id, ]


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

plot_heatmap_split <- function(type_mat, mat_cor, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 0.6, pheno_name){
  
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
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = T, 
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s\nz-statistic %s', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 10)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s%s_correlation_%s_relatedPheno_heatmap.png", outFold, type_mat, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s%s_correlation_%s_relatedPheno_heatmap.pdf", outFold, type_mat, pheno_name), height = height_pl, width =width_pl)
  
}


# pvalue
plot_heatmap_split(type_mat = 'tscore', mat_cor = tscore_mat, pheno_ann = color_pheno, pheno_info = pheno_ann_t, 
                   color_tissues = color_tissues, height_pl = 12, width_pl = 13, pheno_name = pheno_name)

plot_heatmap_split(type_mat = 'pathScore', mat_cor = path_mat, pheno_ann = color_pheno, pheno_info = pheno_ann_t, 
                   color_tissues = color_tissues, height_pl = 12, width_pl = 13, pheno_name = pheno_name)


