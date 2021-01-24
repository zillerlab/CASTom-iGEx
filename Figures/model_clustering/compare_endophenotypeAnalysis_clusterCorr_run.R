# charachterize cluster results: prove on CAD 

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


parser <- ArgumentParser(description="compare cluster correlation and endophenotype analysis")
parser$add_argument("--corr_cl_file", type = "character", help = "")
parser$add_argument("--endopheno_analysis_file", type = "character", help = "")
parser$add_argument("--pval_FDR_pheno", type = "double", default = 0.05, help = "")
parser$add_argument("--pval_pheno_show", type = "double", default = 0.001, help = "")
parser$add_argument("--thr_plot", type = "double", default = 0.2, help = "")
parser$add_argument("--color_pheno_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--outFold", type = "character", help = "")

args <- parser$parse_args()
corr_cl_file <- args$corr_cl_file
endopheno_analysis_file <- args$endopheno_analysis_file
pval_FDR_pheno <- args$pval_FDR_pheno
pval_pheno_show <- args$pval_pheno_show
color_pheno_file <- args$color_pheno_file
pheno_name <- args$pheno_name
thr_plot <- args$thr_plot
outFold <- args$outFold

#####################################################################################################
# tissue <- 'Liver'
# setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/')
# corr_cl_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/PearC_clusterCases_CAD_HARD_groupCorrelation_relatedPheno.RData', tissue)
# endopheno_analysis_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', tissue)
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue)
# pval_FDR_pheno <- 0.05
# pval_pheno_show <- 0.001
# thr_plot <- 0.2
# color_pheno_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
#####################################################################################################

color_pheno <- read.table(color_pheno_file, h=T, stringsAsFactors = F)

endo_res <- read.delim(endopheno_analysis_file, h=T, stringsAsFactors = F, sep = '\t')
endo_res_tot <- endo_res[!is.na(endo_res$pvalue), ]
endo_res <- endo_res_tot[endo_res_tot$pval_corr <= pval_FDR_pheno | endo_res_tot$pvalue <= pval_pheno_show,  ]

corr_cl <- get(load(corr_cl_file))
corr_cl$pval_FDRcorr <- apply(corr_cl$pval[, -ncol(corr_cl$pval)], 2, function(x) p.adjust(p = x, method = 'BH'))

# subset with endophenotypes associations
pheno_common <- rev(intersect(corr_cl$pheno$pheno, endo_res$pheno_id))
corr_cl_red <- list(corr = corr_cl$corr[match(pheno_common, corr_cl$pheno$pheno),], pval = corr_cl$pval[match(pheno_common, corr_cl$pheno$pheno),], 
                    pval_FDRcorr = corr_cl$pval_FDRcorr[match(pheno_common, corr_cl$pheno$pheno),], pheno = corr_cl$pheno[match(pheno_common, corr_cl$pheno$pheno),])

gr_name <- colnames(corr_cl$corr[, -ncol(corr_cl$pval)])
corr_cl_red$pheno$name_plot <- corr_cl_red$pheno$Field
corr_cl_red$pheno$name_plot[!is.na(corr_cl_red$pheno$meaning)] <- paste0(corr_cl_red$pheno$meaning[!is.na(corr_cl_red$pheno$meaning)], '\n', corr_cl_red$pheno$Field[!is.na(corr_cl_red$pheno$meaning)])

####################### plot functions ###########################
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

plot_heatmap_gr <- function(mat_gr, mat_pval_FDR = NULL, pheno_ann, pheno_info, outFold = outFold, 
                            color_groups, height_pl = 17, width_pl=11, cap_val = 1, show_rownames = T){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  coul[50] <- '#ffffff'
  
  # mat_gr[mat_pval>0.05] <- 0
  labels_sign <- matrix('', nrow = nrow(mat_gr), ncol = ncol(mat_gr))
  if(!is.null(mat_pval_FDR)){
    labels_sign[mat_pval_FDR <= 0.05] <- '*'  
  }
  
  tmp_mat <- as.matrix(mat_gr)
  # tmp_mat[is.na(tmp_mat)] <- 0
  
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
  
  mat_col <- data.frame(gr_type = color_groups$type[match(colnames(tmp_mat), color_groups$type)])
  rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_t <- list(gr_type = unique(color_groups$color))
  names(mat_colors_t$gr_type) <- unique(color_groups$type)
  
  new_color = list(pheno = mat_colors[[1]], gr_type = mat_colors_t[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = show_rownames, 
                    display_numbers = labels_sign, fontsize_number = 12, number_color = 'black',
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    na_col = "grey90",
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 16)
  
  
  save_pheatmap_png(hm_pl, sprintf("%sclusterCases_%s_groupCorrelation_relatedPheno.png", outFold, pheno_name), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%sclusterCases_%s_groupCorrelation_relatedPheno.pdf", outFold, pheno_name), height = height_pl, width =width_pl)
  
}

color_pheno$color[color_pheno$pheno_type == 'ICD10_Endocrine'] <- '#7d7d7dff'
color_pheno$color[color_pheno$pheno_type == 'Family_history'] <- '#cd8500ff'
color_groups <- data.frame(type = gr_name, color = pal_d3(palette = 'category20')(length(gr_name)))

##################################################################

plot_heatmap_gr(mat_gr = corr_cl_red$corr[, gr_name],  mat_pval_FDR =  corr_cl_red$pval_FDRcorr[, gr_name], pheno_info = corr_cl_red$pheno, 
                pheno_ann = color_pheno[match(unique(corr_cl_red$pheno$pheno_type), color_pheno$pheno_type),], cap_val = min(c(1, max(abs(corr_cl_red$corr[, gr_name]), na.rm = T))), 
                color_groups = color_groups, height_pl = 8, width_pl = 10, show_rownames = T, outFold = paste0(outFold, 'subset_'))

# plot phenotype significant in at least 1 group
id_sign <- which(rowSums(abs(corr_cl$corr[,gr_name]) >= thr_plot)>0)
corr_cl_sign <- list(corr = corr_cl$corr[id_sign,], pval = corr_cl$pval[id_sign,], 
                     pval_FDRcorr = corr_cl$pval_FDRcorr[id_sign,], pheno = corr_cl$pheno[id_sign,])

corr_cl_sign$pheno$name_plot <- corr_cl_sign$pheno$Field
corr_cl_sign$pheno$name_plot[!is.na(corr_cl_sign$pheno$meaning) & grepl('Diagnoses', corr_cl_sign$pheno$Field)] <- corr_cl_sign$pheno$meaning[!is.na(corr_cl_sign$pheno$meaning) & grepl('Diagnoses', corr_cl_sign$pheno$Field)]
corr_cl_sign$pheno$name_plot[!is.na(corr_cl_sign$pheno$meaning) & !grepl('Diagnoses', corr_cl_sign$pheno$Field)] <- paste(corr_cl_sign$pheno$Field[!is.na(corr_cl_sign$pheno$meaning) & !grepl('Diagnoses', corr_cl_sign$pheno$Field)], 
                                                                                                                          corr_cl_sign$pheno$meaning[!is.na(corr_cl_sign$pheno$meaning) & !grepl('Diagnoses', corr_cl_sign$pheno$Field)])

plot_heatmap_gr(mat_gr = corr_cl_sign$corr[, gr_name],  mat_pval_FDR =  corr_cl_sign$pval_FDRcorr[, gr_name], pheno_info = corr_cl_sign$pheno, 
                pheno_ann = color_pheno[match(unique(corr_cl_sign$pheno$pheno_type), color_pheno$pheno_type),], cap_val = min(c(1, max(abs(corr_cl_sign$corr[, gr_name]), na.rm = T))), 
                color_groups = color_groups, height_pl = 10, width_pl = 10, show_rownames = T, outFold = paste0(outFold, 'thr0.2_'))

# compute overall concordance rate: significant correlations only
thr_pval <- c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
df_repr <- list()

for(j in 1:length(thr_pval)){
  
  df_repr[[j]] <- data.frame(comp = gr_name, n_pheno = NA, n_same_sign = NA, fraction_rep = NA, pval = NA, thr_corr = thr_pval[j])
  
  for(i in 1:length(gr_name)){
    
    id <- corr_cl$pval[, gr_name[i]] <= thr_pval[j] & !is.na(corr_cl$pval[, gr_name[i]])
    pheno_c <- intersect(corr_cl$pheno$pheno[id], endo_res_tot$pheno_id)
    c_gr <- corr_cl$corr[match(pheno_c, corr_cl$pheno$pheno), gr_name[i]]
    beta_endo <- endo_res_tot[endo_res_tot$comp %in% paste0(gr_name[i], '_vs_all'), ]
    beta_endo <- beta_endo$z[match(pheno_c, beta_endo$pheno_id)]
    df_repr[[j]]$n_pheno[i] <- length(pheno_c)
    df_repr[[j]]$n_same_sign[i] <- sum(sign(c_gr*beta_endo) == 1)
    df_repr[[j]]$fraction_rep[i] <- sum(sign(c_gr*beta_endo) == 1)/length(pheno_c)
    df_repr[[j]]$pval[i] <- binom.test(x = df_repr[[j]]$n_same_sign[i], n = df_repr[[j]]$n_pheno[i])$p.value
    
  }
  
}

df_repr <- do.call(rbind, df_repr)

# plot distribution and save table
write.table(file = sprintf('%sclusterCases_%s_groupCorrelation_relatedPheno_concordance.txt', outFold, pheno_name), x = df_repr, quote = F, sep = '\t', row.names = F, col.names = T)

df_repr$new_thr <- -log10(df_repr$thr_corr)
df_repr$new_pval <- -log10(df_repr$pval)
df_repr$comp <- factor(df_repr$comp, levels = gr_name)

scale_pval <- (max(df_repr$new_pval) - min(df_repr$new_pval))/5

pl <- ggplot(data = df_repr, aes(x = new_thr, y = fraction_rep, color = comp, group = comp, size = new_pval))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('correlation p-value threshold')+ylab('Fraction of concordance')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  scale_x_continuous(breaks=-log10(thr_pval), labels = thr_pval)+
  scale_size_continuous(breaks = seq(min(df_repr$new_pval), max(df_repr$new_pval), length.out = 5), 
                        labels = formatC(10^(-seq(min(df_repr$new_pval), max(df_repr$new_pval), length.out = 5)), format = "e", digits = 2))+
  guides(color=guide_legend(title = ''), size = guide_legend(title = 'binomial test\np-value'))+
  scale_color_d3()

ggsave(filename = sprintf('%sclusterCases_%s_groupCorrelation_relatedPheno_concordance.png', outFold, pheno_name), plot = pl, width = 5, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sclusterCases_%s_groupCorrelation_relatedPheno_concordance.pdf', outFold, pheno_name), plot = pl, width = 5, height = 4.5, dpi = 500, compress = F)


