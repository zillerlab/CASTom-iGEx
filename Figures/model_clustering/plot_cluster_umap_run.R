# plot umap structure for clustering, 
# also plot Age/Gender and PCs distribution

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(umap))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clusterin UMAP structure")
parser$add_argument("--clustFile", type = "character", help = "")
parser$add_argument("--tissue_name", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--tscore_featRelFile", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character", help = "")
parser$add_argument("--type_cluster_data", type = "character", help = "")
parser$add_argument("--type_input", type = "character", help = "")
parser$add_argument("--type_sim", type = "character", help = "")
parser$add_argument("--covDat_file", type = "character", help = "")
# parser$add_argument("--pathR_featRelFile", type = "character", default  = 'NA', help = "")
# parser$add_argument("--pathGO_featRelFile", type = "character",  default  = 'NA', help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
type_cluster <- args$type_cluster
type_cluster_data <- args$type_cluster_data
clustFile <- args$clustFile
tissue_name <- args$tissue_name
pheno_name <- args$pheno_name
tscore_featRelFile <- args$tscore_featRelFile
type_input <- args$type_input
type_sim <- args$type_sim
type_cluster_data <- args$type_cluster_data
covDat_file <- args$covDat_file
# pathR_featRelFile <- args$pathR_featRelFile
# pathGO_featRelFile <- args$pathGO_featRelFile
outFold <- args$outFold

#####################################################################################################################
# type_cluster <- 'Cases'
# type_cluster_data <- 'tscore'
# tissue_name <- 'Liver'
# pheno_name <- 'CAD'
# type_input <- 'zscaled'
# type_sim <- 'HK'
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All_phenoAssoc.txt'
# clustFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissue_name)
# tscore_featRelFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscoreOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# # pathR_featRelFile<- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_ReactomeOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# # pathGO_featRelFile<- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GOOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue_name)
# #####################################################################################################################

# get clustering structure
res_cl <- get(load(clustFile))
cl <- res_cl$cl_best
input_data <- res_cl$input_data
covDat <- read.table(covDat_file, h=T, stringsAsFactors = F, sep = '\t')
covDat <- covDat[match(res_cl$samples_id, covDat$Individual_ID), !colnames(covDat) %in% c('Dx', 'Individual_ID'), ]

P <- length(unique(cl$gr))
gr_color <- pal_d3(palette = 'category20')(P)

# umap
n_comp_umap <- 2
n_neigh_umap <- res_cl$best_k
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap

umap_res <- umap::umap(input_data, custom.settings)

df <- data.frame(component_1=umap_res$layout[,1], component_2=umap_res$layout[,2], gr = cl$gr)
df$gr <- factor(df$gr)

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr))+
  geom_point(size = 0.05, alpha = 0.8)+
  xlab('UMAP component 1')+ ylab('UMAP component 2')+
  scale_color_manual(values = gr_color)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4

ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_umap.png', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_umap.pdf', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')

# plot distribution Age, Gender and PCs
df_cov <- cbind(covDat, data.frame(gr = cl$gr))
df_cov$gr <- paste0('gr', df_cov$gr)
df_cov$gr <- factor(df_cov$gr, levels = paste0('gr', 1:P))
df_cov_PC <- data.frame(val = as.vector(as.matrix(df_cov[, paste0('PC', 1:10)])), PC = unlist(lapply(paste0('PC', 1:10), function(x) rep(x, nrow(df_cov)))), gr = rep(df_cov$gr, 10))
df_cov_PC$PC <- factor(df_cov_PC$PC, levels = paste0('PC', 1:10))

p <- ggboxplot(df_cov_PC, x = "gr", y = "val", fill = "gr", color = 'black', legend = 'none', outlier.size = 0.2, alpha = 0.8) + stat_compare_means(label = "p.format", size = 3) 
p <- ggpar(p, palette = gr_color, xlab = '', ylab = '', x.text.angle = 45)
p <- facet(p, facet.by = "PC", short.panel.labs = T, scales = 'free_y', nrow = 1)

ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_PCs.png', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 13, height = 4, plot = p, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_PCs.pdf', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 13, height = 4, plot = p, device = 'pdf')

df_cov_age <- data.frame(Age = df_cov$Age, gr = df_cov$gr)
res_cl$test_cov$pval[res_cl$test_cov$cov_id == 'Age'] = kruskal.test(df_cov_age$Age, g = df_cov_age$gr)$p.value

df_cov_sex <- data.frame(n = as.vector(table(df_cov$Gender, df_cov$gr)), Sex = rep(c('male', 'female'), P), gr = unlist(lapply(paste0('gr', 1:P), function(x) rep(x, 2))))
df_cov_sex$Sex <- factor(df_cov_sex$Sex, levels = c('male', 'female'))
pl_a <- ggplot(df_cov_age, aes(x = gr, y = Age, fill = gr))+
  geom_violin(alpha = 0.8)+
  geom_boxplot(width=0.2, fill="white")+
  xlab('')+ ylab('Age')+
  scale_fill_manual(values = gr_color)+
  annotate("text", x = 1, y = max(df_cov_age$Age)+2, label = sprintf('p=%s', as.character(round(res_cl$test_cov$pval[res_cl$test_cov$cov_id == 'Age'], digits = 2))))+
  theme_bw()+theme(legend.position = 'none')

pl_s <- ggplot(df_cov_sex, aes(x = gr, y = n, color = gr, fill = Sex))+
  geom_bar(size = 1, alpha = 0.8, stat = 'identity', position = position_dodge())+
  xlab('')+ ylab('number of individual')+
  scale_color_manual(values = gr_color)+
  scale_fill_manual(values = c('grey10', 'grey60'))+
  guides(color = FALSE)+
  annotate("text", x = 1, y = max(df_cov_sex$n)+2, label = sprintf('p=%s', as.character(round(res_cl$test_cov$pval[res_cl$test_cov$cov_id == 'Gender'], digits = 2))))+
  theme_bw()+theme(legend.position = 'right')

tot_pl <- ggarrange(plotlist = list(pl_a, pl_s), ncol = 2, nrow = 1, align='h', widths=c(1, 1.4))
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_Age_Sex.png', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 6.5, height = 3, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_Age_Sex.pdf', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 6.5, height = 3, plot = tot_pl, device = 'pdf')


#################
# color based on best genes
tscore_res <- get(load(tscore_featRelFile))
test_feat <- tscore_res$test_feat
test_feat <- test_feat[order(abs(test_feat$estimates), decreasing = T), ]
genes_plot <- sapply(1:P, function(x) test_feat$feat[test_feat$comp == sprintf('gr%i_vs_all',x)][1])
genes_plot <- unique(genes_plot)
if(all(c('CELSR2', 'SORT1') %in% genes_plot)){
  genes_plot <- genes_plot[!genes_plot %in% 'CELSR2']
}

list_pl <- list()
for(i in 1:length(genes_plot)){

  df <- data.frame(component_1=umap_res$layout[,1], component_2=umap_res$layout[,2], col_gene = tscore_res$inputData[, genes_plot[i]])

  list_pl[[i]] <- ggplot(df, aes(x = component_1, y = component_2, color = col_gene))+
    geom_point(size = 0.05, alpha = 0.8)+
    xlab('UMAP component 1')+ ylab('UMAP component 2')+
    theme_bw()+theme(legend.position = 'bottom')+
    scale_colour_gradient2()+
    labs(color=genes_plot[i]) 
}

tot_pl <- ggarrange(plotlist = list_pl, nrow = 1, align='h')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_umap_geneSpec.png', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 2 + 3*(length(genes_plot)-1), height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_%smetric_umap_geneSpec.pdf', outFold, type_cluster_data, type_input, type_cluster, type_sim), width = 2 + 3*(length(genes_plot)-1), height = 4, plot = tot_pl, device = 'pdf')





