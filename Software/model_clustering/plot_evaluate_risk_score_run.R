#!/usr/bin/env Rscript

# plot risk score evaluation metric

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="R2 from phenotype - risk score, plot")
parser$add_argument("--riskScore_eval_file", nargs = '*', type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--tissues", nargs = '*', type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
riskScore_eval_file <- args$riskScore_eval_file
tissues <- args$tissues
color_tissues_file <- args$color_tissues_file
outFold <- args$outFold

###################################################################################################################
# tissues <- c('Liver', 'Heart_Left_Ventricle')
# riskScore_eval_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissues,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_corrThr0.5_relatedPhenotypes_R2_risk_score_phenotype.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
#################################################################################################################

color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues,color_tissues$tissue),]

df_tot <- list()
for(i in 1:length(tissues)){
  
  tmp <- read.delim(riskScore_eval_file[i], h=T, stringsAsFactors = F, sep = '\t')
  tmp$tissue <- tissues[i]
  df_tot[[i]] <- tmp
}
df_tot <- do.call(rbind, df_tot)

df_tot$tissue <- factor(df_tot$tissue, levels = tissues)

pl1 <- ggplot(data =  df_tot, aes(x = nsamples, y = R2_risk, color = tissue, group = tissue))+
  geom_point(alpha = 0.7, size = 0.3)+
  geom_line(alpha = 0.7, size = 0.3)+
  xlab('n. samples')+ylab('R2 risk-score')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''))+
  scale_color_manual(values = color_tissues$color)

pl2 <- ggplot(data =  df_tot, aes(x = nsamples, y = Fstat_risk, color = tissue, group = tissue))+
  geom_point(alpha = 0.7, size = 0.3)+
  geom_line(alpha = 0.7, size = 0.3)+
  xlab('n. samples')+ylab('F-statistic risk-score')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''))+
  scale_color_manual(values = color_tissues$color)

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h', common.legend = T, legend = 'right')
ggsave(filename =  sprintf('%splot_R2Fstat_nsamples_riskScore.png', outFold), plot = tot_pl, width = 8, height = 4, dpi = 500)
ggsave(filename = sprintf('%splot_R2Fstat_nsamples_riskScore.pdf', outFold), plot = tot_pl, width = 8, height = 4, dpi = 500, compress = F)



