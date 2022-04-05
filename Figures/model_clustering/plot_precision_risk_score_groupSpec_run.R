# plot fraction of correct match per tissue (precision)

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(e1071))
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

parser <- ArgumentParser(description="precision risk-score group specific, plot")
parser$add_argument("--riskScore_comp_file", nargs = '*', type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--tissues", nargs = '*', type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
riskScore_comp_file <- args$riskScore_comp_file
tissues <- args$tissues
color_tissues_file <- args$color_tissues_file
outFold <- args$outFold

###################################################################################################################
# tissues <- c('Liver', 'Heart_Left_Ventricle', 'Colon_Sigmoid')
# riskScore_comp_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissues,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/riskScores_clusterCases_CAD_HARD_group_relatedPheno_measureGoodnessPred.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
#################################################################################################################

color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues,color_tissues$tissue),]

df_tot <- list()
for(i in 1:length(tissues)){
  
  tmp <- read.delim(riskScore_comp_file[i], h=T, stringsAsFactors = F, sep = '\t')
  tmp$tissue <- tissues[i]
  df_tot[[i]] <- tmp
}
df_tot <- do.call(rbind, df_tot)

df_tot$tissue <- factor(df_tot$tissue, levels = tissues)
df_tot <- df_tot[df_tot$rs_sign, ]

thr_measure <- seq(min(df_tot$measure),max(df_tot$measure), length.out = 1000)
df_repr <- list()

for(j in 1:length(thr_measure)){
  # print(j)
  df_repr[[j]] <- data.frame(tissue = tissues, n_pheno_thr = NA, 
                             n_same_sign = NA, precision = NA, recall = NA, 
                             thr_measure = thr_measure[j])
  
  for(i in 1:length(tissues)){
    
    tmp <- df_tot[df_tot$tissue == tissues[i],]
    df_repr[[j]]$n_pheno_thr[i] <- sum(tmp$measure >= thr_measure[j])
    df_repr[[j]]$n_same_sign[i] <- sum(tmp$same_sign == 1)
    
    if(sum(tmp$measure >= thr_measure[j]) > 0){
      df_repr[[j]]$precision[i] <- sum(tmp$measure >= thr_measure[j] & tmp$same_sign == 1)/sum(tmp$measure >= thr_measure[j])
    }else{
      df_repr[[j]]$precision[i] <- NA
    }
    
    if(sum(tmp$same_sign == 1) > 0){
      df_repr[[j]]$recall[i] <- sum(tmp$measure >= thr_measure[j] & tmp$same_sign == 1)/sum(tmp$same_sign == 1)
    }else{
      df_repr[[j]]$recall[i] <- NA
    }
  }
  
}

df_repr <- do.call(rbind, df_repr)
df_repr$tissue <- factor(df_repr$tissue, levels = tissues)

# thr <- max(sapply(tissues, function(x) min(df_repr$thr_measure[df_repr$precision == 1 & df_repr$tissue == x], na.rm  = T)))
thr <- sapply(tissues, function(x) min(df_repr$thr_measure[df_repr$precision ==1 & df_repr$tissue == x], 
                                           na.rm  = T))
if(any(!is.finite(thr))){
  tissue_rep <- tissues[!is.finite(thr)]
  new_thr <- sapply(tissue_rep, function(x) min(df_repr$thr_measure[df_repr$precision == 0 & df_repr$tissue == x], 
                                         na.rm  = T))
  thr <- c(new_thr, thr[is.finite(thr)])
}

thr <- max(thr)
print(paste0('THR CRM-measure for plot: ',thr))


tmp <- subset(df_repr, thr_measure <= thr)
pl1 <- ggplot(data =  tmp,aes(x = thr_measure, y = precision, color = tissue, group = tissue))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('F-stat * |beta|')+ylab('Precision\nsame beta sign')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''))+
  scale_color_manual(values = color_tissues$color)

pl2 <- ggplot(data =  tmp,aes(x = thr_measure, y = n_pheno_thr, color = tissue, group = tissue))+
  geom_point(alpha = 0.7)+
  geom_line(alpha = 0.7, size = 0.8)+
  xlab('F-stat * |beta|')+ylab('n. phenotypes')+ 
  theme_bw()+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), 
                   legend.text = element_text(size = 7),  legend.title = element_text(size = 8), legend.key.size = unit(0.3, "cm"))+
  guides(color=guide_legend(title = ''))+
  scale_color_manual(values = color_tissues$color)

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h', common.legend = T, legend = 'right')
ggsave(filename = sprintf('%sriskScores_clusterCases_group_relatedPheno_npheno_and_precision.png', outFold), plot = tot_pl, width = 9, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%sriskScores_clusterCases_group_relatedPheno_npheno_and_precision.pdf', outFold), plot = tot_pl, width = 9, height = 4.5, dpi = 500, compress = F)

write.table(file = sprintf('%sriskScores_clusterCases_group_relatedPheno_npheno_and_precision.txt', outFold), 
            x = df_repr, col.names = T, row.names = F, sep = '\t', quote = F)





