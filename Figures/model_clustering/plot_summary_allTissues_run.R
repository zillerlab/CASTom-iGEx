options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="overal results clusters all tissues")
parser$add_argument("--clust_res", type = "character", nargs = '*', help = "")
parser$add_argument("--tissues_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--type_input", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
clust_res <- args$clust_res
tissues_name <- args$tissues_name
pheno_name <- args$pheno_name
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
color_tissues_file <- args$color_tissues_file
outFold <- args$outFold

########################################################################################################################
# tissues_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# clust_res <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissues_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/'
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# pheno_name <- 'CAD_HARD'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
######################################################################################################################

NMI_mat <- matrix(0, ncol = length(tissues_name), nrow = length(tissues_name))
df_cluster <- data.frame(tissue = tissues_name, n_cluster = NA, mod = NA)
df_div <- data.frame(tissue = c(), n_samples = c(), gr = c() )

for(i in 1:length(tissues_name)){
  
  print(i)
  tmp <- get(load(clust_res[i]))
  cl1 <- tmp$cl_best$gr
  df_cluster$n_cluster[i] <- tmp$info_tune$n_gr
  df_cluster$mod[i] <- tmp$info_tune$mod
  df_div <- rbind(df_div, data.frame(tissue = rep(tissues_name[i], tmp$info_tune$n_gr), n_samples = sapply(sort(unique(cl1)), function(x) sum(cl1 == x)), gr = sort(unique(cl1))))
  
  if(i<length(tissues_name)){
    for(j in (i+1):length(tissues_name)){
      tmp <- get(load(clust_res[j]))
      cl2 <- tmp$cl_best$gr
      NMI_mat[i, j] <- compare(cl1, cl2, method = 'nmi')
    }
  }
}

rownames(NMI_mat) <- colnames(NMI_mat) <- tissues_name
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)

# correlation
NMI_mat <- NMI_mat + t(NMI_mat)
diag(NMI_mat) <- NA

col <- colorRampPalette(brewer.pal(9, 'YlGnBu'))(100)
ord <- corrMatOrder(NMI_mat, order="hclust", hclust.method = 'ward.D')
newcolours <- color_tissues$color[match(tissues_name, color_tissues$tissues)][ord]
title_pl <- sprintf('NMI %s %s', type_data, type_input)

pdf(file = sprintf('%s/NMI_%s_%s.pdf', outFold, type_data, type_input), width = 9, height = 6, compress = F, pointsize = 12)
corrplot(NMI_mat, type="upper", order = 'hclust', hclust.method = 'ward.D',
         tl.col = newcolours, tl.cex=1.2,
         col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
         addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.9, mar = c(0,0,1,0))
dev.off()
png(file = sprintf('%s/NMI_%s_%s.png', outFold, type_data, type_input), width = 9, height = 6, res = 300, units = 'in')
corrplot(NMI_mat, type="upper", order = 'hclust', hclust.method = 'ward.D',
         tl.col = newcolours, tl.cex=1.2,
         col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
         addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.9, mar = c(0,0,1,0))
dev.off()

# plot clustering number distribution
df_div$gr <- paste0('gr',df_div$gr)
# df_div$gr <- factor(df_div$gr, levels = rev(unique(sort(df_div$gr))))
df_div$gr <- factor(df_div$gr)
df_div$tissue <- factor(df_div$tissue, levels = tissues_name)
newcolours <- color_tissues$color[match(tissues_name, color_tissues$tissues)]
# color_gr <- rev(pal_d3("category10")(length(levels(df_div$gr))))
df_cluster$tissue <- factor(df_cluster$tissue, levels = tissues_name)

pl <- ggplot(data = df_div, aes(x = tissue, y = n_samples, fill = gr))+
  geom_bar(alpha = 0.7, width = 0.8, stat = 'identity')+
  ylab('n. Cases')+ 
  theme_bw()+ 
  theme(legend.position = 'bottom', 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours))+
  scale_fill_d3()+
  coord_flip() + scale_y_reverse() 

pl_side <- ggplot(data = df_cluster, aes(x = tissue, y = mod))+
  geom_bar(alpha = 0.7, width = 0.8, stat = 'identity', fill = 'grey20')+
  ylab('Modularity')+ 
  theme_classic()+ 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
  coord_flip()

tot_pl <- ggarrange(plotlist = list(pl, pl_side), ncol=2, nrow=1, widths=c(1, 0.3), align = 'h', common.legend = TRUE)

ggsave(filename =  sprintf('%scluster_nsamples_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 6, height = 3, dpi = 500)
ggsave(filename = sprintf('%scluster_nsamples_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 6, height = 3, dpi = 500, compress = F)

df_div$frac <- df_div$n_samples/sum(df_div$n_samples[df_div$tissue == 'Whole_Blood'])

# save results
res <- list(NMI = NMI_mat, cl_div = df_div, cl_nsamples = df_cluster)
save(res, file = sprintf('%scluster_summary_allTissues_%s_%s_cluster%s.RData', outFold, type_data, type_input, type_cluster))

