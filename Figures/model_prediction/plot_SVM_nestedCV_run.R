# plot barplot SVM results (combined tissue)

options(stringsAsFactors=F)
options(max.print=1000)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pROC))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="plot SVM results")
parser$add_argument("--type_data", type = "character", help = "tscore path_Reactome and/or path_GO")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "input fold, 1 per pval thr")
parser$add_argument("--pval_thr", type = "double", nargs = '*', help = "pvalue threshold used")
parser$add_argument("--pheno_name", type = "character", help = "name phenotype")
parser$add_argument("--outFold", type = "character", help = "output folder")

args <- parser$parse_args()
type_data <- args$type_data
inputFile <- args$inputFold
pval_thr <- args$pval_thr
color_tissue_file <- args$color_tissue_file
pheno_name <- args$pheno_name
outFold <- args$outFold

# ###################################################################################################################################
# type_data <- 'path_Reactome'
# pheno_name <- 'CAD-HARD'
# pval_thr <- c(0.01, 0.05, 0.1, 0.2, 1)
# inputFile <- paste0('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_subset/path_Reactome_SNFUMAPp_SVM_FDRpval', pval_thr,'.RData')
# outFold <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_subset/path_Reactome_SNFUMAPp_SVM_'
# ###################################################################################################################################

# spec_sens_funct <- function(true, pred){
#   val_t <- sort(unique(as.numeric(as.character(true))))
#   val_p <- sort(unique(as.numeric(as.character(pred))))
#   min_t <- val_t[1]
#   min_p <- val_p[1]
#   max_t <- val_t[2]
#   max_p <- val_p[2]
#   
#   true <- as.numeric(as.character(true))
#   pred <- as.numeric(as.character(pred))
#   
#   TN <- length(which(true == min_t & pred == min_p))
#   N <- length(which(true == min_t))
#   
#   TP <- length(which(true == max_p & pred == max_p))
#   P <- length(which(true == max_t))
#   
#   sens <- TP/P
#   spec <- TN/N
#   
#   # tab <- table(true, pred)
#   # sens <- tab[2,2]/sum(tab[2,])
#   # spec <- tab[1,1]/sum(tab[1,])
#   return(c(sens,spec))
# }

df <- data.frame()
for(j in 1:length(pval_thr)){
      p <- pval_thr[j]
      print(p)
      tmp <- get(load(inputFile[j]))
      df_tmp <- data.frame(auc_val_mean = mean(tmp$SVM$performance_nestedCV$auc_val), auc_val_sd = sd(tmp$SVM$performance_nestedCV$auc_val), 
                           bacc_val_mean =mean(tmp$SVM$performance_nestedCV$bacc_val), bacc_val_sd = sd(tmp$SVM$performance_nestedCV$bacc_val), 
                           nfeat = nrow(tmp$feat), pval_thr = p)
      df <- rbind(df, df_tmp)
}
  
# save 
write.table(x = df, file = sprintf('%sperformance_SVM.txt', outFold), quote = F, sep = '\t', col.names = T, row.names = F)

df$pval_thr <- factor(df$pval_thr, levels = pval_thr)
df$pval_col <- -log10(as.numeric(as.character(df$pval_thr)))
  
pl_val <- ggplot(data = df, aes(x = pval_thr, y = auc_val_mean, fill = pval_col))+
  geom_bar(stat = 'identity',  position=position_dodge(), color = 'black', alpha = 0.9, width = 0.7)+
  geom_errorbar(aes(ymin = auc_val_mean - auc_val_sd, ymax = auc_val_mean + auc_val_sd), position=position_dodge(.9), width = 0.2)+
  ylab('AUC\nnested CV')+ylim(0,1)+xlab('pvalue threshold')+
  geom_hline(yintercept = 0.5, size = 0.5, linetype = 5, color = 'black')+
  theme_bw()+theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient(low="blue", high="#f8ebd9")

pl_nfeat <- ggplot(data = df, aes(x = pval_thr, y = nfeat, fill = pval_col))+
  geom_bar(stat = 'identity',  position=position_dodge(), color = 'black', alpha = 0.9, width = 0.7)+
  ylab('number\nof features')+ggtitle(paste(pheno_name, type_data))+
  scale_fill_gradient(low="blue", high="#f8ebd9")+
  theme_bw() + theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                     strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(hjust = 0.5))
tot_pl <- ggarrange(plotlist = list(pl_nfeat, pl_val), heights = c(1,1.7), ncol = 1, nrow = 2, align = "v")
ggsave(filename = sprintf('%sauc.png', outFold), plot = tot_pl, width = 4, height = 5)


pl_val <- ggplot(data = df, aes(x = pval_thr, y = bacc_val_mean, fill = pval_col))+
  geom_bar(stat = 'identity',  position=position_dodge(), color = 'black', alpha = 0.9, width = 0.7)+
  geom_errorbar(aes(ymin = bacc_val_mean - bacc_val_sd, ymax = bacc_val_mean + bacc_val_sd), position=position_dodge(.9), width = 0.2)+
  ylab('Balanced Accuracy\nnested CV')+ylim(0,1)+xlab('pvalue threshold')+
  geom_hline(yintercept = 0.5, size = 0.5, linetype = 5, color = 'black')+
  theme_bw()+theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient(low="blue", high="#f8ebd9")

pl_nfeat <- ggplot(data = df, aes(x = pval_thr, y = nfeat, fill = pval_col))+
  geom_bar(stat = 'identity',  position=position_dodge(), color = 'black', alpha = 0.9, width = 0.7)+
  ylab('number\nof features')+ggtitle(paste(pheno_name, type_data))+
  scale_fill_gradient(low="blue", high="#f8ebd9")+
  theme_bw() + theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                     strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(hjust = 0.5))
tot_pl <- ggarrange(plotlist = list(pl_nfeat, pl_val), heights = c(1,1.7), ncol = 1, nrow = 2, align = "v")
ggsave(filename = sprintf('%sbacc.png', outFold), plot = tot_pl, width = 4, height = 5)



