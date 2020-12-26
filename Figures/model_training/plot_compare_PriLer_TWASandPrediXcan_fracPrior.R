# plots comparison PriLer vs prediXcan and TWAS, fraction of reg snps with prior
# on denbi
# instead of recomputing intersection with GREsbased on TWAS and prediXcan SNPs, consider already done prior matrix from PriLer and intersect with that
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(latex2exp))

tissues_model <- read.csv('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/final_model_gtex.csv', h=F, stringsAsFactors = F)
outFold <- '/mnt/lucia/PriLer_TRAIN_PLOTS/'
colnames(tissues_model) <- c('tissue', 'type')
tissues_model$folder_train <- sapply(tissues_model$tissue, function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/noGWAS/', x))
tissues_model$folder_train[tissues_model$type == 'CAD'] <- sapply(tissues_model$tissue[tissues_model$type == 'CAD'] ,
                                                                  function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/CAD_GWAS_bin5e-2/', x))
tissues_model$folder_train[tissues_model$type == 'PGC'] <- sapply(tissues_model$tissue[tissues_model$type == 'PGC'], 
                                                                  function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/PGC_GWAS_bin1e-2/', x))

prior_names <- read.csv('/mnt/lucia/PriLer_TRAIN_PLOTS/prior_association_TRAIN.csv', h=F, stringsAsFactors = F)
rownames(prior_names) <- prior_names$V1
prior_names <- prior_names[,-1]
prior_names <- lapply(tissues_model$tissue, function(x) prior_names[rownames(prior_names) == x,][prior_names[rownames(prior_names) == x,] != ''])
names(prior_names) <- tissues_model$tissue

# load reg-SNP TWAS and prediXcan
prediXcan_snps <- read.table(gzfile('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/AllTissues/200kb/noGWAS/prediXcan_regSNPs_annotation.txt.gz'), h=T, stringsAsFactors = F, sep = '\t', check.names = F)
prediXcan_new_id <- paste0(prediXcan_snps$chrom, '_',prediXcan_snps$position) # no duplication
prediXcan_snps <- cbind(prediXcan_snps[, 1:3], prediXcan_snps[, match(names(prior_names),colnames(prediXcan_snps))])
TWAS_snps <- read.table(gzfile('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/AllTissues/200kb/noGWAS/TWAS_regSNPs_annotation.txt.gz'), h=T, stringsAsFactors = F, sep = '\t',  check.names = F)
TWAS_new_id <- paste0(TWAS_snps$chrom, '_',TWAS_snps$position) # no diplication
TWAS_snps <- cbind(TWAS_snps[, 1:3], TWAS_snps[, match(names(prior_names),colnames(TWAS_snps))])

# intersect with prior for model
chr <- 1:22
prediXcan_tissue_chr <- matrix(nrow = 22, ncol = length(prior_names))
prediXcan_prior_tissue_chr <- matrix(nrow = 22, ncol = length(prior_names))
TWAS_tissue_chr <- matrix(nrow = 22, ncol = length(prior_names))
TWAS_prior_tissue_chr <- matrix(nrow = 22, ncol = length(prior_names))
for(i in chr){
  
  print(i)
  snp_info <- read.table(sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/hg19_SNPs_chr%i_matched.txt', i), h=T, stringsAsFactors = F, sep = '\t')
  id_pos <- paste0(snp_info$chrom , '_', snp_info$position)
  priorMat <- read.table(gzfile(sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/priorMatrix_chr%i.txt.gz', i)), h=T, stringsAsFactors = F, sep = '\t')
  # remove duplicated positions
  priorMat <- priorMat[!duplicated(id_pos), ]
  snp_info <- snp_info[!duplicated(id_pos), ]
  
  # intersect with prediXcan
  id <- intersect(paste0(snp_info$chrom , '_', snp_info$position), prediXcan_new_id)
  id_pos <- paste0(snp_info$chrom , '_', snp_info$position)
  tmp_info <- snp_info[match(id, id_pos),]
  tmp_mat <- priorMat[match(id, id_pos),]
  
  tmp_prediXcan <- prediXcan_snps[match(id, prediXcan_new_id), ]
  prediXcan_tissue_chr[i, ] <- colSums(tmp_prediXcan[, -(1:3)]!=0)
  prediXcan_prior_tissue_chr[i, ] <- sapply(1:length(prior_names), function(x) sum((tmp_prediXcan[, names(prior_names)[x]] != 0) & (rowSums(tmp_mat[, prior_names[[x]], drop = F] !=0)>0 )))
    
  
  # intersect with TWAS
  id <- intersect(paste0(snp_info$chrom , '_', snp_info$position), TWAS_new_id)
  id_pos <- paste0(snp_info$chrom , '_', snp_info$position)
  tmp_info <- snp_info[match(id, id_pos),]
  tmp_mat <- priorMat[match(id, id_pos),]
  
  tmp_TWAS <- TWAS_snps[match(id, TWAS_new_id), ]
  TWAS_tissue_chr[i, ] <- colSums(tmp_TWAS[, -(1:3)]!=0)
  TWAS_prior_tissue_chr[i, ] <- sapply(1:length(prior_names), function(x) sum((tmp_TWAS[, names(prior_names)[x]] != 0) & (rowSums(tmp_mat[, prior_names[[x]], drop = F] !=0)>0 )))
  
}

df_prediXcan <- data.frame(tissue = names(prior_names), nregSNPS = colSums(prediXcan_tissue_chr), nregSNPs_with_prior = colSums(prediXcan_prior_tissue_chr))
df_prediXcan$frac <- df_prediXcan$nregSNPs_with_prior/df_prediXcan$nregSNPS
df_TWAS <- data.frame(tissue = names(prior_names), nregSNPS = colSums(TWAS_tissue_chr), nregSNPs_with_prior = colSums(TWAS_prior_tissue_chr))
df_TWAS$frac <- df_TWAS$nregSNPs_with_prior/df_TWAS$nregSNPS
write.table(df_TWAS, file = '/mnt/lucia/PriLer_TRAIN_PLOTS/TWAS_fracPrior_intersect_PriLer.txt', quote = F, sep = '\t', col.names = T, row.names = F)
write.table(df_prediXcan, file = '/mnt/lucia/PriLer_TRAIN_PLOTS/prediXcan_fracPrior_intersect_PriLer.txt', quote = F, sep = '\t', col.names = T, row.names = F)
df_PriLer <- read.delim('/mnt/lucia/PriLer_TRAIN_PLOTS/compare_PriLer_elnet.txt', h=T, stringsAsFactors = F, sep = '\t')

### plot
df_snps <- data.frame(tissue = rep(names(prior_names), 3), type = c(rep('PriLer', length(prior_names)), rep('prediXcan', length(prior_names)), rep('TWAS', length(prior_names))), 
                      frac = c(df_PriLer$frac_prior[match(names(prior_names), df_PriLer$tissue)], df_prediXcan$frac, df_TWAS$frac))
df_snps$tissue <- factor(df_snps$tissue, levels = names(prior_names))
df_snps$type <- factor(df_snps$type, levels = c('PriLer', 'TWAS', 'prediXcan'))
color_tissues <- c(rep('#999900',2),'#3CB371', rep('#CD5C5C', 3), rep('#4169E1', 8), '#FF8C00',  rep('#8B4513', 5), rep('#DC143C', 2),   '#8A2BE2','#708090','#A52A2A',
                   '#3CB371', rep('#FA8072', 2), '#8B4513',  '#3CB371', '#8B4513', '#3CB371', '#FF8C00')
color_type <- c('#1F77B4FF', '#FF7F0EFF', '#D62728FF')
outFold <- '/mnt/lucia/PriLer_TRAIN_PLOTS/'

pl_frac <- ggplot(df_snps, aes(x = tissue, y = frac, fill = type, group = type)) + 
  geom_bar(alpha = 0.7, stat = 'identity', width = 0.8, position = position_dodge()) + ggtitle(expression(frac("n.reg-SNPs with prior","n.reg-SNPs")))+
  ylab('fraction')+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.y = element_blank(), legend.title = element_blank(), 
        axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
        axis.text.y=element_text(size = 10, colour = color_tissues), legend.position = 'bottom', plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"))+
  scale_fill_manual(values = color_type)+
  coord_flip()

ggsave(filename = sprintf('%sFracRegSNPs_PriLerVSTWASandPrediXcan_AllTissues.png', outFold), plot = pl_frac, width = 5.5, height = 10, dpi = 500)
ggsave(filename = sprintf('%sFracRegSNPs_PriLerVSTWASandPrediXcan_AllTissues.pdf', outFold), plot = pl_frac, width = 5.5, height = 10, dpi = 500)

print(mean(df_snps$frac[df_snps$type == 'PriLer'] - df_snps$frac[df_snps$type == 'TWAS'])) # 0.06610335
print(sd(df_snps$frac[df_snps$type == 'PriLer'] - df_snps$frac[df_snps$type == 'TWAS'])) # 0.03147366
print(mean(df_snps$frac[df_snps$type == 'PriLer'] - df_snps$frac[df_snps$type == 'prediXcan'])) # 0.1014615
sd(df_snps$frac[df_snps$type == 'PriLer'] - df_snps$frac[df_snps$type == 'prediXcan']) # 0.03323968


