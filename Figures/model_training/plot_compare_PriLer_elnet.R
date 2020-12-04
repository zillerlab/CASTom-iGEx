# plots on denbi, contain both CMC and GTEx res
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(latex2exp)

tissues_model <- read.csv('/mnt/lucia/PriLer_PROJECT_GTEx/INPUT_DATA/final_model_gtex.csv', h=F, stringsAsFactors = F)
outFold <- '/mnt/lucia/PriLer_TRAIN_PLOTS/'
colnames(tissues_model) <- c('tissue', 'type')
tissues_model$folder_train <- sapply(tissues_model$tissue, function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/noGWAS/', x))
tissues_model$folder_train[tissues_model$type == 'CAD'] <- sapply(tissues_model$tissue[tissues_model$type == 'CAD'] ,
                                                                  function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/CAD_GWAS_bin5e-2/', x))
tissues_model$folder_train[tissues_model$type == 'PGC'] <- sapply(tissues_model$tissue[tissues_model$type == 'PGC'], 
                                                                  function(x) sprintf('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/%s/200kb/PGC_GWAS_bin1e-2/', x))


type <- c('Control50', 'Control100', 'Control150', 'ControlAll', 'All')
res_CMC <- paste0('/mnt/lucia/PriLer_TRAIN_PLOTS/res_CMC/train_', type)

df_prior <- list()
df_noprior <- list()
df_both <- list()
df_regSNPs <- list()
df_nrel <- list()

df_sampleSize <- read.table('/mnt/lucia/PriLer_TRAIN_PLOTS/res_CMC/nSamples_nRelGenes.txt', h=T, stringsAsFactors=F)
df_sampleSize <- df_sampleSize[df_sampleSize$type %in% 'All',]
df_sampleSize$type <- 'DLPC_CMC'
df_sampleSize <- rbind(df_sampleSize, read.table('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/nSamples_nRelGenes.txt', h=T, stringsAsFactors=F))

# load CMC (train_All) res 
df_prior[[1]] <- read.delim(sprintf('%s/resPrior_regEval_allchr.txt', res_CMC[5]), h=T, stringsAsFactors = F, sep = '\t')
df_prior[[1]]$tissue <- 'DLPC_CMC'
df_noprior[[1]] <- read.delim(sprintf('%s/resNoPrior_regEval_allchr.txt', res_CMC[5]), h=T, stringsAsFactors = F, sep = '\t')
df_noprior[[1]]$tissue <- 'DLPC_CMC'
identical(df_noprior[[1]]$ensembl_gene_id, df_prior[[1]]$ensembl_gene_id)
df_nrel[[1]] <- data.frame(tissue = 'DLPC_CMC', n_samples = df_sampleSize$n_samples[df_sampleSize$type == 'DLPC_CMC'], prior = length(which(df_prior[[1]]$test_dev_geno>0 & df_prior[[1]]$dev_geno>0.01)),
                           noprior = length(which(df_noprior[[1]]$test_dev_geno>0 & df_noprior[[1]]$dev_geno>0.01)))
df_both[[1]] <- data.frame(geneID = df_prior[[1]]$ensembl_gene_id, 
                           dev_geno_prior = df_prior[[1]]$dev_geno, dev_geno_noprior = df_noprior[[1]]$dev_geno,
                           test_dev_geno_prior = df_prior[[1]]$test_dev_geno, test_dev_geno_noprior = df_noprior[[1]]$test_dev_geno,
                           tissue = df_prior[[1]]$tissue, stringsAsFactors = F)
# only reliable genes for prior
df_both[[1]] <- df_both[[1]][df_both[[1]]$test_dev_geno_prior>0 & df_both[[1]]$dev_geno_prior>0.01, ]
# # avoid big values due to the division
# df_both[[1]] <- df_both[[1]][df_both[[1]]$dev_geno_noprior >= 10^-3, ]
df_both[[1]]$increase <- (df_both[[1]]$dev_geno_prior - df_both[[1]]$dev_geno_noprior)/df_both[[1]]$dev_geno_noprior


tmp <- read.delim(sprintf('%s/nVariants_prior_compare.txt', res_CMC[5]), h=T, stringsAsFactors = F, sep = '\t')
df_regSNPs[[1]] <- data.frame(tissue = 'DLPC_CMC', nregSNPS_prior = sum(tmp$len[tmp$type == 'el-net learned prior']), 
                              nregSNPS_noprior = sum(tmp$len[tmp$type == 'el-net']),stringsAsFactors = F)
df_regSNPs[[1]]$frac_prior <- sum(tmp$len[tmp$type == 'el-net learned prior' & tmp$type_snp %in% c('with prior (common)', 'with prior (unique)')])/df_regSNPs[[1]]$nregSNPS_prior
df_regSNPs[[1]]$frac_noprior <- sum(tmp$len[tmp$type == 'el-net' & tmp$type_snp %in% c('with prior (common)', 'with prior (unique)')])/df_regSNPs[[1]]$nregSNPS_noprior

# load GTEx res 
for(i in 1:nrow(tissues_model)){
  
  print(i)
  df_prior[[i+1]] <- read.delim(sprintf('%s/resPrior_regEval_allchr.txt',  tissues_model$folder_train[i]), h=T, stringsAsFactors = F, sep = '\t')
  df_prior[[i+1]]$tissue <- tissues_model$tissue[i]
  df_noprior[[i+1]] <- read.delim(sprintf('%s/resNoPrior_regEval_allchr.txt',  tissues_model$folder_train[i]), h=T, stringsAsFactors = F, sep = '\t')
  df_noprior[[i+1]]$tissue <- tissues_model$tissue[i]
  print(identical(df_noprior[[i+1]]$ensembl_gene_id, df_prior[[i+1]]$ensembl_gene_id))
  
  df_nrel[[i+1]] <- data.frame(tissue = tissues_model$tissue[i], n_samples = df_sampleSize$n_samples[df_sampleSize$type == tissues_model$tissue[i]], 
                               prior = length(which(df_prior[[i+1]]$test_dev_geno>0 & df_prior[[i+1]]$dev_geno>0.01)),
                               noprior = length(which(df_noprior[[i+1]]$test_dev_geno>0 & df_noprior[[i+1]]$dev_geno>0.01)))
  df_both[[i+1]] <- data.frame(geneID = df_prior[[i+1]]$ensembl_gene_id, 
                               dev_geno_prior = df_prior[[i+1]]$dev_geno, dev_geno_noprior = df_noprior[[i+1]]$dev_geno,
                               test_dev_geno_prior = df_prior[[i+1]]$test_dev_geno, test_dev_geno_noprior = df_noprior[[i+1]]$test_dev_geno,
                               tissue = df_prior[[i+1]]$tissue, stringsAsFactors = F)
  # only reliable genes for prior
  df_both[[i+1]] <- df_both[[i+1]][df_both[[i+1]]$test_dev_geno_prior>0 & df_both[[i+1]]$dev_geno_prior>0.01, ]
  # # avoid big values due to the division
  # df_both[[i+1]] <- df_both[[i+1]][df_both[[i+1]]$dev_geno_noprior >= 10^-3, ]
  df_both[[i+1]]$increase <- (df_both[[i+1]]$dev_geno_prior - df_both[[i+1]]$dev_geno_noprior)/df_both[[i+1]]$dev_geno_noprior
  
  tmp <- read.delim(sprintf('%s/nVariants_prior_compare.txt', tissues_model$folder_train[i]), h=T, stringsAsFactors = F, sep = '\t')
  df_regSNPs[[i+1]] <- data.frame(tissue = tissues_model$tissue[i], nregSNPS_prior = sum(tmp$len[tmp$type == 'el-net learned prior']), 
                                  nregSNPS_noprior = sum(tmp$len[tmp$type == 'el-net']),stringsAsFactors = F)
  df_regSNPs[[i+1]]$frac_prior <- sum(tmp$len[tmp$type == 'el-net learned prior' & tmp$type_snp %in% c('with prior (common)', 'with prior (unique)')])/df_regSNPs[[i+1]]$nregSNPS_prior
  df_regSNPs[[i+1]]$frac_noprior <- sum(tmp$len[tmp$type == 'el-net' & tmp$type_snp %in% c('with prior (common)', 'with prior (unique)')])/df_regSNPs[[i+1]]$nregSNPS_noprior
  
}

color_tissues <- c('#000080', rep('#999900',2),'#3CB371', rep('#CD5C5C', 3), rep('#4169E1', 8), '#FF8C00',  rep('#8B4513', 5), rep('#DC143C', 2),   '#8A2BE2','#708090','#A52A2A',
                   '#3CB371', rep('#FA8072', 2), '#8B4513',  '#3CB371', '#8B4513', '#3CB371', '#FF8C00')
n_relGenes <- data.frame(tissues = c('DLPC_CMC', tissues_model$tissue), col = color_tissues, stringsAsFactors = F)

##############################
### plot n. reliable genes ###
##############################

df_nrel <- do.call(rbind, df_nrel)
df_nrel$tissue <- factor(df_nrel$tissue, levels = n_relGenes$tissues)

pl_nrel <- ggplot(df_nrel, aes(x = noprior, y = prior, color = tissue, size = n_samples)) + 
  geom_point()+
  # xlim(min_v, max_v)+ ylim(min_v, max_v)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.6)+ggtitle('n. reliable genes')+
  xlab('elnet')+ ylab('PriLer')+ theme_classic()+
  geom_text_repel(label = df_nrel$tissue, size = 3, alpha = 0.8, segment.size = 0.2, force = 10)+
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 0, hjust = 1),
        axis.text.y=element_text(size = 11), legend.position = 'bottom')+
  guides(color = FALSE, size = guide_legend(title="n. samples"))+
  scale_color_manual(values = n_relGenes$col)

df_nrel$col <- color_tissues

ggsave(filename = sprintf('%snRelGenes_PriLerVSelnet_AllTissues.png', outFold), plot = pl_nrel, width = 5, height = 5, dpi = 500)
ggsave(filename = sprintf('%snRelGenes_PriLerVSelnet_AllTissues.pdf', outFold), plot = pl_nrel, width = 5, height = 5, dpi = 500)

# cor(df_nrel$n_samples, df_nrel$prior) # 0.8537 

########################
### plot increase R2 ###
########################

df_R2 <- do.call(rbind, df_both)
df_R2$diff <- df_R2$dev_geno_prior - df_R2$dev_geno_noprior
df_R2$ratio <- df_R2$dev_geno_prior/df_R2$dev_geno_noprior
df_R2$diff_test <- df_R2$test_dev_geno_prior - df_R2$test_dev_geno_noprior
df_R2$ratio_test <- df_R2$test_dev_geno_prior/df_R2$test_dev_geno_noprior

n_relGenes <- data.frame(tissues = c('DLPC_CMC', tissues_model$tissue), col = color_tissues, stringsAsFactors = F)
n_relGenes$n_genes <- sapply(df_prior, function(x) length(which(x$test_dev_geno>0 & x$dev_geno>0.01)))
order_t <- order(n_relGenes$n_genes)
df_R2$tissue <- factor(df_R2$tissue, levels = n_relGenes$tissues[order_t])

# pl_box_R2inc <- ggplot(df_R2, aes(x = tissue, y=increase, fill = tissue)) + 
#   geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(frac(R[PriLer]^2 - R[elnet]^2,R[elnet]^2)))+ ggtitle('R2 genotype reliable genes')+
#   xlab('')+ theme_bw()+scale_fill_manual(values = n_relGenes$col)+
#   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
#         axis.text.x=element_text(size = 7, angle = 45, hjust = 1, colour = n_relGenes$col[order_t]),
#         axis.text.y=element_text(size = 10), legend.position = 'none')
# 
# ggsave(filename = sprintf('/mnt/lucia/PriLer/R2geno_increase_allTissues.png'), plot = pl_box_R2inc, width = 8, height = 4, dpi=500)
# 
# # capped version
# pl_box_R2inc_capped <- ggplot(subset(df_R2, increase <=1 & increase >= -1), aes(x = tissue, y=increase, fill = tissue)) + 
#   geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(frac(R[PriLer]^2 - R[elnet]^2,R[elnet]^2)))+ ggtitle('R2 genotype reliable genes')+
#   xlab('')+ theme_bw()+scale_fill_manual(values = n_relGenes$col[order_t])+
#   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
#         axis.text.x=element_text(size = 7, angle = 45, hjust = 1, colour = n_relGenes$col[order_t]),
#         axis.text.y=element_text(size = 10), legend.position = 'none')
# 
# ggsave(filename = sprintf('/mnt/lucia/PriLer/R2geno_increase_capped_allTissues.png'), plot = pl_box_R2inc_capped, width = 8, height = 4, dpi=500)

# difference 
pl_box_R2diff <- ggplot(df_R2, aes(x = tissue, y=diff, fill = tissue)) + 
  geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(R[PriLer]^2 - R[elnet]^2))+ ggtitle('R2 genotype reliable genes')+
  xlab('')+ theme_bw()+ scale_fill_manual(values = n_relGenes$col[order_t])+
  # scale_y_continuous(trans=weird)+
  geom_hline(yintercept = 0, color = 'darkgrey', linetype='dashed')+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
        axis.text.y=element_text(size = 9, angle = 0, hjust = 1, colour = n_relGenes$col[order_t]),
        axis.text.x=element_text(size = 10), legend.position = 'none')+coord_flip()

ggsave(filename = sprintf('%sR2geno_diff_allTissues.png', outFold), plot = pl_box_R2diff, width = 6, height = 7, dpi=500)
ggsave(filename = sprintf('%sR2geno_diff_allTissues.pdf', outFold), plot = pl_box_R2diff, width = 6, height = 7, dpi=500)

# difference (test)
pl_box_R2diff <- ggplot(df_R2, aes(x = tissue, y=diff_test, fill = tissue)) + 
  geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(R[PriLer]^2 - R[elnet]^2))+ ggtitle('CV test R2 genotype reliable genes')+
  xlab('')+ theme_bw()+ scale_fill_manual(values = n_relGenes$col[order_t])+
  geom_hline(yintercept = 0, color = 'darkgrey', linetype='dashed')+
  # scale_y_continuous(trans=weird)+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
        axis.text.y=element_text(size = 9, angle = 0, hjust = 1, colour = n_relGenes$col[order_t]),
        axis.text.x=element_text(size = 10), legend.position = 'none')+coord_flip()

ggsave(filename = sprintf('%sR2geno_diffTest_allTissues.png', outFold), plot = pl_box_R2diff, width = 6, height = 7, dpi=500)
ggsave(filename = sprintf('%sR2geno_diffTest_allTissues.pdf', outFold), plot = pl_box_R2diff, width = 6, height = 7, dpi=500)

# ratio 
pl_box_R2ratio <- ggplot(subset(df_R2, dev_geno_noprior > 0 & dev_geno_prior > 0 & ratio < 2), aes(x = tissue, y=ratio, fill = tissue)) + 
  geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(R[PriLer]^2/R[elnet]^2))+ ggtitle('R2 genotype reliable genes')+
  xlab('')+ theme_bw()+ scale_fill_manual(values = n_relGenes$col[order_t])+
  geom_hline(yintercept = 1, color = 'darkgrey', linetype='dashed')+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
        axis.text.y=element_text(size = 9, angle = 0, hjust = 1, colour = n_relGenes$col[order_t]),
        axis.text.x=element_text(size = 10), legend.position = 'none')+coord_flip()

ggsave(filename = sprintf('%sR2geno_ratio_allTissues.png', outFold), plot = pl_box_R2ratio, width = 6, height = 7, dpi=500)
ggsave(filename = sprintf('%sR2geno_ratio_allTissues.pdf', outFold), plot = pl_box_R2ratio, width = 6, height = 7, dpi=500)

# difference (test)
pl_box_R2ratioTest <- ggplot(subset(df_R2, test_dev_geno_prior > 0 & test_dev_geno_noprior > 0 & ratio_test < 3 & is.finite(test_dev_geno_prior) & is.finite(test_dev_geno_noprior)), aes(x = tissue, y=ratio_test, fill = tissue)) + 
  geom_boxplot(outlier.size = 0.7, alpha = 0.7)+ ylab(expression(R[PriLer]^2/R[elnet]^2))+ ggtitle('CV test R2 reliable genes')+
  xlab('')+ theme_bw()+ scale_fill_manual(values = n_relGenes$col[order_t])+
  geom_hline(yintercept = 1, color = 'darkgrey', linetype='dashed')+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(size = 10),
        axis.text.y=element_text(size = 9, angle = 0, hjust = 1, colour = n_relGenes$col[order_t]),
        axis.text.x=element_text(size = 10), legend.position = 'none')+coord_flip()

ggsave(filename = sprintf('%sR2geno_ratioTest_allTissues.png', outFold), plot = pl_box_R2ratioTest, width = 6, height = 7, dpi=500)
ggsave(filename = sprintf('%sR2geno_ratioTest_allTissues.pdf', outFold), plot = pl_box_R2ratioTest, width = 6, height = 7, dpi=500)


# test for each tissue
df_test_R2 <- data.frame(tissue = n_relGenes$tissue, 
                         wilcoxon_testgeno_pval = sapply(n_relGenes$tissue, function(x) wilcox.test(df_R2$test_dev_geno_prior[df_R2$tissue == x], df_R2$test_dev_geno_noprior[df_R2$tissue == x])$p.value))
df_test_R2$wilcoxon_geno_pval <- sapply(n_relGenes$tissue, function(x) wilcox.test(df_R2$dev_geno_prior[df_R2$tissue == x], df_R2$dev_geno_noprior[df_R2$tissue == x])$p.value)
df_test_R2$binom_testgeno_pval <- sapply(n_relGenes$tissue, function(x) binom.test(x = length(which(df_R2$diff_test[df_R2$tissue == x]>0)),n = length(df_R2$diff_test[df_R2$tissue == x]))$p.value)
df_test_R2$binom_geno_pval <- sapply(n_relGenes$tissue, function(x) binom.test(x = length(which(df_R2$diff[df_R2$tissue == x]>0)),n = length(df_R2$diff[df_R2$tissue == x]))$p.value)
df_test_R2$binom_testgeno_perc <- sapply(n_relGenes$tissue, function(x) sum(df_R2$diff_test[df_R2$tissue == x]>0, na.rm = T)/length(df_R2$diff_test[df_R2$tissue == x]))
df_test_R2$binom_geno_perc <- sapply(n_relGenes$tissue, function(x) sum(df_R2$diff[df_R2$tissue == x]>0, na.rm = T)/length(df_R2$diff[df_R2$tissue == x]))
df_test_R2$col = df_nrel$col
write.table(file = sprintf('%scompare_PriLer_elnet_R2difftest.txt', outFold), df_test_R2, quote = F, col.names = T, row.names = F, sep = '\t')


###############################
### plot number of reg-SNPs ###
###############################

df_snps <- do.call(rbind, df_regSNPs)
df_snps$tissue <- factor(df_snps$tissue, levels = n_relGenes$tissues)
min_v <- min(c(df_snps$nregSNPS_prior, df_snps$nregSNPS_noprior))
max_v <- max(c(df_snps$nregSNPS_prior, df_snps$nregSNPS_noprior))

pl_nreg <- ggplot(df_snps, aes(x = nregSNPS_noprior, y = nregSNPS_prior, color = tissue)) + 
  geom_point()+ xlim(min_v, max_v)+ ylim(min_v, max_v)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.6)+ggtitle('n. reg-SNPs')+
  xlab('elnet')+ ylab('PriLer')+ theme_classic()+
  geom_text_repel(label = df_snps$tissue, size = 3, alpha = 0.8, segment.size = 0.2, force = 10)+
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 0, hjust = 1),
        axis.text.y=element_text(size = 11), legend.position = 'none')+
  scale_color_manual(values = n_relGenes$col)

ggsave(filename = sprintf('%snRegSNPs_PriLerVSelnet_AllTissues.png', outFold), plot = pl_nreg, width = 5, height = 5, dpi = 500)
ggsave(filename = sprintf('%snRegSNPs_PriLerVSelnet_AllTissues.pdf', outFold), plot = pl_nreg, width = 5, height = 5, dpi = 500)

min_v <- min(c(df_snps$frac_noprior, df_snps$frac_prior))
max_v <- max(c(df_snps$frac_noprior, df_snps$frac_prior))

pl_frac <- ggplot(df_snps, aes(x = frac_noprior, y = frac_prior, color = tissue)) + 
  geom_point()+ xlim(min_v, max_v)+ ylim(min_v, max_v)+ggtitle(expression(frac("n.reg-SNPs with prior","n.reg-SNPs")))+
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.6)+
  xlab('elnet')+ ylab('PriLer')+ theme_classic()+
  geom_text_repel(label = df_snps$tissue, size = 3, alpha = 0.8, segment.size = 0.2, force = 10)+
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 0, hjust = 1),
        axis.text.y=element_text(size = 11), legend.position = 'none')+
  scale_color_manual(values = n_relGenes$col)

ggsave(filename = sprintf('%sFracRegSNPs_PriLerVSelnet_AllTissues.png', outFold), plot = pl_frac, width = 5, height = 5, dpi = 500)
ggsave(filename = sprintf('%sFracRegSNPs_PriLerVSelnet_AllTissues.pdf', outFold), plot = pl_frac, width = 5, height = 5, dpi = 500)

df_info <- cbind(df_nrel, df_snps[,!colnames(df_snps) %in% 'tissue'])
colnames(df_info)[3:4] <- c('nRelGenes_prior', 'nRelGenes_noprior')
write.table(file = sprintf('%scompare_PriLer_elnet_ngenes_nsnps.txt', outFold), df_info, quote = F, col.names = T, row.names = F, sep = '\t')


######################
## overall increase ##
######################

# tomodify

df_tab <- data.frame(method = rep(c('PriLer', other_name),2), n_genes = c(table(df_pl$corr2[df_pl$method == 'PriLer'] >= df_pl$corr2[df_pl$method == other_name]),
                                                                          table(df_pl$nsnps[df_pl$method == 'PriLer'] >= df_pl$nsnps[df_pl$method == other_name])),
                     type = c(rep('CV',2), rep('n.reg-SNPs',2)))
df_tab$type <- factor(df_tab$type)
group_name <- c(TeX('CV $corr^2$'),'n.reg-SNPs')
df_tab$method <- factor(df_tab$method, levels = c('PriLer', other_name))

pl_bar <- ggplot(df_tab, aes(x = type, y = n_genes, fill = method)) +
  geom_bar(stat="identity", alpha = 0.6, width = 0.7) + ggtitle('34 tissues combined')+
  xlab('') +  theme_classic() + ylab(TeX('n. genes')) + theme(legend.position = 'bottom', legend.title = element_blank(), legend.direction  = 'vertical', plot.title = element_text(hjust = 0.5), 
                                                              axis.text = element_text(size = 12), axis.title = element_text(size = 10), 
                                                              axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))+
  scale_fill_manual(labels = parse(text = c(TeX(paste0('PriLer $<$', other_name)), TeX(paste0('PriLer $\\geq$', other_name)))), values = rev(c("#1F77B4FF", color_other)))+
  scale_x_discrete(labels=parse(text = group_name))+ coord_flip()

file_name <- sprintf('/mnt/lucia/PriLer/comparison_methods_increaseAllTissues_PriLer_%s.png', other_name_save)
ggsave(file_name, pl_bar, width = 5, height = 3.5, units = "in", dpi=500)
file_name <- sprintf('/mnt/lucia/PriLer/comparison_methods_increaseAllTissues_PriLer_%s.pdf', other_name_save)
ggsave(file_name, pl_bar, width = 5, height = 3.5, units = "in", dpi=500)

