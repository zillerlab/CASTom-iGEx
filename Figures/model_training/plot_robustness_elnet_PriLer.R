# plot robustness analysis

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jaccard))
suppressPackageStartupMessages(library(aricode))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
options(bitmapType = 'cairo', device = 'png')

setwd('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2/Whole_Blood/200kb/')
nrep <- 10

# load baseline results
load('noGWAS/resPrior_regCoeffSnps_allchr.RData')
load('noGWAS/resNoPrior_regCoeffSnps_allchr.RData')
nSNPs <- sum(sapply(beta_snps_p, nrow))

regSNPs <- data.frame(PriLer_n = rep(0, nSNPs), elnet_n = rep(0, nSNPs), PriLer = rep(F, nSNPs), elnet = rep(F, nSNPs))
regSNPs$PriLer_n <- unlist(lapply(beta_snps_p, function(x) rowSums(x != 0)))
regSNPs$elnet_n <- unlist(lapply(beta_snps_nop, function(x) rowSums(x != 0)))
regSNPs$PriLer[unlist(lapply(beta_snps_p, function(x) rowSums(x != 0))) > 0] <- T
regSNPs$elnet[unlist(lapply(beta_snps_nop, function(x) rowSums(x != 0))) > 0] <- T

regSNPs_rep <- list()
for(i in 1:nrep){
  
  print(i)
  
  load(sprintf('robustness_analysis/rep%i/noGWAS/resPrior_regCoeffSnps_allchr.RData', i))
  load(sprintf('robustness_analysis/rep%i/noGWAS/resNoPrior_regCoeffSnps_allchr.RData', i))
  
  regSNPs_rep[[i]] <-  data.frame(PriLer_n = rep(0, nSNPs), elnet_n = rep(0, nSNPs), PriLer = rep(F, nSNPs), elnet = rep(F, nSNPs))
  regSNPs_rep[[i]]$PriLer_n <- unlist(lapply(beta_snps_p, function(x) rowSums(x != 0)))
  regSNPs_rep[[i]]$elnet_n <- unlist(lapply(beta_snps_nop, function(x) rowSums(x != 0)))
  regSNPs_rep[[i]]$PriLer[unlist(lapply(beta_snps_p, function(x) rowSums(x != 0))) > 0] <- T
  regSNPs_rep[[i]]$elnet[unlist(lapply(beta_snps_nop, function(x) rowSums(x != 0))) > 0] <- T
  
}

# compare with baseline
df_base <- data.frame(cor = rep(NA, nrep*2), jaccard = rep(NA, nrep*2), type = c(rep('PriLer', nrep),rep('el-net', nrep)))
df_base$cor <- c(sapply(regSNPs_rep, function(x) cor(x$PriLer_n, regSNPs$PriLer_n)), sapply(regSNPs_rep, function(x) cor(x$elnet_n, regSNPs$elnet_n)))
df_base$jaccard <- c(sapply(regSNPs_rep, function(x) jaccard(x$PriLer, regSNPs$PriLer)), sapply(regSNPs_rep, function(x) jaccard(x$elnet, regSNPs$elnet)))
# df_base$nmi <- c(sapply(regSNPs_rep, function(x) NMI(x$PriLer, regSNPs$PriLer)), sapply(regSNPs_rep, function(x) NMI(x$elnet, regSNPs$elnet)))
df_base$type <- factor(df_base$type, levels = c('el-net', 'PriLer'))

pl <- ggpaired(df_base, x = "type", y = "cor", fill = "type", alpha = 0.7, palette = c('white', 'grey60'),
               line.color = "gray", line.size = 0.4, add.params = list(fill = "white"))+ 
  stat_compare_means(method.args = list(alternative = "less"), paired = TRUE, label = "p.format", label.y = max(df_base$cor) + 0.01, label.x.npc = 'center')
pl <- ggpar(pl, legend = "none", xlab = '', ylab = 'correlation base VS repetition\nn. of genes regulated')
file_name <- sprintf('robustness_analysis/corr_compare_base')
ggsave(pl, filename = paste0(file_name, '.png'), width = 3, height = 5, dpi = 320, device = 'png')
ggsave(pl, filename = paste0(file_name, '.pdf'), width = 3, height = 5,  device = 'pdf')

pl <- ggpaired(df_base, x = "type", y = "jaccard", fill = "type", alpha = 0.7, palette = c('white', 'grey60'),
               line.color = "gray", line.size = 0.4, add.params = list(fill = "white"))+ 
  stat_compare_means(method.args = list(alternative = "less"), paired = TRUE, label = "p.format", label.y = max(df_base$jaccard) + 0.005, label.x.npc = 'center')
pl <- ggpar(pl, legend = "none", xlab = '', ylab = 'jaccard index base VS repetition\nreg-SNPs')
file_name <- sprintf('robustness_analysis/jaccard_compare_base')
ggsave(pl, filename = paste0(file_name, '.png'), width = 3, height = 5, dpi = 320, device = 'png')
ggsave(pl, filename = paste0(file_name, '.pdf'), width = 3, height = 5,  device = 'pdf')

############
# compare between repetitions
mat_cor <- mat_jacc <- matrix(NA, nrow = nrep, ncol = nrep)
for(i in 1:nrep){
  print(i)
  mat_cor[i,] <- sapply(regSNPs_rep, function(x) cor(x$PriLer_n, regSNPs_rep[[i]]$PriLer_n))
  mat_jacc[i,] <- sapply(regSNPs_rep, function(x) jaccard(x$PriLer, regSNPs_rep[[i]]$PriLer))
}

df_comp <- data.frame(cor = mat_cor[lower.tri(mat_cor)], jaccard = mat_jacc[lower.tri(mat_jacc)], type = rep('PriLer', sum(lower.tri(mat_jacc))))

mat_cor <- mat_jacc <- matrix(NA, nrow = nrep, ncol = nrep)
for(i in 1:nrep){
  print(i)
  mat_cor[i,] <- sapply(regSNPs_rep, function(x) cor(x$elnet_n, regSNPs_rep[[i]]$elnet_n))
  mat_jacc[i,] <- sapply(regSNPs_rep, function(x) jaccard(x$elnet, regSNPs_rep[[i]]$elnet))
}

df_comp <- rbind(df_comp, data.frame(cor = mat_cor[lower.tri(mat_cor)], jaccard = mat_jacc[lower.tri(mat_jacc)], type = rep('el-net', sum(lower.tri(mat_jacc)))))
df_comp$type <- factor(df_comp$type, levels = c( 'el-net','PriLer'))

pl <- ggboxplot(df_comp, x = "type", y = "cor", fill = "type", alpha = 0.7, palette = c('white', 'grey60'),
                add = "jitter", add.params = list(size = 0.7, jitter = 0.2),  
                outlier.shape = NA)+ 
  stat_compare_means(method.args = list(alternative = "less"), paired = TRUE, label = "p.format", label.y = max(df_comp$cor) + 0.01, label.x.npc = 'center')
pl <- ggpar(pl, legend = "none", xlab = '', ylab = 'correlation repetitions\nn. of genes regulated')
file_name <- sprintf('robustness_analysis/corr_compare_rep')
ggsave(pl, filename = paste0(file_name, '.png'), width = 3, height = 5, dpi = 320, device = 'png')
ggsave(pl, filename = paste0(file_name, '.pdf'), width = 3, height = 5,  device = 'pdf')

pl <- ggboxplot(df_comp, x = "type", y = "jaccard", fill = "type", alpha = 0.7, palette = c('white', 'grey60'),
                add = "jitter",  add.params = list(size = 0.7, jitter = 0.2),  
                outlier.shape = NA)+ 
  stat_compare_means(method.args = list(alternative = "less"), paired = TRUE, label = "p.format", label.y = max(df_comp$jaccard) + 0.005, label.x.npc = 'center')
pl <- ggpar(pl, legend = "none", xlab = '', ylab = 'jaccard index repetitions\nreg-SNPs')
file_name <- sprintf('robustness_analysis/jaccard_compare_rep')
ggsave(pl, filename = paste0(file_name, '.png'), width = 3, height = 5, dpi = 320, device = 'png')
ggsave(pl, filename = paste0(file_name, '.pdf'), width = 3, height = 5,  device = 'pdf')

# save df
write.table(df_base, file = 'robustness_analysis/compare_base_repetitions.txt', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(df_comp, file = 'robustness_analysis/compare_repetitions.txt', col.names = T, row.names = F, sep = '\t', quote = F)


                      
                      
                      

