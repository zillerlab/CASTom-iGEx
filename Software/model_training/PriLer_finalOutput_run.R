#!/usr/bin/env Rscript

#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggExtra))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="Combain all results and produce quality plot checks")

parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--priorDat_file", type = "character", help = "prior matrix file, common path (ending chr specific)")
parser$add_argument("--priorInf", type="integer", nargs = '*' , default = 0, help = "index prior feature to be used, if 0 all column are used [default %(default)s]")
parser$add_argument("--functR", type="character", help = "Rscript with functions to be used, complete path")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps [default %(default)s]")
parser$add_argument("--covDat_file", type = "character", help = "covariance file, complete path")
parser$add_argument("--Dx", type = "logical", default = FALSE, help = "if true, Dx included as covariate [default %(default)s]")
parser$add_argument("--part1Res_fold", type = "character", help = "path folder reuslt part1")
parser$add_argument("--part2Res_fold", type = "character", help = "path folder reuslt part2")
parser$add_argument("--part3Res_fold", type = "character", help = "path folder reuslt part3")
parser$add_argument("--part4Res_fold", type = "character", help = "path folder reuslt part4")
parser$add_argument("--convert_par", type="double",  default = 0.25, help = "parameter to rescale all priors, needed to start E par search not too low [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
InfoFold <- args$InfoFold
priorDat_file <- args$priorDat_file
priorInf <- args$priorInf
functR <- args$functR
covDat_file <- args$covDat_file
Dx <- args$Dx
cis_thres <- args$cis_thres
part1Res_fold <- args$part1Res_fold
part2Res_fold <- args$part2Res_fold
part3Res_fold <- args$part3Res_fold
part4Res_fold <- args$part4Res_fold
convert_par <- args$convert_par
outFold <- args$outFold

################################################################
# cis_thres = 200000
# part1Res_fold = './'
# part2Res_fold = part3Res_fold = part4Res_fold =  './'
# outFold <- './'
# InfoFold <- '../'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_training/PriLer_functions.R'
# priorDat_file <- '../priorMatrix_'
# priorInf <- c(2,3)
# covDat_file <- '../../INPUT_DATA/Covariates/Covariates_PEERfact_PCs.txt'
################################################################

covDat <- read.table(covDat_file, header = T, sep = '\t', stringsAsFactors = F, check.names=F)
sampleAnn <- covDat[, colnames(covDat) %in% c('Individual_ID', 'genoSample_ID', 'RNASample_ID')]

if(!Dx){col_exclude <- c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx')}
covDat <- covDat[, !colnames(covDat) %in% col_exclude]

source(functR)
all_Chroms <- 1:22

# total prior and seed/nfold CV
err_test <- read.table(sprintf('%scv_test_Epar_allchr.txt', part2Res_fold), header = T)
E_tot <- sapply(colnames(err_test), function(x) as.numeric(strsplit(x, 'E_h.')[[1]][2]))
names(E_tot) <- NULL
if(which.min(colMeans(err_test))<ncol(err_test)){
  id_opt <- which.min(colMeans(err_test))
  E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
  E_par <- -1
  print(sprintf('optimal E parameter: %.2f at position %i', E_hat, id_opt))
  res_p_her_tot_allchr <- get(load(sprintf('%sresPrior_EOpt_HeritableGenes_allchr.RData', part3Res_fold)))
  res_p_her_allchr <- get(load(sprintf('%sresPrior_EOpt_NestedCV_HeritableGenes_allchr.RData',part2Res_fold)))
}else{
  id_opt <- which(abs(diff(colMeans(err_test))) < 0.5)[1] + 1
  E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
  E_par <- E_hat
  print(sprintf('minimum not reached in the interval, E parameter based on convergence: %.2f at position %i', E_par, id_opt))
  res_p_her_tot_allchr <- get(load(sprintf('%sresPrior_EFixed%.2f_HeritableGenes_allchr.RData', part3Res_fold, E_hat)))
  res_p_her_allchr <- get(load(sprintf('%sresPrior_EFixed%.2f_NestedCV_HeritableGenes_allchr.RData',part2Res_fold,E_hat)))
}

res_nop_her_tot_allchr <- get(load(sprintf('%sresNoPrior_HeritableGenes_allchr.RData',part3Res_fold)))

# combine all results together, prediction performance and final model
## res without prior
res_nop <- vector(mode = 'list', length = length(all_Chroms))
beta_snps_nop <- vector(mode = 'list', length = length(all_Chroms))
beta_cov_nop <- vector(mode = 'list', length = length(all_Chroms))

res_p <- vector(mode = 'list', length = length(all_Chroms))
beta_snps_p <- vector(mode = 'list', length = length(all_Chroms))
beta_cov_p <- vector(mode = 'list', length = length(all_Chroms))

for(i in all_Chroms){
  
  print(i)
  
  ################################
  #### whithout prior results #### 
  ################################
  
  res_nop_her <- get(load(sprintf('%sresNoPrior_NestedCV_HeritableGenes_chr%i.RData',part1Res_fold, i)))
  res_nop_nother <- get(load(sprintf('%sresNoPrior_NestedCV_NotHeritableGenes_chr%i.RData',part4Res_fold, i)))
  res_nop_nother_tot <- get(load(sprintf('%sresNoPrior_NotHeritableGenes_chr%i.RData',part4Res_fold, i)))
  
  id <- sapply(res_nop_her$geneAnn$ensembl_gene_id, function(x) which(res_nop_her_tot_allchr$geneAnn$ensembl_gene_id==x))
  res_nop_her_tot <- list(geneAnn = res_nop_her_tot_allchr$geneAnn[id,], tot = res_nop_her_tot_allchr$tot[id,], 
                          beta_snps = res_nop_her_tot_allchr$beta_snps[[i]], beta_cov = res_nop_her_tot_allchr$beta_cov[[i]], seed = res_nop_nother_tot$seed)
  
  res_nop[[i]] <- as.data.frame(rbind(res_nop_her$geneAnn, res_nop_nother$geneAnn))
  
  res_nop[[i]]$train_dev <- c(rowMeans(sapply(res_nop_her$train, function(x) x$dev)), rowMeans(sapply(res_nop_nother$train, function(x) x$dev)))
  res_nop[[i]]$train_dev_geno <- c(rowMeans(sapply(res_nop_her$train, function(x) x$dev_geno)), rowMeans(sapply(res_nop_nother$train, function(x) x$dev_geno)))
  res_nop[[i]]$train_dev_cov <- c(rowMeans(sapply(res_nop_her$train, function(x) x$dev_cov)), rowMeans(sapply(res_nop_nother$train, function(x) x$dev_cov)))
  res_nop[[i]]$train_dev_geno_cov <- c(rowMeans(sapply(res_nop_her$train, function(x) x$dev_geno_cov)), rowMeans(sapply(res_nop_nother$train, function(x) x$dev_geno_cov)))
  res_nop[[i]]$train_cor <- c(rowMeans(sapply(res_nop_her$train, function(x) x$cor)), rowMeans(sapply(res_nop_nother$train, function(x) x$cor)))
  res_nop[[i]]$train_cor_pvalue <- c(rowMeans(sapply(res_nop_her$train, function(x) x$cor_pval)), rowMeans(sapply(res_nop_nother$train, function(x) x$cor_pval)))
  res_nop[[i]]$train_cor_noadj <- c(rowMeans(sapply(res_nop_her$train, function(x) x$cor_noadj)), rowMeans(sapply(res_nop_nother$train, function(x) x$cor_noadj)))
  res_nop[[i]]$train_cor_noadj_pvalue <- c(rowMeans(sapply(res_nop_her$train, function(x) x$cor_noadj_pval)), rowMeans(sapply(res_nop_nother$train, function(x) x$cor_noadj_pval)))
  
  res_nop[[i]]$test_dev <- c(rowMeans(sapply(res_nop_her$test, function(x) x$dev)), rowMeans(sapply(res_nop_nother$test, function(x) x$dev)))
  res_nop[[i]]$test_dev_geno <- c(rowMeans(sapply(res_nop_her$test, function(x) x$dev_geno)), rowMeans(sapply(res_nop_nother$test, function(x) x$dev_geno)))
  res_nop[[i]]$test_dev_cov <- c(rowMeans(sapply(res_nop_her$test, function(x) x$dev_cov)), rowMeans(sapply(res_nop_nother$test, function(x) x$dev_cov)))
  res_nop[[i]]$test_dev_geno_cov <- c(rowMeans(sapply(res_nop_her$test, function(x) x$dev_geno_cov)), rowMeans(sapply(res_nop_nother$test, function(x) x$dev_geno_cov)))
  res_nop[[i]]$test_cor <- c(rowMeans(sapply(res_nop_her$test, function(x) x$cor)), rowMeans(sapply(res_nop_nother$test, function(x) x$cor)))
  res_nop[[i]]$test_cor_pvalue <- c(rowMeans(sapply(res_nop_her$test, function(x) x$cor_pval)), rowMeans(sapply(res_nop_nother$test, function(x) x$cor_pval)))
  res_nop[[i]]$test_cor_noadj <- c(rowMeans(sapply(res_nop_her$test, function(x) x$cor_noadj)), rowMeans(sapply(res_nop_nother$test, function(x) x$cor_noadj)))
  res_nop[[i]]$test_cor_noadj_pvalue <- c(rowMeans(sapply(res_nop_her$test, function(x) x$cor_noadj_pval)), rowMeans(sapply(res_nop_nother$test, function(x) x$cor_noadj_pval)))
  
  res_nop[[i]]$test_comb_cor <- c(res_nop_her$cor_comb_test$cor, res_nop_nother$cor_comb_test$cor)
  res_nop[[i]]$test_comb_cor_pvalue <- c(res_nop_her$cor_comb_test$cor_pval, res_nop_nother$cor_comb_test$cor_pval)
  res_nop[[i]]$test_comb_cor_noadj <- c(res_nop_her$cor_comb_noadj_test$cor, res_nop_nother$cor_comb_noadj_test$cor)
  res_nop[[i]]$test_comb_cor_noadj_pvalue <- c(res_nop_her$cor_comb_noadj_test$cor_pval, res_nop_nother$cor_comb_noadj_test$cor_pval)
  
  res_nop[[i]] <- cbind(res_nop[[i]], rbind(res_nop_her_tot$tot[, !colnames(res_nop_her_tot$tot) %in% c('lambda', 'alpha', 'dev_lmgeno')], 
                                            res_nop_nother_tot$tot[, !colnames(res_nop_nother_tot$tot) %in% c('lambda', 'alpha', 'dev_lmgeno')]))
  
  ord_genes <- order(res_nop[[i]]$start_position)
  res_nop[[i]] <- res_nop[[i]][order(res_nop[[i]]$start_position), ]
  
  beta_snps_nop[[i]] <- cbind(res_nop_her_tot$beta_snps, res_nop_nother_tot$beta_snps)
  beta_snps_nop[[i]] <-  beta_snps_nop[[i]][, ord_genes]
  beta_cov_nop[[i]] <- cbind(res_nop_her_tot$beta_cov, res_nop_nother_tot$beta_cov)
  beta_cov_nop[[i]] <-  beta_cov_nop[[i]][, ord_genes]
  
  
  ############################
  #### whit prior results #### 
  ############################
  
  res_p_nother <- get(load(sprintf('%sresPrior_NestedCV_NotHeritableGenes_chr%i.RData',part4Res_fold, i)))
  res_p_nother_tot <- get(load(sprintf('%sresPrior_NotHeritableGenes_chr%i.RData', part4Res_fold, i)))
  
  id_chr <- which(res_p_her_allchr$geneAnn$chrom == paste0('chr', i))
  print(identical(res_p_her_allchr$geneAnn[id_chr,], res_p_her_tot_allchr$geneAnn[id_chr,]))
  
  res_p_her_tot <- list(geneAnn = res_p_her_tot_allchr$geneAnn[id_chr,], tot = res_p_her_tot_allchr$tot[id_chr,], 
                        beta_snps = res_p_her_tot_allchr$beta_snps[[i]], beta_cov = res_p_her_tot_allchr$beta_cov[[i]], seed = res_p_her_tot_allchr$seed)
  colnames(res_p_her_tot$tot)[which(colnames(res_p_her_tot$tot) == 'dev_genocov')] = 'dev_geno_cov'
  
  res_p_her <- list(geneAnn = res_p_her_allchr$geneAnn[id_chr,], train = lapply(res_p_her_allchr$train_opt, function(x) x[id_chr,]), 
                    test = lapply(res_p_her_allchr$test_opt, function(x) x[id_chr,]), 
                    cor_comb_test = res_p_her_allchr$cor_comb_test_opt[id_chr,], 
                    cor_comb_noadj_test = res_p_her_allchr$cor_comb_test_noadj_opt[id_chr,])
  
  res_p[[i]] <- as.data.frame(rbind(res_p_her_allchr$geneAnn[id_chr,], res_p_nother$geneAnn))
  
  res_p[[i]]$train_dev <- c(rowMeans(sapply(res_p_her$train, function(x) x$dev)), rowMeans(sapply(res_p_nother$train, function(x) x$dev)))
  res_p[[i]]$train_dev_geno <- c(rowMeans(sapply(res_p_her$train, function(x) x$dev_geno)), rowMeans(sapply(res_p_nother$train, function(x) x$dev_geno)))
  res_p[[i]]$train_dev_cov <- c(rowMeans(sapply(res_p_her$train, function(x) x$dev_cov)), rowMeans(sapply(res_p_nother$train, function(x) x$dev_cov)))
  res_p[[i]]$train_dev_geno_cov <- c(rowMeans(sapply(res_p_her$train, function(x) x$dev_geno_cov)), rowMeans(sapply(res_p_nother$train, function(x) x$dev_geno_cov)))
  res_p[[i]]$train_cor <- c(rowMeans(sapply(res_p_her$train, function(x) x$cor)), rowMeans(sapply(res_p_nother$train, function(x) x$cor)))
  res_p[[i]]$train_cor_pvalue <- c(rowMeans(sapply(res_p_her$train, function(x) x$cor_pval)), rowMeans(sapply(res_p_nother$train, function(x) x$cor_pval)))
  res_p[[i]]$train_cor_noadj <- c(rowMeans(sapply(res_p_her$train, function(x) x$cor_noadj)), rowMeans(sapply(res_p_nother$train, function(x) x$cor_noadj)))
  res_p[[i]]$train_cor_noadj_pvalue <- c(rowMeans(sapply(res_p_her$train, function(x) x$cor_noadj_pval)), rowMeans(sapply(res_p_nother$train, function(x) x$cor_noadj_pval)))
  
  
  res_p[[i]]$test_dev <- c(rowMeans(sapply(res_p_her$test, function(x) x$dev)), rowMeans(sapply(res_p_nother$test, function(x) x$dev)))
  res_p[[i]]$test_dev_geno <- c(rowMeans(sapply(res_p_her$test, function(x) x$dev_geno)), rowMeans(sapply(res_p_nother$test, function(x) x$dev_geno)))
  res_p[[i]]$test_dev_cov <- c(rowMeans(sapply(res_p_her$test, function(x) x$dev_cov)), rowMeans(sapply(res_p_nother$test, function(x) x$dev_cov)))
  res_p[[i]]$test_dev_geno_cov <- c(rowMeans(sapply(res_p_her$test, function(x) x$dev_geno_cov)), rowMeans(sapply(res_p_nother$test, function(x) x$dev_geno_cov)))
  res_p[[i]]$test_cor <- c(rowMeans(sapply(res_p_her$test, function(x) x$cor)), rowMeans(sapply(res_p_nother$test, function(x) x$cor)))
  res_p[[i]]$test_cor_pvalue <- c(rowMeans(sapply(res_p_her$test, function(x) x$cor_pval)), rowMeans(sapply(res_p_nother$test, function(x) x$cor_pval)))
  res_p[[i]]$test_cor_noadj <- c(rowMeans(sapply(res_p_her$test, function(x) x$cor_noadj)), rowMeans(sapply(res_p_nother$test, function(x) x$cor_noadj)))
  res_p[[i]]$test_cor_noadj_pvalue <- c(rowMeans(sapply(res_p_her$test, function(x) x$cor_noadj_pval)), rowMeans(sapply(res_p_nother$test, function(x) x$cor_noadj_pval)))
  
  
  res_p[[i]]$test_comb_cor <- c(res_p_her$cor_comb_test$cor, res_p_nother$cor_comb_test$cor)
  res_p[[i]]$test_comb_cor_pvalue <- c(res_p_her$cor_comb_test$cor_pval, res_p_nother$cor_comb_test$cor_pval)
  res_p[[i]]$test_comb_cor_noadj <- c(res_p_her$cor_comb_noadj_test$cor, res_p_nother$cor_comb_noadj_test$cor)
  res_p[[i]]$test_comb_cor_noadj_pvalue <- c(res_p_her$cor_comb_noadj_test$cor_pval, res_p_nother$cor_comb_noadj_test$cor_pval)
  
  res_p[[i]] <- cbind(res_p[[i]], rbind(res_p_her_tot$tot[, !colnames(res_p_her_tot$tot) %in% c('lambda', 'alpha', 'dev_lmgeno')], 
                                        res_p_nother_tot$tot[, !colnames(res_p_nother_tot$tot) %in% c('lambda', 'alpha', 'dev_lmgeno')]))
  
  
  ord_genes <- order(res_p[[i]]$start_position)
  res_p[[i]] <- res_p[[i]][order(res_p[[i]]$start_position), ]
  
  beta_snps_p[[i]] <- cbind(res_p_her_tot$beta_snps, res_p_nother_tot$beta_snps)
  beta_snps_p[[i]] <-  beta_snps_p[[i]][, ord_genes]
  beta_cov_p[[i]] <- cbind(res_p_her_tot$beta_cov, res_p_nother_tot$beta_cov)
  beta_cov_p[[i]] <-  beta_cov_p[[i]][, ord_genes]
  
  
}

res_nop <- do.call(rbind, res_nop)
beta_cov_nop <- do.call(cbind, beta_cov_nop)
beta_cov_nop <- as.data.frame(as.matrix(t(beta_cov_nop)))
colnames(beta_cov_nop) <- c(colnames(covDat), 'intercept')

res_p <- do.call(rbind, res_p)
beta_cov_p <- do.call(cbind, beta_cov_p)
beta_cov_p <- as.data.frame(as.matrix(t(beta_cov_p)))
colnames(beta_cov_p) <- c(colnames(covDat), 'intercept')


# save results
write.table(file = sprintf('%sresNoPrior_regEval_allchr.txt', outFold), x = res_nop, col.names = T, row.names = F, sep = '\t', quote = F)
write.table(file = sprintf('%sresPrior_regEval_allchr.txt', outFold), x = res_p, col.names = T, row.names = F, sep = '\t', quote = F)

write.table(file = sprintf('%sresNoPrior_regCoeffCov_allchr.txt', outFold), x = beta_cov_nop, col.names = T, row.names = F, sep = '\t', quote = F)
write.table(file = sprintf('%sresPrior_regCoeffCov_allchr.txt', outFold), x = beta_cov_p, col.names = T, row.names = F, sep = '\t', quote = F)

save(beta_snps_nop, file = sprintf('%sresNoPrior_regCoeffSnps_allchr.RData', outFold))
save(beta_snps_p, file = sprintf('%sresPrior_regCoeffSnps_allchr.RData', outFold))

##########################################################################################################################################################################
##### plots 

###################
#### CV error #####
###################

err_train <- read.table(sprintf('%scv_train_Epar_allchr.txt', part2Res_fold), header = T)
df_train <- data.frame(mean_err = colMeans(err_train), E_par = E_tot, sd_err = apply(err_train,2,sd))
df_test <- data.frame(mean_err = colMeans(err_test), E_par = E_tot, sd_err = apply(err_test,2,sd), opt = rep('not', length(E_tot)), stringsAsFactors = F)

if(E_par > 0){
  df_test$opt[id_opt] <- 'conv'
  df_test$opt[which.min(df_test$mean_err)] <- 'opt'
}else{
  df_test$opt[id_opt] <- 'opt'
}
df_test$opt <- factor(df_test$opt, level = c('not', 'opt', 'conv'))
df_train$opt <- df_test$opt

df_tot <- rbind(df_train, df_test)  
df_tot$type <- c(rep('Train Set', nrow(df_train)),  rep('Test Set', nrow(df_test)))
df_tot$type <- factor(df_tot$type, levels = c('Train Set', 'Test Set'))

plot_errbars<- ggplot(data = df_tot, aes(x = E_par, y = mean_err)) +
  geom_line(position=position_dodge(0.1)) +theme_bw()+
  geom_point(aes(color = opt), position=position_dodge(0.1), size = 2)+
  geom_errorbar(aes(ymin=mean_err-sd_err, ymax=mean_err+sd_err), width=.02,
                position=position_dodge(0.1))+
  facet_wrap(.~type, scales = 'free_y', nrow = 1)+
  scale_color_manual(values = c('black', 'red', 'blue'))+
  xlab('E par')+ylab('CV final MSE')+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16),legend.position = 'bottom',
        axis.text.x=element_text(size = 11),
        axis.text.y=element_text(size = 11))+
  scale_x_continuous("E par", labels = as.character(E_tot), breaks = E_tot)
ggsave(sprintf("%splot_CVerr_comp.png", outFold), plot = plot_errbars, width = 10, height = 6, units = "in")
ggsave(sprintf("%splot_CVerr_comp.pdf", outFold), plot = plot_errbars, width = 10, height = 6, units = "in")

#######################
#### CV iteration ##### 
#######################

res_cv_test <- get(load(sprintf('%sevalf_Epar_cvtest_allchr.RData', part2Res_fold)))
res_cv_train <- get(load(sprintf('%sevalf_Epar_cvtrain_allchr.RData', part2Res_fold)))
id_Eopt <- id_opt

res_cv_test <- res_cv_test[[id_Eopt]]
res_cv_train <- res_cv_train[[id_Eopt]]
nfold_out <- length(res_cv_test)

df_cv_test <- data.frame(fold = unlist(lapply(1:nfold_out, function(x) rep(sprintf('fold_%i', x), length(res_cv_test[[x]])))))
df_cv_test$fold <- factor(df_cv_test$fold)
df_cv_test$error <- unlist(res_cv_test)
df_cv_test$it <- unlist(lapply(1:nfold_out, function(x) 1:length(res_cv_test[[x]])))

df_cv_train <- data.frame(fold = unlist(lapply(1:nfold_out, function(x) rep(sprintf('fold_%i', x), length(res_cv_train[[x]])))))
df_cv_train$fold <- factor(df_cv_test$fold)
df_cv_train$error <- unlist(res_cv_train)
df_cv_train$it <- unlist(lapply(1:nfold_out, function(x) 1:length(res_cv_train[[x]])))


plot_cverr_test<- ggplot(data = df_cv_test, aes(x = it, y = error)) +
  geom_line() +theme_bw()+
  geom_point()+
  facet_wrap(.~fold, ncol=1, scale = 'free')+
  xlab('iteration')+ylab('CV iteration MSE')+ ggtitle('Test set')+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16),
        axis.text.x=element_text(size = 11),
        axis.text.y=element_text(size = 11))

if(E_par>0){
  file_name <- sprintf("%splot_CVmse_TEST_Eopt", outFold)
}else{
  file_name <- sprintf("%splot_CVmse_TEST_Efixed%s", outFold, as.character(E_hat))
}
ggsave(sprintf('%s.pdf',file_name), plot = plot_cverr_test, width = 8, height = 8, units = "in")
ggsave(sprintf('%s.png',file_name), plot = plot_cverr_test, width = 8, height = 8, units = "in")


plot_cverr_train <- ggplot(data = df_cv_train, aes(x = it, y = error)) +
  geom_line() +theme_bw()+
  geom_point()+
  facet_wrap(.~fold, ncol=1, scale = 'free')+
  xlab('iteration')+ylab('CV iteration MSE')+ ggtitle('Train set')+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16),
        axis.text.x=element_text(size = 11),
        axis.text.y=element_text(size = 11))

if(E_par>0){
  file_name <- sprintf("%splot_CVmse_TRAIN_Eopt", outFold)
}else{
  file_name <- sprintf("%splot_CVmse_TRAIN_Efixed%s", outFold, as.character(E_hat))
}
ggsave(sprintf('%s.pdf',file_name), plot = plot_cverr_train, width = 8, height = 8, units = "in")
ggsave(sprintf('%s.png',file_name), plot = plot_cverr_train, width = 8, height = 8, units = "in")


# plot objective function on the train and test set
res_objcv_train <- get(load(sprintf('%sobj_Epar_cvtrain_allchr.RData', outFold)))
res_objcv_train <- res_objcv_train[[id_Eopt]]

df_objcv_train <- data.frame(fold = unlist(lapply(1:nfold_out, function(x) rep(sprintf('fold_%i', x), length(res_objcv_train[[x]])))))
df_objcv_train$fold <- factor(df_objcv_train$fold)
df_objcv_train$error <- unlist(res_objcv_train)
df_objcv_train$it <- unlist(lapply(1:nfold_out, function(x) 1:length(res_objcv_train[[x]])))

plot_cvobjerr_train<- ggplot(data = df_objcv_train, aes(x = it, y = error)) +
  geom_line() +theme_bw()+
  geom_point()+
  facet_wrap(.~fold, ncol=1, scale = 'free')+
  xlab('iteration')+ylab('CV iteration obj')+ ggtitle('Train set')+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16),
        axis.text.x=element_text(size = 11),
        axis.text.y=element_text(size = 11))

if(E_par>0){
  file_name <- sprintf("%splot_CVobj_TRAIN_Eopt", outFold)
}else{
  file_name <- sprintf("%splot_CVobj_TRAIN_Efixed%s", outFold, as.character(E_hat))
}
ggsave(sprintf('%s.pdf',file_name), plot = plot_cvobjerr_train, width = 8, height = 8, units = "in")
ggsave(sprintf('%s.png',file_name), plot = plot_cvobjerr_train, width = 8, height = 8, units = "in")

res_objcv_test <- get(load(sprintf('%sobj_Epar_cvtest_allchr.RData', outFold)))
res_objcv_test <- res_objcv_test[[id_Eopt]]

df_objcv_test <- data.frame(fold = unlist(lapply(1:nfold_out, function(x) rep(sprintf('fold_%i', x), length(res_objcv_test[[x]])))))
df_objcv_test$fold <- factor(df_objcv_test$fold)
df_objcv_test$error <- unlist(res_objcv_test)
df_objcv_test$it <- unlist(lapply(1:nfold_out, function(x) 1:length(res_objcv_test[[x]])))

plot_cvobjerr_test<- ggplot(data = df_objcv_test, aes(x = it, y = error)) +
  geom_line() +theme_bw()+
  geom_point()+
  facet_wrap(.~fold, ncol=1, scale = 'free')+
  xlab('iteration')+ylab('CV iteration obj')+ ggtitle('Test set')+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16),
        axis.text.x=element_text(size = 11),
        axis.text.y=element_text(size = 11))

if(E_par>0){
  file_name <- sprintf("%splot_CVobj_TEST_Eopt", outFold)
}else{
  file_name <- sprintf("%splot_CVobj_TEST_Efixed%s", outFold, as.character(E_hat))
}
ggsave(sprintf('%s.pdf',file_name), plot = plot_cvobjerr_test, width = 8, height = 8, units = "in")
ggsave(sprintf('%s.png',file_name), plot = plot_cvobjerr_test, width = 8, height = 8, units = "in")

#############################
#### CV dev (all genes) #####
#############################

dev_train_test <- data.frame(train = res_p$train_dev, test = res_p$test_dev, type = res_p$type) 
id_na <-  is.na(rowSums(dev_train_test[, 1:2])) 
dev_train_test <- dev_train_test[!id_na, ]
id_zero <- dev_train_test$train == 0 # correspond to NA reuslts (no snps)
dev_train_test <- dev_train_test[!id_zero, ]
dev_train_test$type <- factor(dev_train_test$type)

pl_devcv <- ggplot(dev_train_test, aes(x = train, y = test, color = type, group = type)) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed")+ylim(min(dev_train_test$test), 1)+
  xlim(min(dev_train_test$train), 1)+
  geom_point(size = 0.5) + ggtitle(bquote('Mean CV R'^2))+
  xlab('Train Set') +  ylab('Test Set')+ theme_bw() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = c('black', 'darkgrey'))

file_name <- sprintf('%splot_CVdev_trainVStest_allgenes', outFold)
ggsave(sprintf('%s.png',file_name), ggMarginal(pl_devcv, groupColour = TRUE, groupFill = TRUE), width = 5, height = 5, units = "in")
ggsave(sprintf('%s.pdf',file_name), ggMarginal(pl_devcv, groupColour = TRUE, groupFill = TRUE), width = 5, height = 5, units = "in")

##################################
#### CV dev geno (all genes) #####
##################################

devgeno_train_test <- data.frame(train = res_p$train_dev_geno, test = res_p$test_dev_geno, type = res_p$type) 
id_na <-  is.na(rowSums(devgeno_train_test[, 1:2])) 
devgeno_train_test <- devgeno_train_test[!id_na, ]
id_zero <- devgeno_train_test$train == 0 # correspond to NA reuslts (no snps)
devgeno_train_test <- devgeno_train_test[!id_zero, ]
devgeno_train_test$type <- factor(devgeno_train_test$type)

pl_devgenocv <- ggplot(devgeno_train_test, aes(x = train, y = test, color = type, group = type)) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed")+ylim(min(devgeno_train_test$test), 1)+
  xlim(min(devgeno_train_test$train), 1)+
  geom_point(size = 0.5) + ggtitle(bquote('Mean CV genotype R'^2))+
  xlab('Train Set') +  ylab('Test Set')+ theme_bw() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = c('black', 'darkgrey'))

file_name <- sprintf('%splot_CVdevgeno_trainVStest_allgenes', outFold)
ggsave(sprintf('%s.png',file_name), ggMarginal(pl_devgenocv, groupColour = TRUE, groupFill = TRUE), width = 5, height = 5, units = "in")
ggsave(sprintf('%s.pdf',file_name), ggMarginal(pl_devgenocv, groupColour = TRUE, groupFill = TRUE), width = 5, height = 5, units = "in")

##########################################################
#### CV dev geno (all genes) diff with/without prior #####
##########################################################

# put together both with and without ptiot, consider only reliable genes
tot_df <- res_p[, colnames(res_p) %in% c('type', 'ensembl_gene_id')]
tot_df$train_dev_geno_nop <- res_nop$train_dev_geno
tot_df$test_dev_geno_nop <- res_nop$test_dev_geno
tot_df$train_dev_geno_p <- res_p$train_dev_geno
tot_df$test_dev_geno_p <- res_p$test_dev_geno

id_filt <- which(res_p$dev_geno >= 0.01 & res_p$test_dev_geno >0)
tot_df <- tot_df[id_filt,]

tot_df$perc_train <- (tot_df$train_dev_geno_p-tot_df$train_dev_geno_nop)/tot_df$train_dev_geno_nop
tot_df$perc_test <- (tot_df$test_dev_geno_p-tot_df$test_dev_geno_nop)/tot_df$test_dev_geno_nop
tot_df$train_col <- 'red'
tot_df$test_col <- 'red'

tot_df$train_col[tot_df$train_dev_geno_p<tot_df$train_dev_geno_nop] <- 'blue'
tot_df$test_col[tot_df$test_dev_geno_p<tot_df$test_dev_geno_nop] <- 'blue'
tot_df$train_col <- factor(tot_df$train_col, levels = c('red', 'blue'))
tot_df$test_col <- factor(tot_df$test_col, levels = c('red', 'blue')) 

pos_text <- (max(tot_df$train_dev_geno_nop[tot_df$train_dev_geno_nop>=0.001])-min(tot_df$train_dev_geno_nop[tot_df$train_dev_geno_nop>=0.001]))/2

pl_genonew_train <- ggplot(subset(tot_df, train_dev_geno_nop>=0.001), aes(x = train_dev_geno_nop, y = perc_train, color = train_col)) + 
  geom_point(size = 0.1) + ggtitle('Train sets (reliable genes)')+ 
  geom_hline(yintercept =0, color = 'black', size = 0.5)+
  geom_vline(xintercept = 0.0009, color = 'black', size = 0.5, linetype =2, alpha = 0.6)+
  xlab(expression (R[el-net]^2)) +  xlab(expression (R[el-net]^2))+
  theme_classic() + ylab(expression(frac(R[el-net-prior]^2 - R[el-net]^2,R[el-net]^2))) + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x=pos_text, y=max(tot_df$perc_train[tot_df$train_dev_geno_nop>=0.001])-0.1, 
           label=paste("% improved genes\n", round(length(which(tot_df$train_col == 'red'))/nrow(tot_df),digits=3)), color="black", size = 3)+
  scale_color_manual(values =c("red", "blue"))

file_name <- sprintf('%splot_devtrain_compareWithPrior_increase', outFold)
ggsave(sprintf('%s.png',file_name), pl_genonew_train, width = 4, height = 4, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name), pl_genonew_train, width = 4, height = 4, units = "in", dpi = 500)

pos_text <- (max(tot_df$test_dev_geno_nop[tot_df$test_dev_geno_nop>=0.001])-min(tot_df$test_dev_geno_nop[tot_df$test_dev_geno_nop>=0.001]))/2

pl_genonew_test <- ggplot(subset(tot_df, test_dev_geno_nop>=0.001), aes(x = test_dev_geno_nop, y = perc_test, color = test_col)) + 
  geom_point(size = 0.1) + ggtitle('Test sets (reliable genes)')+ 
  geom_hline(yintercept =0, color = 'black', size = 0.5)+
  geom_vline(xintercept = 0.0009, color = 'black', size = 0.5, linetype =2, alpha = 0.6)+
  xlab(expression (R[el-net]^2)) +  xlab(expression (R[el-net]^2))+
  theme_classic() + ylab(expression(frac(R[el-net-prior]^2 - R[el-net]^2,R[el-net]^2))) + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x=pos_text, y=max(tot_df$perc_test[tot_df$test_dev_geno_nop>=0.001])-0.1,
           label=paste("% improved genes\n", round(length(which(tot_df$test_col == 'red'))/nrow(tot_df),digits=3)), color="black", size = 3)+
  scale_color_manual(values =c("red", "blue"))

file_name <- sprintf('%splot_devtest_compareWithPrior_increase', outFold)
ggsave(sprintf('%s.png',file_name), pl_genonew_test, width = 4, height = 4, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name), pl_genonew_test, width = 4, height = 4, units = "in", dpi = 500)

##################################################################
#### Final dev decomposition (all genes)  with/without prior #####
##################################################################

df_tot <- data.frame(dev_type = rep(c(rep('total', nrow(res_p)), rep('geno', nrow(res_p)), rep('cov', nrow(res_p)), rep('cov_geno', nrow(res_p))),2), 
                     reg_type = c(rep('el-net-prior', nrow(res_p)*4), rep('el-net', nrow(res_p)*4)), 
                     dev = c(res_p$dev, res_p$dev_geno, res_p$dev_cov, res_p$dev_geno_cov, res_nop$dev, res_nop$dev_geno, res_nop$dev_cov, res_nop$dev_geno_cov))

df_tot$reg_type <- factor(df_tot$reg_type, levels = c('el-net', 'el-net-prior'))

# boxplot
pl_box <- ggplot(df_tot, aes(x = dev_type, y = dev, fill = reg_type)) + 
  geom_boxplot(color = 'black') + 
  xlab('') +  theme_bw() + ylab(expression(R^2)) + theme(legend.position = "bottom")

file_name <- sprintf('%sboxplot_devainceDecomposition_elVSelP_allgenes', outFold)
ggsave(plot = pl_box, file = sprintf('%s.png',file_name), width = 5, height = 6)
ggsave(plot = pl_box, file = sprintf('%s.pdf',file_name), width = 5, height = 6)

###########################################################
#### Final dev sorted (all genes)  with/without prior #####
###########################################################

tmp <- cbind(res_p$dev, res_nop$dev)
id_rm <- is.na(rowSums(tmp))
ngenes <- ifelse(length(which(id_rm))>0, nrow(res_p) -length(which(id_rm)), nrow(res_p) )

df_tot <- data.frame(type = c(rep('el-net-prior', ngenes), rep('el-net',ngenes)), 
                     deviance = c(sort(res_p$dev[!id_rm], decreasing = T), sort(res_nop$dev[!id_rm], decreasing = T)))

df_tot$type <- factor(df_tot$type,  levels = c('el-net', 'el-net-prior'))
df_tot$gene <- c(rep(1:ngenes, 2))

pl_dev <- ggplot(df_tot, aes(x = gene, y = deviance, color = type)) + 
  geom_line() + ggtitle('Evaluation entire set')+ ylab(expression(R^2))+
  xlab(bquote('Genes sorted by R'^2)) +  theme_bw() +
  theme(legend.position = "bottom",axis.text.x=element_text(size = 11), axis.text.y=element_text(size = 11))+
  scale_color_brewer(palette="Set1")

file_name <- sprintf('%splot_compare_dev_elVSelP_allgenes', outFold)
ggsave(sprintf('%s.png',file_name), pl_dev, width = 4.5, height = 4.5, units = "in")
ggsave(sprintf('%s.pdf',file_name), pl_dev, width = 4.5, height = 4.5, units = "in")

# geno
tmp <- cbind(res_p$dev_geno, res_nop$dev_geno)
id_rm <- is.na(rowSums(tmp))
ngenes <- ifelse(length(which(id_rm))>0, nrow(res_p) -length(which(id_rm)), nrow(res_p) )

df_tot <- data.frame(type = c(rep('el-net-prior', ngenes), rep('el-net',ngenes)), 
                     deviance = c(sort(res_p$dev_geno[!id_rm], decreasing = T), sort(res_nop$dev_geno[!id_rm], decreasing = T)))

df_tot$type <- factor(df_tot$type,  levels = c('el-net', 'el-net-prior'))
df_tot$gene <- c(rep(1:ngenes, 2))

pl_dev <- ggplot(df_tot, aes(x = gene, y = deviance, color = type)) + 
  geom_line() + ggtitle('Evaluation entire set')+ ylab(expression(genotype~R^2))+
  xlab(bquote('Genes sorted by R'^2)) +  theme_bw() +
  theme(legend.position = "bottom",axis.text.x=element_text(size = 11), axis.text.y=element_text(size = 11))+
  scale_color_brewer(palette="Set1")

file_name <- sprintf('%splot_compare_devgeno_elVSelP_allgenes', outFold)
ggsave(sprintf('%s.png',file_name), pl_dev, width = 4.5, height = 4.5, units = "in")
ggsave(sprintf('%s.pdf',file_name), pl_dev, width = 4.5, height = 4.5, units = "in")


#############################################################
#### number of reg-SNPs (all genes)  with/without prior #####
#############################################################

id_snps_gene_nop <- lapply(beta_snps_nop, function(x) rowSums(x)!=0)
id_snps_gene_p <- lapply(beta_snps_p, function(x) rowSums(x)!=0)

pDat <- list()
snpPos <- list()
for(i in 1:length(all_Chroms)){
  
  # genotype
  chr <- paste0('chr',all_Chroms[i])
  print(chr)
  
  ## prior matrix
  pDat[[i]] <- read.table(gzfile(sprintf('%s%s.txt.gz', priorDat_file, chr)), header = T, sep = '\t')
  pNames <- colnames(pDat[[i]])[priorInf]
  pDat[[i]] <- as.matrix(pDat[[i]][,priorInf]) # only cell type of interest
  pDat[[i]] <- pDat[[i]]*convert_par
  snpPos[[i]] <- read.table(sprintf('%shg19_SNPs_%s_matched.txt', InfoFold, chr), header = T, sep = '\t')
  
  
}

tmp <- do.call(rbind, pDat)
# pNames <- colnames(tmp)
print(pNames)

# N_k <- apply(tmp, 2, function(x) mean(x))
# pDat <- lapply(pDat, function(x)  sweep(x, 2, N_k/min(N_k), '/'))

id_snps_prior <- sapply(pDat, function(x) rowSums(x)>0)

id_snps_gene_and_prior_p <- mapply(function(x,y) x & y, x = id_snps_gene_p, y = id_snps_prior)
id_snps_gene_and_prior_nop <- mapply(function(x,y) x & y, x = id_snps_gene_nop, y = id_snps_prior)

################
tmp_nop <- data.frame(id = which(unlist(id_snps_gene_nop)), prior = unlist(id_snps_gene_and_prior_nop)[which(unlist(id_snps_gene_nop))])
tmp_p <- data.frame(id = which(unlist(id_snps_gene_p)), prior = unlist(id_snps_gene_and_prior_p)[which(unlist(id_snps_gene_p))])

df_len <- data.frame(type = c(rep('el-net',4),rep('el-net learned prior',4)), type_snp = rep(c('without prior (common)', 'without prior (unique)', 'with prior (common)', 'with prior (unique)'),2))

val = length(intersect(tmp_nop$id[!tmp_nop$prior], tmp_p$id[!tmp_p$prior]))
val = c(val,length(setdiff(tmp_nop$id[!tmp_nop$prior], intersect(tmp_nop$id[!tmp_nop$prior], tmp_p$id[!tmp_p$prior]))))
val = c(val, length(intersect(tmp_nop$id[tmp_nop$prior], tmp_p$id[tmp_p$prior])))
val = c(val,length(setdiff(tmp_nop$id[tmp_nop$prior], intersect(tmp_nop$id[tmp_nop$prior], tmp_p$id[tmp_p$prior]))))

val = c(val, length(intersect(tmp_nop$id[!tmp_nop$prior], tmp_p$id[!tmp_p$prior])))
val = c(val,length(setdiff(tmp_p$id[!tmp_p$prior], intersect(tmp_nop$id[!tmp_nop$prior], tmp_p$id[!tmp_p$prior]))))
val = c(val, length(intersect(tmp_nop$id[tmp_nop$prior], tmp_p$id[tmp_p$prior])))
val = c(val,length(setdiff(tmp_p$id[tmp_p$prior], intersect(tmp_nop$id[tmp_nop$prior], tmp_p$id[tmp_p$prior]))))
df_len$len <- val
df_len$type_snp <- factor(df_len$type_snp, levels = (c('with prior (unique)', 'with prior (common)', 'without prior (unique)', 'without prior (common)')))

pl_hist <- ggplot(df_len, aes(x = type, y = val, fill = type_snp)) + 
  geom_bar(stat="identity", color = 'black', alpha = 0.6) + 
  xlab('') +  theme_classic() + ylab('n. reg-SNPs') + theme(legend.position = "right")+
  labs(fill = "")+
  scale_fill_manual(values =c("#1F78B4", "#A6CEE3", "#E31A1C", "#FB9A99"))

file_name <- sprintf('%splot_numberSNPs_elVSelP_gene', outFold)
ggsave(sprintf('%s.png',file_name), pl_hist, width = 5.1, height = 5.1, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name), pl_hist, width = 5.1, height = 5.1, units = "in", dpi = 500)
write.table(file = sprintf('%snVariants_prior_compare.txt', outFold),x = df_len, col.names = T, row.names = F, sep ='\t', quote = F)

###############################################
#### prior coeff distribution across snps #####
###############################################

coeffPrior <- res_p_her_tot_allchr$prior_coeff
df_pr <- data.frame(prior_coeff = unlist(coeffPrior), reg_snps = unlist(id_snps_gene_p))

# correct position
ends <- sapply(snpPos, function(x) max(x$position))/1000
starts <- sapply(snpPos, function(x) min(x$position))/1000
corr_pos <- c(0,cumsum(ends)[-length(all_Chroms)])
new_pos <- mapply(function(x,y) x$position/1000+y, x = snpPos, y = corr_pos)
df_pr$pos = unlist(new_pos)

df_pr$type <- 'reg-SNPs'
df_pr$type[!df_pr$reg_snps] <- 'not reg-SNPs'
df_pr$type <- factor(df_pr$type)

totdensity_p <- ggplot(df_pr) + geom_density(mapping = aes(x = prior_coeff, fill = type, group = type), position = "stack", alpha = 0.5)+ 
  ylab('density') +  xlab('Prior coefficient')+ theme_bw() + theme(legend.position = "bottom")

# pl_priordist <- ggplot(df_pr[1:1000000,], aes(x = pos, y = prior_coeff, color = type, group = type)) + 
#   geom_point(size = 0.01) + 
#   xlab('') +  ylab('Prior coefficient')+ theme_bw() + theme(legend.position = "bottom")+
#   scale_colour_manual(values = c('black', 'darkgrey'))

file_name <- sprintf('%splot_priorCoeff_dist', outFold)
ggsave(sprintf('%s.png',file_name), totdensity_p, width = 7, height = 5, units = "in")
ggsave(sprintf('%s.pdf',file_name), totdensity_p, width = 7, height = 5, units = "in")

#######################################################################
#### contribution of prior features to prior coeff, only reg-SNPs #####
#######################################################################

weights <- res_p_her_tot_allchr$weights
contr_prior <- lapply(pDat, function(x) t(t(x)*weights))

df_feat <- do.call(rbind, contr_prior)
df_feat <- as.matrix(df_feat[unlist(id_snps_gene_p),], ncol = length(pNames))
nreg_snps <- nrow(df_feat)
df_feat <- data.frame(val = as.vector(df_feat))
df_feat$feature <- unlist(lapply(pNames, function(x) rep(x, nreg_snps)))
df_feat$feature <- factor(df_feat$feature, levels = pNames)

totdensity_feat <- ggplot(df_feat) + geom_density(mapping = aes(x = val, fill = feature, group = feature), position = "stack", alpha = 0.5)+ 
  ylab('density') +  xlab('Contribution to prior coeff')+ theme_bw() + theme(legend.position = "none")+scale_fill_d3('category20')

file_name <- sprintf('%splot_priorFeatContr_dist', outFold)
ggsave(sprintf('%s.png',file_name), totdensity_feat, width = 7, height = 5, units = "in")
ggsave(sprintf('%s.pdf',file_name), totdensity_feat, width = 7, height = 5, units = "in")

################################################
#### bar plot max contribution and weigths #####
################################################

tmp <- do.call(rbind, pDat)
mean_feat <- apply(tmp, 2, mean)
max_feat <- apply(tmp, 2, max)
df_contr <- data.frame(feature = pNames, max = max_feat*weights, mean = mean_feat*weights, w = weights)
df_contr$feature <- factor(df_contr$feature, levels = pNames)

pl_w <- ggplot(df_contr, aes(x = feature, y = max, fill = feature) ) + 
  geom_bar(stat="identity", color = 'black') + 
  coord_flip()+
  xlab('') +  theme_bw() + ylab('') + theme(legend.position = "none") + ggtitle('Maximum contribution to prior coefficient')+
  scale_fill_d3('category20')

file_name <- sprintf('%s/plot_maxcontr_weights.png', outFold)
ggsave(file_name,  pl_w, width = 6, height = 2+0.1*length(pNames), units = "in")

pl_w <- ggplot(df_contr, aes(x = feature, y = w, fill = feature) ) + 
  geom_bar(stat="identity", color = 'black') + 
  coord_flip()+
  xlab('') +  theme_bw() + ylab('') + theme(legend.position = "none") + ggtitle('Weights associated to prior feature')+
  scale_fill_d3('category20')

file_name <- sprintf('%splot_weights_priorFeat', outFold)
ggsave(sprintf('%s.png',file_name),  pl_w, width = 6, height = 2+0.1*length(pNames), units = "in")
ggsave(sprintf('%s.pdf',file_name),  pl_w, width = 6, height = 2+0.1*length(pNames), units = "in")

##############################################
#### iteration weigths for prior feature #####
##############################################

if(id_con == id_opt){
  res_p_her_tot_allchr_it <- get(load(sprintf('%sresPrior_EOpt_Iteration_HeritableGenes_allchr.RData', part3Res_fold)))
}else{
  res_p_her_tot_allchr_it <- get(load(sprintf('%sresPrior_EFixed%.2f_Iteration_HeritableGenes_allchr.RData', part3Res_fold, Epar[id_con])))
}

it_weights <- as.matrix(res_p_her_tot_allchr_it$pWeight, ncol = length(pNames))
n_it <- nrow(it_weights)

df_it <- data.frame(w = as.vector(it_weights))
df_it$it <- rep(1:n_it, ncol(it_weights))
df_it$prior <- as.vector(sapply(pNames, function(x) rep(x, n_it)))
df_it$prior <- factor(df_it$prior, levels = pNames )


pl_it <- ggplot(df_it, aes(x = it, y = w, color = prior)) + 
  geom_point(size=0.7, alpha=0.8) + geom_line(alpha = 0.8)+
  theme_classic() + xlab('iteration')+  
  ylab('weights') + theme(legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  scale_color_d3('category20')+labs(color = "")

file_name <- sprintf('%splot_iteration_weigths_priorFeat', outFold)
ggsave(filename=sprintf('%s.png',file_name), plot = pl_it, width = 6, height = 4, dpi=500)
ggsave(filename=sprintf('%s.pdf',file_name), plot = pl_it, width = 6, height = 4, dpi=500)


##########################################
#### fraction of reg SNPs with prior #####
##########################################

id_snps_prior_sep <- apply(tmp, 2, function(x) x!=0)

id_snps_gene_and_prior_p_sep <- list()
id_snps_gene_and_prior_nop_sep <- list()

for(i in 1:length(pNames)){
  print(i)
  id_snps_gene_and_prior_p_sep[[i]] <- id_snps_prior_sep[,i] & unlist(id_snps_gene_p)
  id_snps_gene_and_prior_nop_sep[[i]] <- id_snps_prior_sep[,i] & unlist(id_snps_gene_nop)
}


df_prior <- data.frame(Prior = c('Union', pNames), n_snps_reg = 0, n_snps_with_prior_and_gene = 0, frac = 0)

df_prior[1,-1] <- c(length(which(unlist(id_snps_gene_p))), length(which(unlist(id_snps_gene_and_prior_p))), 
                    length(which(unlist(id_snps_gene_and_prior_p)))/length(which(unlist(id_snps_gene_p))))

df_prior[-1, 2] <- length(which(unlist(id_snps_gene_p)))
df_prior[-1, 3] <- sapply(id_snps_gene_and_prior_p_sep, function(x) length(which(x)))
df_prior[-1, 4] <- df_prior[-1, 3]/df_prior[-1, 2]

df_prior$weights <- c(NA, max_feat*weights)

### same but nop regression 

df_prior$n_snps_reg_nop <- length(which(unlist(id_snps_gene_nop)))
df_prior$n_snps_with_prior_and_gene_nop <- c(length(which(unlist(id_snps_gene_and_prior_nop))), sapply(id_snps_gene_and_prior_nop_sep, function(x) length(which(x))))
df_prior$frac_nop <-df_prior$n_snps_with_prior_and_gene_nop/df_prior$n_snps_reg_nop 
df_prior <- df_prior[-1,]
df_prior$Prior <- factor(df_prior$Prior, levels = pNames)
write.table(file = sprintf('%sfractionPrior_perFeature.txt', outFold),x = df_prior, col.names = T, row.names = F, sep ='\t', quote = F)

library(latex2exp)
library(cowplot)
library(ggpubr)
library(ggrepel)

# make plot
pl_prior <- ggplot(df_prior, aes(x = frac_nop, y = frac, color = Prior, size = weights)) + 
  geom_point(alpha = 0.6) + geom_text_repel(label =df_prior$Prior, size=3, segment.size = 0.1)+
  xlim(0, max(df_prior$frac))+ylim(0,max(df_prior$frac))+
  geom_abline(slope = 1, intercept = 0, linetype=2, alpha = 0.6)+
  geom_abline(slope = 2, intercept = 0, linetype=2, alpha = 0.6, color = 'red')+
  xlab(TeX('$\\frac{n. reg-SNPs\\,with\\,prior}{n. reg-SNPs}$ el-net'))  +  theme_classic() + 
  ylab(TeX('$\\frac{n. reg-SNPs\\,with\\,prior}{n. reg-SNPs}$ el-net-prior')) + theme(legend.position = c(0.75, 0.3), plot.title = element_text(hjust = 0.5))+
  scale_color_d3('category20')+labs(size = "Max contribution \nto prior coefficient")+
  guides(color =FALSE)

tmp <- ggplot(df_prior, aes(x = frac_nop, y = frac, color = Prior, size = weights)) + 
  geom_point(alpha = 0.6) + 
  xlim(0, max(df_prior$frac))+ylim(0,max(df_prior$frac))+
  geom_abline(slope = 1, intercept = 0, linetype=2, alpha = 0.6)+
  geom_abline(slope = 2, intercept = 0, linetype=2, alpha = 0.6, color = 'red')+
  xlab(TeX('$\\frac{n. reg-SNPs\\,with\\,prior}{n. reg-SNPs}$ el-net')) +  theme_classic() +
  ylab(TeX('$\\frac{n. reg-SNPs\\,with\\,prior}{n. reg-SNPs}$ el-net-prior')) + theme(legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  scale_color_d3('category20')+labs(size = "Max contribution \nto prior coefficient")+
  guides(size =FALSE)

leg <- get_legend(tmp)
pl_leg <- as_ggplot(leg)

file_name <- sprintf('%splot_priorFraction_elVSelP', outFold)
ggsave(sprintf('%s.png',file_name),   pl_prior, width = 5.1, height = 5.1, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name),   pl_prior, width = 5.1, height = 5.1, units = "in", dpi = 500)

file_name <- sprintf('%splot_priorFraction_elVSelP_legend', outFold)
ggsave(sprintf('%s.png',file_name),   pl_leg, width = 5.1, height = 5.1, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name),   pl_leg, width = 5.1, height = 5.1, units = "in", dpi = 500)

###################################################
#### distribution of regSNPs with prior coeff #####
###################################################

abs_beta_reg_nop <- unlist(lapply(beta_snps_nop, function(x) rowSums(abs(x))))
abs_beta_reg_p <- unlist(lapply(beta_snps_p, function(x) rowSums(abs(x))))

abs_beta_reg_prior_sep_nop <- lapply(id_snps_gene_and_prior_nop_sep, function(x) abs_beta_reg_nop[x])
abs_beta_reg_prior_sep_p <- lapply(id_snps_gene_and_prior_p_sep, function(x) abs_beta_reg_p[x])

# plot distribution in boxplot format
df_abs <- data.frame(val = log10(c(unlist(abs_beta_reg_prior_sep_nop), unlist(abs_beta_reg_prior_sep_p))), 
                     prior = c(unlist(lapply(1:length(pNames), function(x) rep(pNames[x], length(abs_beta_reg_prior_sep_nop[[x]])))), 
                               unlist(lapply(1:length(pNames), function(x) rep(pNames[x], length(abs_beta_reg_prior_sep_p[[x]]))))),
                     type = c(rep('el-net', length(unlist(abs_beta_reg_prior_sep_nop))), rep('el-net-prior', length(unlist(abs_beta_reg_prior_sep_p)))))

df_abs$prior <- factor(df_abs$prior, levels = pNames)
df_abs$type <- factor(df_abs$type)

# make plot
pl_abs <- ggplot(df_abs, aes(x = prior, y = val, fill = type)) + 
  geom_boxplot(outlier.size = 0.2) +
  ylab('log10 abs(beta) reg-SNPs with prior')  +  
  theme_classic() + scale_fill_grey(start = 0.7, end = 0.4) +
  theme(legend.position = 'right', plot.title = element_text(hjust = 0.5),axis.text.x=element_text(angle = 45, hjust = 1))

file_name <- sprintf('%splot_distAbsBeta_regSNPsPrior_elANDelP', outFold)
ggsave(sprintf('%s.png',file_name),   pl_abs, width = 8, height = 5.1, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name),   pl_abs, width = 8, height = 5.1, units = "in", dpi = 500)

# plot increase in term of median
df_abs_med <- data.frame(nop = sapply(abs_beta_reg_prior_sep_nop, function(x) median(log10(x))), p = sapply(abs_beta_reg_prior_sep_p, function(x) median(log10(x))), prior = pNames)
df_abs_med$weights <- max_feat*weights
df_abs_med$prior <- factor(df_abs_med$prior, levels = pNames)


# make plot
pl_med <- ggplot(df_abs_med, aes(x = nop, y = p, color = prior, size = weights)) + 
  geom_point(alpha = 0.6) + geom_text_repel(label =df_abs_med$prior, size=3, segment.size = 0.1)+
  xlim(min(c(df_abs_med$nop,df_abs_med$p)), max(c(df_abs_med$nop,df_abs_med$p)))+ylim(min(c(df_abs_med$nop,df_abs_med$p)), max(c(df_abs_med$nop,df_abs_med$p)))+
  ggtitle('log10 abs(beta) reg-SNPs with prior (median)')+
  geom_abline(slope = 1, intercept = 0, linetype=2, alpha = 0.6)+
  xlab('el-net')  +  theme_classic() + 
  ylab('el-net-prior') + theme(legend.position = c(0.75, 0.3), plot.title = element_text(hjust = 0.5))+
  scale_color_d3('category20')+labs(size = "Max contribution \nto prior coefficient")+
  guides(color =FALSE)

file_name <- sprintf('%splot_AbsBeta_regSNPsPrior_elVSelP', outFold)
ggsave(sprintf('%s.png',file_name), pl_med, width = 5.1, height = 5.1, units = "in", dpi = 500)
ggsave(sprintf('%s.pdf',file_name), pl_med, width = 5.1, height = 5.1, units = "in", dpi = 500)


