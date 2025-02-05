#!/usr/bin/env Rscript
# cluster, correct for PCs, single cohort

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(clustAnalytics))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering using PG method")
parser$add_argument("--inputFile", type = "character", default = NULL, help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--covDatFile", type = "character", default = NULL, help = "additional cov to test")
parser$add_argument("--cluster_method", type = "character", default = "leiden", help = "leidein or louvain, community detection method")
parser$add_argument("--kNN_par", type = "integer", nargs = '*', default = 20, help = "parameter used for PG method")
parser$add_argument("--type_cluster", type = "character", help = "All, Cases, Controls")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character", help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 1, help = "id to be used on pvalue file")
parser$add_argument("--pval_thr", type = "double", default = 1, help = "threshold to filter features")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", default = "tscore", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--exclude_MHC", type = "logical", default = F, help = "if true, MHC region excluded (only ossible for tscore)")
parser$add_argument("--capped_zscore", type = "logical", default = F, help = "if true, zstat is capped based on distribution")
parser$add_argument("--geneRegionFile", type = "character", default = NULL, help = "used if tscore and exclude_MHC")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pvalresFile <- args$pvalresFile
cluster_method <- args$cluster_method
tissues_name <- args$tissues_name
capped_zscore <- args$capped_zscore
pval_id <- args$pval_id
inputFile <- args$inputFile
geneRegionFile <- args$geneRegionFile
exclude_MHC <- args$exclude_MHC
covDatFile <- args$covDatFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
split_tot <- args$split_tot
functR <- args$functR
type_data <- args$type_data
corr_thr <- args$corr_thr
min_genes_path <- args$min_genes_path
type_sim <- args$type_sim
type_input <- args$type_input
pval_thr <- args$pval_thr
kNN_par <- args$kNN_par
outFold <- args$outFold

#####################################################################################################################
# inputFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/predictedTscores_splitGenes'
# sampleAnnFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All.txt'
# covDatFile <- NULL
# split_tot <- 100
# pvalresFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData'
# pval_id <- 1
# min_genes_path <- 2
# type_data <- 'tscore'
# type_cluster <- 'Cases'
# type_sim <- 'HK'
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software//model_clustering/clustering_functions.R'
# corr_thr <- 0.9
# type_input <- 'zscaled'
# kNN_par <- 20
# tissues_name <- 'Liver'
# cluster_method <- 'leiden'
#####################################################################################################################

source(functR)
print(paste("kNN parameters:", kNN_par))

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)

if(type_cluster == 'Cases'){
  sampleAnn <- sampleAnn[sampleAnn$Dx == 1,]
}else{
  if(type_cluster == 'Controls'){
    sampleAnn <- sampleAnn[sampleAnn$Dx == 0,]
  }else{
    if(type_cluster != 'All')
      stop('type_cluster must be either Cases or Controls or All')
  }
}

if(!is.null(covDatFile)){
  covDat <- read.table(covDatFile, h=T, stringsAsFactors = F, sep = '\t')
  covDat <- covDat[match(sampleAnn$Individual_ID, covDat$Individual_ID), ]
  sampleAnn$Batch <- covDat$Batch
  sampleAnn$Array <- covDat$Array
  sampleAnn$initial_assessment_centre <- covDat$initial_assessment_centre
}  

# load pval res
res_pval <- get(load(pvalresFile))
if(type_data == 'tscore'){
  res_pval <- res_pval$tscore[[pval_id]]
  id_pval <- 8
  id_info <- 2
  id_geno_summ <- 3
}else{
  if(type_data == 'path_Reactome'){
    res_pval <- res_pval$pathScore_reactome[[pval_id]]
    id_pval <- 13
    id_info <- 1
    id_geno_summ <- 4
  }else{
    if(type_data == 'path_GO'){
      res_pval <- res_pval$pathScore_GO[[pval_id]]
      id_pval <- 15
      id_info <- 1
      id_geno_summ <- 6
    }else{
      stop('unknown pathway called')
    }
  }
}

# recompute pvalue if ngenes_tscore > 1
if(min_genes_path > 1 & grepl('path',type_data)){
  res_pval <- res_pval[res_pval$ngenes_tscore >= min_genes_path, ]
  res_pval[,id_pval+1] <- qvalue(res_pval[,id_pval])$qvalues
  res_pval[,id_pval+2] <- p.adjust(res_pval[,id_pval], method = 'BH')
}

if(exclude_MHC & type_data == 'tscore'){
  res_pval$start_position <- NA
  res_pval$chrom <- NA
  tmp <- read.table(geneRegionFile, h=T,stringsAsFactors = F)
  tmp <- tmp[match(res_pval$ensembl_gene_id, tmp$ensembl_gene_id),]
  res_pval$start_position <- tmp$start_position
  res_pval$chrom <- tmp$chrom
  HLA_reg <- c(26000000, 34000000)
  res_pval <- res_pval[!(res_pval$chrom %in% 'chr6' & res_pval$start_position <=HLA_reg[2] & res_pval$start_position >= HLA_reg[1]) , ]
}

# filter based on pvalue thr:
res_pval <- res_pval[res_pval[,id_pval+2] <= pval_thr,]
print(dim(res_pval))

#### load input matrix ####
load_output <- load_input_matrix(inputFile = inputFile, 
                              sampleAnn = sampleAnn, 
                              res_pval = res_pval, 
                              split_tot = split_tot, 
                              id_info = id_info)

scoreMat <- load_output$scoreMat
res_pval <- load_output$res_pval
sampleAnn <- load_output$sampleAnn

print(identical(colnames(scoreMat), res_pval[, id_info]))

### clumping: sort according best SNP association ###
if(corr_thr < 1){
  cor_score <- cor(scoreMat)
  element_rm <- clumping_features(res_pval=res_pval, 
                                id_info = id_info, 
                                corr_feat = cor_score, 
                                id_pval = id_pval, 
                                corr_thr = corr_thr)

  print(paste(length(element_rm),'features removed due to high correlation'))
  scoreMat <- scoreMat[,!colnames(scoreMat) %in% element_rm]
  res_pval <- res_pval[match(colnames(scoreMat), res_pval[, id_info]),]
}else{
  print('All features considered')
}

input_data_notcorr <- scale(scoreMat)
attr(input_data_notcorr, "scaled:scale") <- NULL
attr(input_data_notcorr, "scaled:center") <- NULL

# remove PCs1-10 for each genes
input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
rownames(input_data) <- rownames(input_data_notcorr)
colnames(input_data) <- colnames(input_data_notcorr)
name_cov <- setdiff(colnames(sampleAnn),
                    c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex',
                      'Age', 'Temp_ID', 'cohort', 'cohort_id'))
fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))
for(i in 1:ncol(input_data_notcorr)){
  # print(i)
  tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn[, name_cov])
  reg <- lm(fmla, data = tmp)
  input_data[,i] <- reg$residuals
}
print("corrected for PCs")

if(type_input == 'zscaled'){
  
  # cap zscores (cap values)
  if(capped_zscore){
    tmp_val <- res_pval[, id_pval-1]
    thr_val <- quantile(tmp_val, probs=c(0.005, 0.995))
    tmp_val[tmp_val<=thr_val[1]] <- thr_val[1]
    tmp_val[tmp_val>=thr_val[2]] <- thr_val[2]
  }else{
    tmp_val <- res_pval[, id_pval-1]
  }
  input_data <- sapply(1:ncol(input_data), function(x) 
    input_data[, x]*tmp_val[res_pval[,id_info] == colnames(input_data)[x]])
  
  colnames(input_data) <- colnames(scoreMat)
}

###############################################################
## pheno graph clustering:
ed_dist <- as.matrix(dist(input_data, method = 'euclidian'))^2

if(type_sim == 'HK'){fun_cl <- clust_PGmethod_HKsim}
if(type_sim == 'ED'){fun_cl <- clust_PGmethod_EDdist}
if(!type_sim %in% c('ED', 'HK')){stop('type similarity nust be ED or HK')}

# NOTE: for cycle can be parallelized
PG_cl <- vector(mode = 'list', length = length(kNN_par))
test_cov <- vector(mode = 'list', length = length(kNN_par))
for(i in 1:length(kNN_par)){
  
  PG_cl[[i]] <- fun_cl(kNN = kNN_par[i], 
                       score = input_data, 
                       type_Dx = type_cluster, 
                       sample_info=sampleAnn,
                       euclDist = ed_dist, 
                       cluster_method = cluster_method)
  
  print(PG_cl[[i]]$info)
  # cluster depend on PC?
  id <- PG_cl[[i]]$cl$membership
  df <- cbind(data.frame(cl = id), sampleAnn[,! colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID')])
  test_cov[[i]] <- data.frame(cov_id = colnames(df)[-(1)])
  test_cov[[i]]$test_type <- test_cov[[i]]$statistic <- test_cov[[i]]$pval <- NA
  for(j in 1:(ncol(df)-1)){
    
    if(is.integer(df[, j+1]) | is.character(df[, j+1])){
      tmp <- chisq.test(table(df$cl, df[, j+1]))
      test_cov[[i]]$pval[j] <- tmp$p.value
      test_cov[[i]]$statistic[j] <-tmp$statistic
      test_cov[[i]]$test_type[j] <- 'chisq'
    }else{
      tmp <- kruskal_test(df[, j+1] ~ factor(df$cl))
      test_cov[[i]]$pval[j] <- pvalue(tmp)
      test_cov[[i]]$statistic[j] <- tmp@statistic@teststatistic
      test_cov[[i]]$test_type[j] <- 'kruskal'
    }
  }
  test_cov[[i]]$pval_BHcorr <- p.adjust(test_cov[[i]]$pval, method = 'BH')
  test_cov[[i]]$kNN <- kNN_par[i]
  print(test_cov[[i]])
}

test_cov <- do.call(rbind, test_cov)
info_hyperParam <- do.call(rbind, lapply(PG_cl, function(x) x$info))
opt_k <- kNN_par[which.max(info_hyperParam$coverage_and_conductance)]

# if type_clster == 'All' compute percentage for each group
df_perc <- df_perc_test <- list()
if(type_cluster == 'All'){
  
  for(i in 1:length(kNN_par)){
    perc <- table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)/rowSums(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))
    # test fisher for each group
    cl_id <- sort(unique(PG_cl[[i]]$cl$membership))
    df_perc[[i]] <- data.frame(gr = rep(cl_id,2), Dx = c(rep(0,nrow(perc)), rep(1,nrow(perc))), 
                               perc = as.vector(perc), count = as.vector(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx)))
    df_perc_test[[i]] <- data.frame(gr = c(cl_id, 'all'), 
                                    fisher_test = c(sapply(cl_id, function(x) fisher.test(table(PG_cl[[i]]$cl$membership == x, sampleAnn$Dx))$p.value), 
                                                    chisq.test(table(PG_cl[[i]]$cl$membership, sampleAnn$Dx))$p.value))
  }
}

output <- list(best_k = opt_k, cl_res = PG_cl, test_cov = test_cov, 
               info_tune = info_hyperParam, feat = colnames(input_data),
               res_pval = res_pval,
               cl_best = data.frame(id = sampleAnn$Individual_ID, 
               gr = PG_cl[[which.max(info_hyperParam$coverage_and_conductance)]]$cl$membership))
output$Dx_perc <- list(perc = df_perc, test = df_perc_test)
output$samples_id <- rownames(input_data)

# most significant elements
test_diff <- data.frame(id = colnames(input_data), pval = apply(input_data, 2, function(x) kruskal.test(x = x, g = factor(output$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
output$test_diff_gr <- test_diff
# compute mean for each gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))

df_gr_mean[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) mean(x[output$cl_best$gr == y]) )))

df_gr_sd[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) sd(x[output$cl_best$gr == y]) )))

df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(input_data)

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
output$input_data <- input_data
output$ed_dist <- ed_dist

# save results:
save(output, file = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric.RData', 
                            outFold, type_data, type_input, type_cluster, type_sim))

### plot: UMAP ###
n_comp_umap <- 2
n_neigh_umap <- opt_k
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap

umap_res <- umap::umap(input_data, custom.settings)

df <- data.frame(component_1=umap_res$layout[,1], 
                 component_2=umap_res$layout[,2], 
                 gr = output$cl_best$gr)
df$gr <- factor(df$gr)
P <- length(unique(df$gr))
gr_color <- pal_d3(palette = 'category20')(P)

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr))+
  geom_point(size = 0.05, alpha = 0.8)+
  xlab('UMAP component 1')+ ylab('UMAP component 2')+
  scale_color_manual(values = gr_color)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4

ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap.png', outFold, type_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_umap.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')

# plot distribution PCs
df_cov <- cbind(sampleAnn, data.frame(gr = df$gr))
df_cov$gr <- paste0('gr', df_cov$gr)
df_cov$gr <- factor(df_cov$gr, levels = paste0('gr', 1:P))
df_cov_PC <- data.frame(val = as.vector(as.matrix(df_cov[, paste0('PC', 1:10)])), 
                        PC = unlist(lapply(paste0('PC', 1:10), function(x) rep(x, nrow(df_cov)))), 
                        gr = rep(df_cov$gr, 10))
df_cov_PC$PC <- factor(df_cov_PC$PC, levels = paste0('PC', 1:10))

p <- ggboxplot(df_cov_PC, x = "gr", y = "val", fill = "gr", color = 'black', legend = 'none', outlier.size = 0.2, alpha = 0.8) + stat_compare_means(label = "p.format", size = 3) 
p <- ggpar(p, palette = gr_color, xlab = '', ylab = '', x.text.angle = 45)
p <- facet(p, facet.by = "PC", short.panel.labs = T, scales = 'free_y', nrow = 1)

ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_PCs.png', outFold, type_data, type_input, type_cluster, type_sim), width = 13, height = 4, plot = p, device = 'png')
ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_PCs.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 13, height = 4, plot = p, device = 'pdf')

# plot distribution AGE, Sex if the information is available
if ("Age" %in% colnames(df_cov) & "Gender" %in% colnames(df_cov)) {
  df_cov_age <- data.frame(Age = df_cov$Age, gr = df_cov$gr)
  output$test_cov$pval[output$test_cov$cov_id == 'Age'] <- kruskal.test(df_cov_age$Age, g = df_cov_age$gr)$p.value

  df_cov_sex <- data.frame(n = as.vector(table(df_cov$Gender, df_cov$gr)), Sex = rep(c('male', 'female'), P), gr = unlist(lapply(paste0('gr', 1:P), function(x) rep(x, 2))))
  df_cov_sex$Sex <- factor(df_cov_sex$Sex, levels = c('male', 'female'))

  pl_a <- ggplot(df_cov_age, aes(x = gr, y = Age, fill = gr))+
    geom_violin(alpha = 0.8)+
    geom_boxplot(width=0.2, fill="white")+
    xlab('')+ ylab('Age')+
    scale_fill_manual(values = gr_color)+
    annotate("text", x = 1, y = max(df_cov_age$Age)+2,
      label = sprintf('p=%s',	as.character(round(output$test_cov$pval[output$test_cov$cov_id == 'Age'], digits = 2))))+
    theme_bw()+theme(legend.position = 'none')

  pl_s <- ggplot(df_cov_sex, aes(x = gr, y = n, color = gr, fill = Sex))+
    geom_bar(size = 1, alpha = 0.8, stat = 'identity', position = position_dodge())+
    xlab('')+ ylab('number of individual')+
    scale_color_manual(values = gr_color)+
    scale_fill_manual(values = c('grey10', 'grey60'))+
    guides(color = FALSE)+
    annotate("text", x = 1, y = max(df_cov_sex$n)+2,
      label = sprintf('p=%s', as.character(round(output$test_cov$pval[output$test_cov$cov_id == 'Gender'], digits = 2))))+
    theme_bw()+theme(legend.position = 'right')

  tot_pl <- ggarrange(plotlist = list(pl_a, pl_s), ncol = 2, nrow = 1, align='h', widths=c(1, 1.4))
  ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_Age_Sex.png', outFold, type_data, type_input, type_cluster, type_sim), width = 6.5, height = 3, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%s%s_corrPCs_%s_cluster%s_PGmethod_%smetric_Age_Sex.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 6.5, height = 3, plot = tot_pl, device = 'pdf')
}
