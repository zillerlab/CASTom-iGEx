# cluster 

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
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlist))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering using PG method, merge loci")
parser$add_argument("--inputFile", type = "character", default = NULL, help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--color_file", type = "character", help = "file with color based on phenotype")
parser$add_argument("--covDatFile", type = "character", default = NULL, help = "additional cov to test")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character", help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--pval_thr", type = "double", default = 1, help = "threshold to filter features")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED")
parser$add_argument("--kNN_par", type = "integer", nargs = '*', default = 30, help = "parameter used for PG method")
parser$add_argument("--geneRegionFile", type = "character", default=NULL, help = "used if tscore and exclude_MHC")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pvalresFile <- args$pvalresFile
tissues_name <- args$tissues_name
color_file <- args$color_file
pval_id <- args$pval_id
inputFile <- args$inputFile
geneRegionFile <- args$geneRegionFile
covDatFile <- args$covDatFile
sampleAnnFile <- args$sampleAnnFile
type_cluster <- args$type_cluster
split_tot <- args$split_tot
functR <- args$functR
type_sim <- args$type_sim
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
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software//model_clustering/clustering_functions.R'
# corr_thr <- 0.9
# type_input <- 'zscaled'
# kNN_par <- 30
# color_file <- '/psycl/g/mpsziller/lucia/castom-igex/Figures/color_tissues.txt'
# tissues_name <- 'Liver'
# pval_thr <- 0.05
# capped_zscore <- F
####################################################################################################################

source(functR)

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
res_pval <- res_pval$tscore[[pval_id]]
id_pval <- 8
id_info <- 2

# filter based on pvalue thr:
res_pval <- res_pval[res_pval[,id_pval] <= pval_thr,]
print(dim(res_pval))

# load input matrix 
if(split_tot == 0){
  
  scoreMat <- get(load(inputFile))
  # filter out based on samples and ids
  id_el <- intersect(scoreMat[,1], res_pval[, id_info])
  scoreMat <- scoreMat[match(id_el,scoreMat[,1]), ]
  
  common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
  
  rownames(scoreMat) <- common_samples
  colnames(scoreMat) <- id_el
  res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
  
}else{
  
  ###### load score Mat #######
  scoreMat_list <- vector(mode = 'list', length = split_tot)
  samplesID <- vector(mode = 'list', length = split_tot)
  elementID <- NULL
  
  for(i in 1:split_tot){
    
    print(i)
    if(file.exists(sprintf('%s%i.RData', inputFile, i))){
      tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
      elementID <- c(elementID,tmp[,1])
      samplesID[[i]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
      scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
    }else{
      print(sprintf('split %i does not exist', i))
      split_tot <- split_tot - 1
    }
  }
  
  print(split_tot)
  # check samplesID always the same
  if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
  
  scoreMat <- do.call(cbind, scoreMat_list)
  colnames(scoreMat) <- elementID
  rm(scoreMat_list)
  
  # filter out elements that are repeated twice:
  id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
  scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
  
  id_el <- intersect(colnames(scoreMat),  res_pval[, id_info])
  scoreMat <- scoreMat[, match(id_el, colnames(scoreMat))]
  
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
  res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
  
}
print(identical(colnames(scoreMat), res_pval[, id_info]))

input_data <- scale(scoreMat)
attr(input_data, "scaled:scale") <- NULL
attr(input_data, "scaled:center") <- NULL

  
tmp_val <- res_pval[, id_pval-1]
input_data <- sapply(1:ncol(input_data), 
                       function(x) input_data[, x]*tmp_val[res_pval[,id_info] == colnames(input_data)[x]])
colnames(input_data) <- colnames(scoreMat)
  
# merge loci
gene_info <- fread(geneRegionFile, data.table = F, stringsAsFactors = F, h=T)
gene_info <- gene_info[gene_info$ensembl_gene_id %in% res_pval$ensembl_gene_id, 1:9]
gene_loci <- merge_loci_genes(gene_table = gene_info, cis_size = 200000, 
                              bp_loci = 1000000, tissue = tissues_name)

input_data_loci <- matrix(nrow = nrow(input_data), ncol = nrow(gene_loci))
for(i in 1:ncol(input_data_loci)){
  genes <- strsplit(gene_loci$gene[i], split = ",")[[1]]
  input_data_loci[,i] <- rowMeans(input_data[,genes, drop=F])
}
rownames(input_data_loci) <- rownames(input_data)
colnames(input_data_loci) <- gene_loci$loci_id

###############################################################
## pheno graph clustering:
ed_dist <- as.matrix(dist(input_data_loci, method = 'euclidian'))^2

if(type_sim == 'HK'){fun_cl <- clust_PGmethod_HKsim}
if(type_sim == 'ED'){fun_cl <- clust_PGmethod_EDdist}
if(!type_sim %in% c('ED', 'HK')){stop('type similarity nust be ED or HK')}

# NOTE: for cycle can be parallelized
PG_cl <- vector(mode = 'list', length = length(kNN_par))
test_cov <- vector(mode = 'list', length = length(kNN_par))
for(i in 1:length(kNN_par)){
  
  PG_cl[[i]] <- fun_cl(kNN = kNN_par[i], score = input_data_loci, 
                       type_Dx = type_cluster, sample_info=sampleAnn, euclDist=ed_dist)
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
opt_k <- kNN_par[which.max(info_hyperParam$DB_mean)]

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

output <- list(best_k = opt_k, cl_res = PG_cl, test_cov = test_cov, info_tune = info_hyperParam, feat = colnames(input_data), res_pval = res_pval,
               cl_best = data.frame(id = sampleAnn$Individual_ID, gr = PG_cl[[which.max(info_hyperParam$DB_mean)]]$cl$membership))
output$Dx_perc <- list(perc = df_perc, test = df_perc_test)
output$samples_id <- rownames(input_data_loci)
# most significant elements
test_diff <- data.frame(id = colnames(input_data_loci), pval = apply(input_data_loci, 2, function(x) kruskal.test(x = x, g = factor(output$cl_best$gr))$p.value))
test_diff$pval_corr <- p.adjust(test_diff$pval, method = 'BH')
output$test_diff_gr <- test_diff
# compute mean for each gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_best$gr)), nrow = nrow(test_diff))

df_gr_mean[,] <- t(apply(input_data_loci, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) mean(x[output$cl_best$gr == y]) )))

df_gr_sd[,] <- t(apply(input_data_loci, 2, function(x) 
  sapply(sort(unique(output$cl_best$gr)), function(y) sd(x[output$cl_best$gr == y]) )))

df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_best$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(input_data_loci)

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)
output$input_data <- input_data_loci
output$gene_loci <- gene_loci

# save results:
save(output, file = sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric.RData', 
                            outFold, type_cluster, type_sim))

### plot: UMAP
n_comp_umap <- 2
n_neigh_umap <- opt_k
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap

umap_res <- umap::umap(input_data_loci, custom.settings)

df <- data.frame(component_1=umap_res$layout[,1], component_2=umap_res$layout[,2], gr = output$cl_best$gr)
df$gr <- factor(df$gr)

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr))+
  geom_point(size = 0.05)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4
if(type_cluster == 'All'){
  df$Dx <- factor(sampleAnn$Dx)
  pl_extra <- ggplot(df, aes(x = component_1, y=component_2, color = Dx))+
    geom_point(size = 0.05)+
    theme_bw()+theme(legend.position = 'right')
  tot_pl <- ggarrange(plotlist = list(tot_pl, pl_extra), ncol = 2, nrow = 1, align='h')
  width_pl <- 8
}

ggsave(filename = sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric_umap.png', 
                          outFold, type_cluster, type_sim), 
       width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric_umap.pdf', 
                          outFold, type_cluster, type_sim), 
       width = width_pl, height = 4, plot = tot_pl, device = 'pdf')

df$Gender <- factor(sampleAnn$Gender)
df$Batch <- factor(sampleAnn$Batch)
df$Array <- factor(sampleAnn$Array)
df$Centre <- factor(sampleAnn$initial_assessment_centre)
df$Age <- sampleAnn$Age

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(sampleAnn$Age), max(sampleAnn$Age)))

pl1 <- ggplot(df, aes(x = component_1, y=component_2, color = Gender))+
  geom_point(size = 0.03)+ggtitle('Gender')+
  theme_bw()+theme(legend.position = 'none')
pl2 <- ggplot(df, aes(x = component_1, y=component_2, color = Batch))+
  geom_point(size = 0.03)+ggtitle('Batch')+
  theme_bw()+theme(legend.position = 'none')
pl3 <- ggplot(df, aes(x = component_1, y=component_2, color = Array))+
  geom_point(size = 0.03)+ggtitle('Array')+
  theme_bw()+theme(legend.position = 'none')
pl4 <- ggplot(df, aes(x = component_1, y=component_2, color = Centre))+
  geom_point(size = 0.03)+ggtitle('Centre')+
  theme_bw()+theme(legend.position = 'none')
pl5 <- ggplot(df, aes(x = component_1, y=component_2, color = Age))+
  geom_point(size = 0.03)+ggtitle('Age')+sc+
  theme_bw()+theme(legend.position = 'right')
tot_pl <- ggarrange(plotlist = list(pl1, pl2, pl3, pl4, pl5), ncol = 3, nrow = 2)
ggsave(filename = sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric_umap_cov.png', 
                          outFold, type_cluster, type_sim), 
       width = 12, height = 8, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric_umap_cov.pdf',
                          outFold, type_cluster, type_sim), 
       width = 12, height = 8, plot = tot_pl, device = 'pdf')

### plot: heatmap
file_name <- sprintf('%stscore_zscaled_loci_cluster%s_PGmethod_%smetric', 
                     outFold,type_cluster, type_sim)

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name, color_tissues$tissues),]

tmp <- test_diff[test_diff$pval_corr < 0.05,]
if(nrow(tmp)>50){
  keep_feat <- test_diff[order(test_diff$pval)[1:50],]
}else{
  keep_feat <- tmp
}

mat <- output$gr_input$cv[output$gr_input$cv$id %in%  keep_feat$id, ]
width_pl <- 10

# remove duplicated
mat <- mat[!duplicated(mat$id),]
mat$tissue <- tissues_name

pheat_pl_gr(mat, type_mat = 'tscore', height_pl = 7, 
            width_pl = width_pl, color_df = color_tissues, outFile = paste0(file_name, '_heatmap_gr'))


