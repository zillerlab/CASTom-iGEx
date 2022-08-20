# predict cluster probability based on phenograp approach

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
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(SNFtool))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="predict cluster probability for new samples")
parser$add_argument("--sampleAnnNew_file", type = "character", help = "")
parser$add_argument("--name_cohort", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--inputFile", type = "character", help = "input files (scores)")
parser$add_argument("--clustFile", type = "character", help = "file cluster results")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
name_cohort <- args$name_cohort
sampleAnnNew_file <- args$sampleAnnNew_file
clustFile <- args$clustFile
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
split_tot <- args$split_tot
type_cluster <- args$type_cluster
tissues_name <- args$tissues_name
outFold <- args$outFold


###################################################################################################################
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All.txt'
# sampleAnnNew_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/CG/covariateMatrix.txt'
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_clustering/clustering_functions.R'
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# clustFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscore_corrPCs_zscaled_clusterCases_PGmethod_HKmetric.RData'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/CG/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/'
# tissues_name <- 'Liver'
# color_file <- '/psycl/g/mpsziller/lucia/castom-igex/Figures/color_tissues.txt'
# split_tot = 0
# inputFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/CG/devgeno0.01_testdevgeno0/predictedTscores.txt'
# ###################################################################################################################

source(functR)

sampleAnn_new <- read.table(sampleAnnNew_file, h=T, stringsAsFactors = F)
name_cov <- setdiff(colnames(sampleAnn_new),c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex', 'Age', 'Array'))
if(!is.null(name_cohort)){
  sampleAnn_new$cohort <- name_cohort
}

# load model
clust_res <- get(load(clustFile))
cl <- clust_res$cl_best
sampleAnn <- clust_res$sampleInfo

if(type_cluster == 'Cases'){
  sampleAnn <- sampleAnn[sampleAnn$Dx == 1,]
  sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 1,]
}else{
  if(type_cluster == 'Controls'){
    sampleAnn <- sampleAnn[sampleAnn$Dx == 0,]
    sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 0,]
  }else{
    if(type_cluster != 'All')
      stop('type_cluster must be either Cases or Controls or All')
  }
}

sampleAnn_tot <- rbind(sampleAnn[, intersect(colnames(sampleAnn), colnames(sampleAnn_new))], 
                       sampleAnn_new[, intersect(colnames(sampleAnn), colnames(sampleAnn_new))])

res_pval <- clust_res$res_pval
if(type_data == 'tscore'){
  id_pval <- 8
  id_info <- 2
  id_geno_summ <- 3
}else{
  if(type_data == 'path_Reactome'){
    id_pval <- 13
    id_info <- 1
    id_geno_summ <- 4
  }else{
    if(type_data == 'path_GO'){
      id_pval <- 15
      id_info <- 1
      id_geno_summ <- 6
    }else{
      stop('unknown pathway called')
    }
  }
}

# load input matrix (new data, old data input saved in clustering result)
#### load input matrix ####
load_output <- load_input_matrix(inputFile = inputFile, 
                                 sampleAnn = sampleAnn_new, 
                                 res_pval = res_pval, 
                                 split_tot = split_tot, 
                                 id_info = id_info)

scoreMat <- load_output$scoreMat
sampleAnn_new <- load_output$sampleAnn

print(identical(colnames(scoreMat), res_pval[, id_info]))

### scale data ###
input_data_notcorr <- scale(scoreMat)
attr(input_data_notcorr, "scaled:scale") <- NULL
attr(input_data_notcorr, "scaled:center") <- NULL

# remove PCs1-10 for each genes
input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
rownames(input_data) <- rownames(input_data_notcorr)
colnames(input_data) <- colnames(input_data_notcorr)
fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))
for(i in 1:ncol(input_data_notcorr)){
  # print(i)
  tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn_new[, name_cov])
  reg <- lm(fmla, data = tmp)
  input_data[,i] <- reg$residuals
}
print("corrected for PCs")

if(type_input == 'zscaled'){
  input_data <- sapply(1:ncol(input_data), 
                       function(x) input_data[, x]*res_pval[res_pval[,id_info] == colnames(input_data)[x], id_pval-1])
  colnames(input_data) <- colnames(scoreMat)
}

data_tot <- rbind(clust_res$input_data, input_data)
print(identical(rownames(data_tot), sampleAnn_tot$Individual_ID))

#### project new data clustering ####
# compute total ED
ed_dist <-  dist2(as.matrix(data_tot),as.matrix(data_tot))
res_pred <- project_clust_PGmethod_HKsim(kNN = clust_res$best_k, score = data_tot, 
                                         data_mod = clust_res$input_data, 
                                         sample_info = sampleAnn_new, euclDist = ed_dist, 
                                         cl_mod = clust_res$cl_best)

gr_id <- colnames(res_pred$probability)[-(1:2)]
gr_id <- as.numeric(sapply(gr_id, function(x) strsplit(x, split='gr_')[[1]][2]))
cl_new <- data.frame(id = res_pred$probability$Individiual_ID, gr = gr_id[apply(res_pred$probability[, -(1:2)], 1, function(x) which.max(x))])
output <- list(probability = res_pred$probability, tot_W_sNN = res_pred$tot_W_sNN, sampleAnn = sampleAnn_new, data_new = input_data, cl_new = cl_new, res_pval=res_pval)

# compute mean for each gr
test_diff <- clust_res$test_diff_gr
df_gr_mean <- matrix(ncol = length(unique(output$cl_new$gr)), nrow = nrow(test_diff))
df_gr_sd <- matrix(ncol = length(unique(output$cl_new$gr)), nrow = nrow(test_diff))

df_gr_mean[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_new$gr)), function(y) mean(x[output$cl_new$gr == y]) )))

df_gr_sd[,] <- t(apply(input_data, 2, function(x) 
  sapply(sort(unique(output$cl_new$gr)), function(y) sd(x[output$cl_new$gr == y]) )))

df_gr_cv <- df_gr_mean/df_gr_sd
colnames(df_gr_cv) <- colnames(df_gr_mean) <- colnames(df_gr_sd) <- paste0('gr_', sort(unique(output$cl_new$gr)))
df_gr_cv <- as.data.frame(df_gr_cv)
df_gr_mean <- as.data.frame(df_gr_mean)
df_gr_sd <- as.data.frame(df_gr_sd)
df_gr_mean$id <- df_gr_sd$id <- df_gr_cv$id <- colnames(input_data)

output$gr_input <- list(mean = df_gr_mean, sd = df_gr_sd, cv = df_gr_cv)

# save
save(output, file = sprintf('%s%s_corrPCs_%s_predictCluster%s_PGmethod_HKmetric.RData',  
                            outFold, type_data, type_input, type_cluster))

### plot UMAP and project new samples ###
n_comp_umap <- 2
n_neigh_umap <- clust_res$best_k
min_dist_umap <- 0.01
seed_umap <- 67

custom.settings = umap.defaults
custom.settings$min_dist = min_dist_umap
custom.settings$n_components = n_comp_umap
custom.settings$n_neighbors = n_neigh_umap
custom.settings$random_state <- seed_umap
custom.settings$transform_state <- seed_umap+10

umap_res <- umap::umap(clust_res$input_data, custom.settings)

df_old <- data.frame(component_1=umap_res$layout[,1], 
                     component_2=umap_res$layout[,2],
                     gr = clust_res$cl_best$gr)

df_old$gr <- factor(df_old$gr)
df_old$type <- 'model'
# predict
umap_pred <- predict(umap_res, input_data)

df_new <- data.frame(component_1=umap_pred[,1], component_2=umap_pred[,2], gr = output$cl_new$gr)
df_new$gr <- factor(df_new$gr)
df_new$type <- 'predict'
df <- rbind(df_old, df_new)
df$type <- factor(df$type, levels = c('predict', 'model'))
P <- length(unique(df_old$gr)) # new could be not complete in prediction
gr_color <- pal_d3(palette = 'category20')(P)

tot_pl <- ggplot(df, aes(x = component_1, y = component_2, color = gr, alpha = type))+
  geom_point(size = 0.05)+
  scale_alpha_manual(values = c(1, 0.05))+
  scale_color_manual(values = gr_color)+
  theme_bw()+theme(legend.position = 'right')
width_pl <- 4.5

ggsave(filename = sprintf('%s%s_corrPCs_%s_predictCluster%s_PGmethod_HKmetric_umap.png',  outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'png')
ggsave(filename = sprintf('%s%s_corrPCs_%s_predictCluster%s_PGmethod_HKmetric_umap.pdf',  outFold, type_data, type_input, type_cluster), width = width_pl, height = 4, plot = tot_pl, device = 'pdf')


