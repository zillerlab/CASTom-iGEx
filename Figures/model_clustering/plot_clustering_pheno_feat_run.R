# plot for tissue specific clustering: 
# genome wide genes
# associated pathways (> 5 genes)
# comparison gr vs all: endophenotype association
# treatment response

# script working with R/4.0.2

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="clustering using PG method, show genes/pathways distribution and endophenotype division")
parser$add_argument("--clustFile", type = "character", help = "")
parser$add_argument("--tissue_name", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--tscore_featRelFile", type = "character", help = "")
parser$add_argument("--pathR_featRelFile", type = "character", default  = 'NA', help = "")
parser$add_argument("--pathGO_featRelFile", type = "character",  default  = 'NA', help = "")
parser$add_argument("--endophenoFile", type = "character",  nargs = '*', default  = 'NA', help = "")
parser$add_argument("--endophenoPairwiseFile", type = "character",  nargs = '*', default  = 'NA', help = "")
parser$add_argument("--treatmentResponseFile", type = "character", default  = 'NA', help = "")
parser$add_argument("--treatmentResponsePairwiseFile", type = "character", default  = 'NA', help = "")
parser$add_argument("--geneInfoFile", type = "character",  help = "")
parser$add_argument("--pval_feat", type = "double",  help = "")
parser$add_argument("--pval_pheno", type = "double",  help = "")
parser$add_argument("--pval_treat", type = "double",  help = "")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_cluster_data", type = "character", help = "")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--tscore_info_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")
parser$add_argument("--color_file", type = "character", help = "file with color based on phenotype")
# parser$add_argument("--covDatFile", type = "character", default = 'NA', help = "additional cov to test")
# parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
# parser$add_argument("--pvalresFile", type = "character", default = 'NA', help = "file with pvalue results")
# parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
# parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
# parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
# parser$add_argument("--type_sim", type = "character", default = 'HK', help = "HK or ED")
# parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
# parser$add_argument("--kNN_par", type = "integer", nargs = '*', default = 30, help = "parameter used for PG method")

args <- parser$parse_args()
color_file <- args$color_file
clustFile <- args$clustFile
tissue_name <- args$tissue_name
tscore_featRelFile <- args$tscore_featRelFile
pathR_featRelFile <- args$pathR_featRelFile
pathGO_featRelFile <- args$pathGO_featRelFile
functR <- args$functR
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
geneInfoFile <- args$geneInfoFile
min_genes_path <- args$min_genes_path
type_cluster_data <- args$type_cluster_data
type_cluster <- args$type_cluster
pval_feat <- args$pval_feat
pval_pheno <- args$pval_pheno
tscore_info_file <- args$tscore_info_file
endophenoFile <- args$endophenoFile
endophenoPairwiseFile <- args$endophenoPairwiseFile
treatmentResponseFile <- args$treatmentResponseFile
treatmentResponsePairwiseFile <- args$treatmentResponsePairwiseFile
pval_treat <- args$pval_treat
pheno_name <- args$pheno_name
outFold <- args$outFold

#####################################################################################################################
# color_file <-  '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software//model_clustering/clustering_functions.R'
# type_cluster <- 'Cases'
# pval_feat <- 0.05
# pval_pheno <- 0.001
# pval_treat <- 0.001
# type_cluster_data <- 'tscore'
# tissue_name <- 'Liver'
# min_genes_path <- 5
# pheno_name <- 'CAD'
# geneInfoFile <- sprintf('OUTPUT_GTEx/train_GTEx/%s/200kb/CAD_GWAS_bin5e-2/resPrior_regEval_allchr.txt', tissue_name)
# clustFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_clusterCases_PGmethod_HKmetric.RData', tissue_name)
# tscore_featRelFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscoreOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# pathR_featRelFile<- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_ReactomeOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# pathGO_featRelFile<- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/path_GOOriginal_tscoreClusterCases_featAssociation.RData', tissue_name)
# outFold <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/', tissue_name)
# GOterms_file <- '/psycl/g/mpsziller/lucia/priler_project/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/psycl/g/mpsziller/lucia/priler_project/refData/ReactomePathways.gmt'
# tscore_info_file <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/tscore_info.RData', tissue_name)
# endophenoFile <- c(sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', tissue_name),
#                    sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withoutMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData', tissue_name))
# endophenoPairwiseFile <- c(sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', tissue_name),
#                    sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withoutMedication_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLMpairwise.RData', tissue_name))
# treatmentResponseFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_TreatResponse.txt', tissue_name)
# treatmentResponsePairwiseFile <- sprintf('OUTPUT_GTEx/predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/withMedication_tscore_zscaled_clusterCases_TreatResponse_pairwise.txt', tissue_name)
# #####################################################################################################################

source(functR)

# get clustering structure
res_cl <- get(load(clustFile))
cl <- res_cl$cl_best

##################################################
######## block T-score + path_Reactome ###########
##################################################

tmp <- get(load(tscore_featRelFile))
tscore_input <- tmp$scaleData
keep_gene <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
geneInfo <- read.table(geneInfoFile, h=T, stringsAsFactors = F, sep ='\t')
geneInfo <- geneInfo[geneInfo$dev_geno > 0.01 & geneInfo$test_dev_geno > 0, ]
# remove duplicates (already done for input data)
id_dup <- geneInfo$external_gene_name[duplicated(geneInfo$external_gene_name)]
if(length(id_dup)>0){
  geneInfo <- geneInfo[!geneInfo$external_gene_name %in% id_dup, ]
}
geneInfo_keep <- geneInfo[match(keep_gene, geneInfo$external_gene_name),]
geneInfo_keep$Zstat <- tmp$res_pval[match(keep_gene,tmp$res_pval$external_gene_name),7]
test_feat_tscore <- tmp$test_feat

# remove genes that are from the same locus (use 250kb region) and correlated (LD structure)
# matrix of distance based on start position
chr <- sort(as.numeric(unique(sapply(geneInfo_keep$chrom, function(x) strsplit(x, split  = 'chr')[[1]][2]))))
geneInfo_keep_chr <- list()
for(i in 1:length(chr)){
  print(i)
  tmp_gene <- geneInfo_keep[geneInfo_keep$chrom == paste0('chr', chr[i]),]
  cor_mat <- cor(tscore_input[, match(tmp_gene$external_gene_name, colnames(tscore_input)), drop = F])
  
  if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
    locus_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
    len_gene <- c()
    keep_gene <- c()
    for(j in 1:nrow(locus_list)){
      tmp_sel <-  tmp_gene[locus_list[j,],]
      tmp_sel <- tmp_sel[!tmp_sel$external_gene_name %in% len_gene, ]
      len_gene <- unique(c(len_gene, tmp_sel$external_gene_name))
      keep_gene <- unique(c(keep_gene, tmp_sel$external_gene_name[which.max(abs(tmp_sel$Zstat))]))
    }
     #print(keep_gene)
     geneInfo_keep_chr[[i]] <- tmp_gene[tmp_gene$external_gene_name %in% keep_gene, ]
     dist_mat <- sapply( geneInfo_keep_chr[[i]]$start_position, function(x) abs(x -  geneInfo_keep_chr[[i]]$start_position))
     
     if(any(dist_mat[upper.tri(dist_mat)] < 250000)){
       locus_list <- apply(dist_mat, 1, function(x) x<250000)
       len_gene <- c()
       keep_gene <- c()
       for(j in 1:nrow(locus_list)){
         tmp_sel <-  geneInfo_keep_chr[[i]][locus_list[j,],]
         tmp_sel <- tmp_sel[!tmp_sel$external_gene_name %in% len_gene, ]
         len_gene <- unique(c(len_gene, tmp_sel$external_gene_name))
         keep_gene <- unique(c(keep_gene, tmp_sel$external_gene_name[which.max(abs(tmp_sel$Zstat))]))
       }
       #print(keep_gene)
       geneInfo_keep_chr[[i]] <- geneInfo_keep_chr[[i]][geneInfo_keep_chr[[i]]$external_gene_name %in% keep_gene, ]
     }
     
  }else{
     geneInfo_keep_chr[[i]] <- tmp_gene
  }
}

geneInfo_keep <- do.call(rbind, geneInfo_keep_chr)
# save total feat-gene:
write.table(x = tmp$test_feat, file = sprintf('%stscoreOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)
# save gene info
geneInfo_save <- geneInfo
geneInfo_save$Zstat <- tmp$res_pval[match(geneInfo_save$external_gene_name,tmp$res_pval$external_gene_name),7]
geneInfo_save$relevant_cluster <- F
geneInfo_save$relevant_cluster[geneInfo_save$external_gene_name %in% unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])] <- T
geneInfo_save$relevant_cluster_plot <- F
geneInfo_save$relevant_cluster_plot[geneInfo_save$external_gene_name %in% geneInfo_keep$external_gene_name] <- T
write.table(x = geneInfo_save, file = sprintf('%stscoreOriginal_%sCluster%s_infoGenes.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)


#### pathR ####
tmp <- get(load(pathR_featRelFile))
pathR_input <- tmp$scaleData
feat_rel <- tmp$test_feat
keep_path <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
pathInfo <- tmp$res_pval
pathInfo_keep <- pathInfo[match(keep_path, pathInfo[,1]),]
pathInfo_keep <- pathInfo_keep[pathInfo_keep$ngenes_tscore >= min_genes_path,]
pathInfo_keep$Zstat <- pathInfo_keep[,12]
test_feat_pathR <- tmp$test_feat
# filter out pathway based on shared percentage
gs <- readGmt(reactome_file)
tscore_info <- get(load(tscore_info_file))
# annotate gs with the used genes
for(i in 1:length(gs)){
  gs[[i]]@data <- gs[[i]]@ids[gs[[i]]@ids %in% tscore_info$external_gene_name]
}
pathInfo_save <- pathInfo[, colnames(pathInfo) %in% c('path',"ngenes_tscore","ngenes_path", 
                                                      "mean_dev_geno","sd_dev_geno","mean_test_dev_geno",
                                                      "sd_test_dev_geno","mean_gene_corr","sd_gene_corr")]
pathInfo_save$genes_id <- NA
for(i in 1:nrow(pathInfo_save)){
  id_tmp <- which(sapply(gs, function(x) x@reference == pathInfo_save$path[i]))
  pathInfo_save$genes_id[i] <- paste0(unique(gs[[id_tmp]]@data), collapse = ',')
}
pathInfo_save$Zstat <- pathInfo[match(pathInfo_save$path,pathInfo$path),12]

gs <- gs[which(sapply(gs, function(x) x@reference %in% pathInfo_keep[,1]))]
cor_mat <- cor(pathR_input[, match(pathInfo_keep[,1], colnames(pathR_input)), drop = F])
# remove higly correlated (same gene or genes in LD) and pathways that share high percentage of genes
if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
  path_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
  len_path <- c()
  keep_path <- c()
  for(j in 1:nrow(path_list)){
    tmp_sel <- pathInfo_keep[pathInfo_keep[,1] %in% pathInfo_keep[path_list[j,],1],]
    tmp_sel <-  tmp_sel[!tmp_sel[,1] %in% len_path, ]
    len_path <- unique(c(len_path, tmp_sel[,1]))
    keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
  }
  #print(keep_gene)
  pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
  gs <- gs[which(sapply(gs, function(x) x@reference %in% pathInfo_keep[,1]))]
  perc_mat <- sapply(gs, function(y) sapply(gs, function(x) length(intersect(y@data, x@data))/length(union(y@data, x@data))))
  if(any(perc_mat[upper.tri(perc_mat)] > 0.1)){
    path_list <- apply(perc_mat, 1, function(x) x>0.1)
    len_path <- c()
    keep_path <- c()
    for(j in 1:nrow(path_list)){
      tmp_sel <-  pathInfo_keep[pathInfo_keep[,1] %in% sapply(gs[path_list[j,]], function(x) x@reference),]
      tmp_sel <- tmp_sel[!tmp_sel[,1] %in% len_path, ]
      len_path <- unique(c(len_path, tmp_sel[,1]))
      keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
    }
    #print(keep_path)
    pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
  }
}

# any pathway have a better performance than the single genes reported?
P <- length(unique(cl$gr))
comp <- unique(test_feat_tscore$comp)
impr_path <- c()
for(i in 1:nrow(pathInfo_keep)){
  tmp <- gs[[which(sapply(gs, function(x) x@reference == pathInfo_keep[i,1]))]]
  test_feat_tscore_tmp <- test_feat_tscore[test_feat_tscore$feat %in% tmp@data,]
  test_feat_path_tmp <- test_feat_pathR[test_feat_pathR$feat == pathInfo_keep[i,1], ]
  test_feat_path_tmp <- test_feat_path_tmp[match(comp, test_feat_path_tmp$comp),]
  tmp_gr_tscore <- lapply(comp, function(x) test_feat_tscore_tmp$pval[test_feat_tscore_tmp$comp == x])
  tscore_path_comp <- mapply(function(x, y, z) all(x < y & z <= 0.05), y = tmp_gr_tscore, x = test_feat_path_tmp$pval, z = test_feat_path_tmp$pval_corr)
  if(any(tscore_path_comp)){
    print(i)
    impr_path <- c(impr_path, pathInfo_keep[i,1])
  }
}
pathInfo_keep$impr <- F
pathInfo_keep$impr[pathInfo_keep[,1] %in% impr_path] <- T

# save total feat-path:
write.table(x = feat_rel, file = sprintf('%spath_ReactomeOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

# save path info
pathInfo_save$relevant_cluster <- F
pathInfo_save$relevant_cluster[pathInfo_save$path %in%  unique(feat_rel$feat[feat_rel$pval_corr <= pval_feat])] <- T
pathInfo_save$relevant_cluster_plot <- F
pathInfo_save$relevant_cluster_plot[pathInfo_save$path %in% pathInfo_keep$path] <- T
write.table(x = pathInfo_save, file = sprintf('%spath_ReactomeOriginal_%sCluster%s_infopath.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

### plot ###
pheat_pl_tot(pheno_name = pheno_name, mat_tscore = tscore_input, mat_path = pathR_input, pval_thr_est = 0.05,
             info_feat_tscore = geneInfo_keep, info_feat_path = pathInfo_keep, cl = cl, 
             test_feat_tscore = test_feat_tscore, test_feat_path = test_feat_pathR, 
             cap = 3, res_pl = 250, width_pl = 14 + round(length(unique(cl$gr))*0.25), 
             height_pl = round(7 + (nrow(geneInfo_keep)+nrow(pathInfo_keep))*0.1), 
             outFile = sprintf('%stscoreOriginal_path_ReactomeOriginal_%sCluster%s',outFold,type_cluster_data,type_cluster))

print('tscore-path_Reactome heatmap finished')

############################################
######## block T-score + path_GO ###########
############################################

tmp <- get(load(pathGO_featRelFile))
pathGO_input <- tmp$scaleData
test_feat_pathGO <- tmp$test_feat
print(identical(tmp$res_pval$path_id, colnames(pathGO_input)))
print(identical(tmp$res_pval$path_id, tmp$test_feat$feat))
colnames(pathGO_input) <- tmp$res_pval$path
test_feat_pathGO$feat <- tmp$res_pval$path
  
keep_path <- unique(tmp$test_feat$feat[tmp$test_feat$pval_corr <= pval_feat])
pathInfo <- tmp$res_pval
pathInfo_keep <- pathInfo[match(keep_path, pathInfo[,1]),]
pathInfo_keep <- pathInfo_keep[pathInfo_keep$ngenes_tscore >= min_genes_path + 5,]
pathInfo_keep$Zstat <- pathInfo_keep[,14]
pathInfo_keep[,1] <- pathInfo_keep$path

# filter out pathway based on shared percentage
go <- get(load(GOterms_file))
# annotate go with the used genes
for(i in 1:length(go)){
  go[[i]]$new <- go[[i]]$geneIds[go[[i]]$geneIds %in% tscore_info$external_gene_name]
}

pathInfo_save <- pathInfo[, colnames(pathInfo) %in% c('path_id','path',"ngenes_tscore","ngenes_path", 
                                                      "mean_dev_geno","sd_dev_geno","mean_test_dev_geno",
                                                      "sd_test_dev_geno","mean_gene_corr","sd_gene_corr")]
pathInfo_save$genes_id <- NA
for(i in 1:nrow(pathInfo_save)){
  id_tmp <- which(sapply(go, function(x) x$Term == pathInfo_save$path[i]))
  pathInfo_save$genes_id[i] <- paste0(unique(go[[id_tmp]]$new), collapse = ',')
}
pathInfo_save$Zstat <- pathInfo[match(pathInfo_save$path,pathInfo$path),14]

go <- go[which(sapply(go, function(x) x$Term %in% pathInfo_keep[,1]))]
cor_mat <- cor(pathGO_input[, match(pathInfo_keep[,1], colnames(pathGO_input)), drop = F])
# remove higly correlated (same gene or genes in LD) and pathways that share high percentage of genes
if(any(abs(cor_mat[upper.tri(cor_mat)]) > 0.6)){
  path_list <- apply(cor_mat, 1, function(x) abs(x)>0.6)
  len_path <- c()
  keep_path <- c()
  for(j in 1:nrow(path_list)){
    tmp_sel <- pathInfo_keep[pathInfo_keep[,1] %in% pathInfo_keep[path_list[j,],1],]
    tmp_sel <-  tmp_sel[!tmp_sel[,1] %in% len_path, ]
    len_path <- unique(c(len_path, tmp_sel[,1]))
    keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
  }
  #print(keep_gene)
  pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
  go <- go[which(sapply(go, function(x) x$Term %in% pathInfo_keep[,1]))]
  
  perc_mat <- sapply(go, function(y) sapply(go, function(x) length(intersect(y$new, x$new))/length(union(y$new, x$new))))
  if(any(perc_mat[upper.tri(perc_mat)] > 0.1)){
   path_list <- apply(perc_mat, 1, function(x) x>0.1)
   len_path <- c()
   keep_path <- c()
    for(j in 1:nrow(path_list)){
      tmp_sel <-  pathInfo_keep[pathInfo_keep[,1] %in% sapply(go[path_list[j,]], function(x) x$Term),]
      tmp_sel <- tmp_sel[!tmp_sel[,1] %in% len_path, ]
      len_path <- unique(c(len_path, tmp_sel[,1]))
      keep_path <- unique(c(keep_path, tmp_sel[which.max(abs(tmp_sel$Zstat)), 1]))
   }
   #print(keep_path)
    pathInfo_keep <- pathInfo_keep[pathInfo_keep[,1] %in% keep_path, ]
  }
}

# any pathway have a better performance than the single genes reported?
P <- length(unique(cl$gr))
comp <- unique(test_feat_tscore$comp)
impr_path <- c()
for(i in 1:nrow(pathInfo_keep)){
  tmp <- go[[which(sapply(go, function(x) x$Term == pathInfo_keep[i,1]))]]
  test_feat_tscore_tmp <- test_feat_tscore[test_feat_tscore$feat %in% tmp$new,]
  test_feat_path_tmp <- test_feat_pathGO[test_feat_pathGO$feat == pathInfo_keep[i,1], ]
  test_feat_path_tmp <- test_feat_path_tmp[match(comp, test_feat_path_tmp$comp),]
  tmp_gr_tscore <- lapply(comp, function(x) test_feat_tscore_tmp$pval[test_feat_tscore_tmp$comp == x])
  tscore_path_comp <- mapply(function(x, y, z) all(x < y & z <= 0.05), y = tmp_gr_tscore, x = test_feat_path_tmp$pval, z = test_feat_path_tmp$pval_corr)
  if(any(tscore_path_comp)){
    print(i)
    impr_path <- c(impr_path, pathInfo_keep[i,1])
  }
}
pathInfo_keep$impr <- F
pathInfo_keep$impr[pathInfo_keep[,1] %in% impr_path] <- T

# save total feat-path:
write.table(x = test_feat_pathGO, file = sprintf('%spath_GOOriginal_%sCluster%s_featAssociation.txt', outFold, type_cluster_data, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

# save path info
pathInfo_save$relevant_cluster <- F
pathInfo_save$relevant_cluster[pathInfo_save$path %in%  unique(test_feat_pathGO$feat[test_feat_pathGO$pval_corr <= pval_feat])] <- T
pathInfo_save$relevant_cluster_plot <- F
pathInfo_save$relevant_cluster_plot[pathInfo_save$path %in% pathInfo_keep$path] <- T
write.table(x = pathInfo_save, file = sprintf('%spath_GOOriginal_%sCluster%s_infopath.txt', outFold, type_cluster_data,type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

### plot ###
pheat_pl_tot(pheno_name = pheno_name, mat_tscore = tscore_input, mat_path = pathGO_input, pval_thr_est = 0.05,
             info_feat_tscore = geneInfo_keep, info_feat_path = pathInfo_keep, cl = cl, 
             test_feat_tscore = test_feat_tscore, test_feat_path = test_feat_pathGO, 
             cap = 3, res_pl = 250,  width_pl = round(13 + length(unique(cl$gr))*0.25), 
             height_pl = round(7 + (nrow(geneInfo_keep)+nrow(pathInfo_keep))*0.1), 
             outFile = sprintf('%stscoreOriginal_path_GOOriginal_%sCluster%s', outFold,type_cluster_data,type_cluster))

print('tscore-path_GO heatmap finished')

# pheat_pl_tscore(mat = tscore_input, test_feat = tmp$test_feat, info_feat = geneInfo_keep, cl = cl ,
#                 outFile = sprintf('%stscoreOriginal_%sCluster%s',outFold,type_cluster_data,type_cluster), cap = 3, res_pl = 300, 
#                 width_pl = round(9 + length(unique(cl$gr))*0.25), height_pl = round(5 + nrow(geneInfo_keep)*0.1))
# 
# pheat_pl_path(mat = pathR_input, test_feat = tmp$test_feat, info_feat = pathInfo_keep, cl = cl ,
#                   outFile = sprintf('%spath_ReactomeOriginal_%sCluster%s',outFold,type_cluster_data,type_cluster), cap = 3, res_pl = 300, 
#                   width_pl = round(9 + length(unique(cl$gr))*0.25), height_pl = round(5 + nrow(pathInfo_keep)*0.1))
#   

##############################
#### endophenotype plots #####
##############################

pheno_ann <- read.delim(color_file, header = T, stringsAsFactors = F)
pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4', 'brown', 'brown','chartreuse4', 'brown'), pheno_type = c('ICD9-10_OPCS4', 'Medications', 'Medication', 
                                                                                                                           'Medical_conditions', 'Alcohol', 'Asthma_related_drugs')))
pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'
pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
pheno_ann$color[pheno_ann$pheno_type %in% c('ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')] <- 'grey40'

if(file.exists(endophenoFile[1])){

  res_pheno <- list()
  for(i in 1:length(endophenoFile)){
    tmp <- get(load(endophenoFile[[i]]))
    res_pheno[[i]] <- tmp$bin_reg
    if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
      tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
      tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
      tmp$phenoInfo$pheno_type <- tmp_name
      tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
    }
    res_pheno[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(res_pheno[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
  }
  res_pheno <- do.call(rbind, res_pheno)
  if(length(endophenoFile)>1){
    comp <- unique(res_pheno$comp)
    tmp <- list()
    for(i in 1:length(comp)){
      tmp[[i]] <- res_pheno[res_pheno$comp == comp[i],]
      tmp[[i]]$pval_corr <- p.adjust(tmp[[i]]$pvalue, method = 'BH')
    }
    res_pheno <- do.call(rbind, tmp)
    # save results
    write.table(x = res_pheno, 
                file = sprintf('%s%s_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_combined.txt', outFold, type_cluster_data , type_cluster), 
                col.names = T, row.names = F, sep = '\t', quote = F)
  }

  id_keep <- unique(res_pheno$pheno_id[res_pheno$pvalue <= pval_pheno | res_pheno$pval_corr <= 0.05])
  df_red <- res_pheno[res_pheno$pheno_id %in% id_keep, ]
  
  df_red$new_id <- df_red$Field
  df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
  df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$sign <- 'no'
  df_red$sign[df_red$pval_corr <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$type_pheno!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  pheno_ann_red1 <- pheno_ann[match(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS'], pheno_ann$pheno_type), ]
  pheno_ann_red2 <- pheno_ann[match(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS'], pheno_ann$pheno_type), ]
  # pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
  # pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]
  # pheno_ann_red <- pheno_ann[match(unique(df_red$pheno_type), pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red$comp))
  len_h <- length(unique(df_red$pheno_id))
  # change labels 
  labs_new <- sapply(as.character(unique(df_red$comp)), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  names(labs_new) <- as.character(unique(df_red$comp))
  
  # # cannot work!  
  # pl_tot <-  ggplot(df_red, aes(x = new_id, y = OR_or_Beta, color = pheno_type, shape = sign))+
  #   geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
  #   theme_bw()+ 
  #   ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1)+
  #   facet_grid(type_res~comp, scales = 'free')+
  #   theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(),
  #         axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7), strip.text = element_text(size=7.5))+
  #   scale_shape_manual(values=c(1, 19))+
  #   scale_color_manual(values=pheno_ann_red$color)+
  #   scale_y_continuous(trans='log')+
  #   coord_flip()
  # ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.png', outFold, type_cluster), width = len_w+4, height = len_h*0.2+1.5, plot = pl_tot, device = 'png')
  # ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.pdf', outFold, type_cluster), width = len_w+4, height = len_h*0.2+1.5, plot = pl_tot, device = 'pdf')

  P <- length(unique(df_red$comp))
  gr_color <- pal_d3(palette = 'category20')(P)
  
  if(any(df_red$type_pheno != 'CONTINUOUS')){
    
    pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
      geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
      theme_bw()+ 
      ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
      facet_wrap(comp~.,  nrow = 1, strip.position="top", labeller = labeller(comp = labs_new))+
      theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7, colour = pheno_ann_red1$color),
            strip.text = element_text(size=8, color = 'white', face = 'bold'))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=pheno_ann_red1$color)+
      scale_y_continuous(trans='log', labels = scales::number_format(accuracy = 0.01))+
      coord_flip()
    
    pl_OR <- ggplot_gtable(ggplot_build(pl_OR))
    stripr <- which(grepl('strip-t', pl_OR$layout$name))
    fills <- gr_color
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', pl_OR$grobs[[i]]$grobs[[1]]$childrenOrder))
      pl_OR$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  
  if(any(df_red$type_pheno == 'CONTINUOUS')){
    
    pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
      geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(comp~., nrow = 1, strip.position="top",  labeller = labeller(comp = labs_new))+
      theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red2$color), 
            strip.text = element_text(size=8, color = 'white', face = 'bold'))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=pheno_ann_red2$color)+
      coord_flip()
    ratio_OR_beta <- sum(df_red$type_pheno == 'CONTINUOUS')/sum(df_red$type_pheno != 'CONTINUOUS')
    
    pl_beta <- ggplot_gtable(ggplot_build(pl_beta))
    stripr <- which(grepl('strip-t', pl_beta$layout$name))
    fills <- gr_color
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', pl_beta$grobs[[i]]$grobs[[1]]$childrenOrder))
      pl_beta$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  
  if(any(df_red$type_pheno == 'CONTINUOUS') & any(df_red$type_pheno != 'CONTINUOUS')){
    tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, ratio_OR_beta))
  }else{
    if(any(df_red$type_pheno == 'CONTINUOUS')){
      tot_pl <- pl_beta
    }else{
      tot_pl <- pl_OR
    }
  }
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.png', outFold, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLM_betaOR.pdf', outFold, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'pdf')

}

if(file.exists(endophenoPairwiseFile[1])){

  res_pheno <- list()
  for(i in 1:length(endophenoPairwiseFile)){
    tmp <- get(load(endophenoPairwiseFile[[i]]))
    res_pheno[[i]] <- tmp$bin_reg
    if(!'pheno_type' %in% colnames(tmp$phenoInfo)){
      tmp_name <- sapply(tmp$phenoInfo$Path, function(x) strsplit(x, split = '> ')[[1]][length(strsplit(x, split = '> ')[[1]])])
      tmp_name <- sapply(tmp_name, function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))
      tmp$phenoInfo$pheno_type <- tmp_name
      tmp$phenoInfo$pheno_type[tmp$phenoInfo$pheno_type == 'Summary_Information_(diagnoses)'] <- 'ICD9-10_OPCS4'
    }
    res_pheno[[i]]$pheno_type <- tmp$phenoInfo$pheno_type[match(res_pheno[[i]]$pheno_id,tmp$phenoInfo$pheno_id)]
  }
  res_pheno <- do.call(rbind, res_pheno)
  
  if(length(endophenoPairwiseFile)>1){
    comp <- unique(res_pheno$comp)
    tmp <- list()
    for(i in 1:length(comp)){
      tmp[[i]] <- res_pheno[res_pheno$comp == comp[i],]
      tmp[[i]]$pval_corr <- p.adjust(tmp[[i]]$pvalue, method = 'BH')
    }
    res_pheno <- do.call(rbind, tmp)
    # save results
    write.table(x = res_pheno, 
                file = sprintf('%s%s_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLMpairwise_combined.txt', outFold, type_cluster_data , type_cluster), 
                col.names = T, row.names = F, sep = '\t', quote = F)
  }

  id_keep <- unique(res_pheno$pheno_id[res_pheno$pvalue <= pval_pheno | res_pheno$pval_corr <= 0.05])
  df_red <- res_pheno[res_pheno$pheno_id %in% id_keep, ]
  
  df_red$new_id <- df_red$Field
  df_red$new_id[!is.na(df_red$meaning)] <- paste(df_red$meaning[!is.na(df_red$meaning)], df_red$Field[!is.na(df_red$meaning)], sep = '\n')
  df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$sign <- 'no'
  df_red$sign[df_red$pval_corr <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$type_pheno!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  pheno_ann_red1 <- pheno_ann[match(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS'], pheno_ann$pheno_type), ]
  pheno_ann_red2 <- pheno_ann[match(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS'], pheno_ann$pheno_type), ]
  # pheno_ann_red1 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno != 'CONTINUOUS']), pheno_ann$pheno_type), ]
  # pheno_ann_red2 <- pheno_ann[match(unique(df_red$pheno_type[df_red$type_pheno == 'CONTINUOUS']), pheno_ann$pheno_type), ]
  # pheno_ann_red <- pheno_ann[match(unique(df_red$pheno_type), pheno_ann$pheno_type), ]
  
  len_w <- length(unique(df_red$comp))
  len_h <- length(unique(df_red$pheno_id))

  pl_OR <-  ggplot(subset(df_red, type_pheno != 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
    geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
    theme_bw()+ 
    ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
    facet_wrap(comp~.,  nrow = 1, strip.position="top")+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7, colour = pheno_ann_red1$color),
          strip.text = element_text(size=8, color = 'black'))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=pheno_ann_red1$color)+
    scale_y_continuous(trans='log', labels = scales::number_format(accuracy = 0.01))+
    coord_flip()
  
  pl_beta <-  ggplot(subset(df_red, type_pheno == 'CONTINUOUS'), aes(x = new_id, y = OR_or_Beta, shape = sign))+
    geom_point()+geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(0.05))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+ geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(comp~., nrow = 1, strip.position="top")+
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7,  colour = pheno_ann_red2$color), 
          strip.text = element_text(size=8, color = 'black'))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=pheno_ann_red2$color)+
    coord_flip()
  ratio_OR_beta <- sum(df_red$type_pheno == 'CONTINUOUS')/sum(df_red$type_pheno != 'CONTINUOUS')
  
  if(any(df_red$type_pheno == 'CONTINUOUS') & any(df_red$type_pheno != 'CONTINUOUS')){
    tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1, heights = c(1, ratio_OR_beta))
  }else{
    if(any(df_red$type_pheno == 'CONTINUOUS')){
      tot_pl <- pl_beta
    }else{
      tot_pl <- pl_OR
    }
  }
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.png', outFold, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_phenoAssociation_GLMpairwise_betaOR.pdf', outFold, type_cluster), width = len_w+3, height = len_h*0.2+2, plot = tot_pl, device = 'pdf')
  
}

###################################
#### treatment response plots #####
###################################


if(file.exists(treatmentResponseFile)){
  
  res_treat <- fread(treatmentResponseFile, h=T, stringsAsFactors = F, data.table = F)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= pval_treat | res_treat$pvalue_corr_diff <= 0.05])
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning)
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_all')[[1]][1])
  df_red$gr <- factor(df_red$gr, levels = unique(df_red$gr))
  
  len_w <- round(nrow(df_red)/length(unique(df_red$gr)))
  len_h <- 4
    
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_red$gr)))
  
  if(any(df_red$pheno_type == 'CONTINUOUS')){
    pl_beta <-  ggplot(subset(df_red, pheno_type == 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~., nrow = 1, strip.position="right", scales = 'free_x')+
      theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), 
            plot.title = element_text(size=9), axis.title.x = element_blank(),  axis.title.y = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=8))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+guides(shape=FALSE)
     # coord_flip()
  }
  if(any(df_red$pheno_type!= 'CONTINUOUS')){
    pl_OR <-  ggplot(subset(df_red, pheno_type != 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~.,  nrow = 1, strip.position="right", scales = 'free_x')+
      theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), axis.title.y = element_text(size = 7),
            plot.title = element_text(size=9), axis.title.x = element_blank(),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7),
            strip.text = element_text(size=8))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+
      scale_y_continuous(trans='log')+guides(shape=FALSE)
      # coord_flip()
    ratio_OR_beta <- sum(df_red$pheno_type == 'CONTINUOUS')/sum(df_red$pheno_type != 'CONTINUOUS')
  }
    
    if(any(df_red$pheno_type == 'CONTINUOUS') & any(df_red$pheno_type != 'CONTINUOUS')){
      tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'h', nrow = 1,  common.legend = T, widths = c(1, ratio_OR_beta*0.6))
    }else{
      if(any(df_red$pheno_type == 'CONTINUOUS')){
        tot_pl <- pl_beta
      }else{
        tot_pl <- pl_OR
      }
    }
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse.png', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse.pdf', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}


if(file.exists(treatmentResponsePairwiseFile)){
  
  res_treat <- fread(treatmentResponsePairwiseFile, h=T, stringsAsFactors = F, data.table = F)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= pval_treat | res_treat$pvalue_corr_diff <= 0.05])
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning)
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr1 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  # df_red$gr1 <- factor(df_red$gr1, levels = unique(df_red$gr1))
  df_red$gr2 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  # df_red$gr2 <- factor(df_red$gr2, levels = unique(df_red$gr2))
  
  tot_gr <- sort(unique(c(as.character(df_red$gr1),  as.character(df_red$gr2))))
  
  df_tot_gr <- list()
  for(i in 1:length(tot_gr)){
    
    tmp <- df_red[df_red$gr1 == tot_gr[i] | df_red$gr2 == tot_gr[i], ]
    feat_name <- unique(tmp$comb_name)
    df_tot_gr[[i]] <- data.frame(comb_name =c(),gr =c(), sign =c(),pheno_class =c(),new_id =c(),treat_meaning  =c(),
                         gr_ORorBeta =c(),gr_CI_low =c(), gr_CI_up =c())
    for(j in 1:length(feat_name)){
      df_new <- data.frame(comb_name = feat_name[j], gr = tot_gr[i], sign = any(tmp$sign[tmp$comb_name == feat_name[j]] == 'yes'), 
                           pheno_class = tmp$pheno_class[tmp$comb_name == feat_name[j]][1], 
                           new_id = tmp$new_id[tmp$comb_name == feat_name[j]][1], 
                           treat_meaning = tmp$treat_meaning[tmp$comb_name == feat_name[j]][1], 
                           pheno_type = tmp$pheno_type[tmp$comb_name == feat_name[j]][1])
      name_gr <- ifelse(any(tot_gr[i] %in% tmp$gr1[tmp$comb_name == feat_name[j]]), 'gr1', 'gr2')
      df_new$gr_ORorBeta <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_ORorBeta')][1]
      df_new$gr_CI_low <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_low')][1]
      df_new$gr_CI_up <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_up')][1]
      df_tot_gr[[i]] <- rbind(df_tot_gr[[i]], df_new)
    }
  }
  df_tot_gr <- do.call(rbind, df_tot_gr)
  df_tot_gr$gr <- factor(df_tot_gr$gr, levels = unique(df_tot_gr$gr))
  
  len_w <- round(nrow(df_tot_gr)/length(unique(df_tot_gr$gr)))
  len_h <- 4
  
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_tot_gr$gr)))
  
  if(any(df_tot_gr$pheno_type == 'CONTINUOUS')){
    pl_beta <-  ggplot(subset(df_tot_gr, pheno_type == 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~., nrow = 1, strip.position="right", scales = 'free_x')+
      theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), 
            plot.title = element_text(size=9), axis.title.x = element_blank(),  axis.title.y = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=8))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+guides(shape=FALSE)
    # coord_flip()
  }
  
  if(any(df_tot_gr$pheno_type!= 'CONTINUOUS')){
    pl_OR <-  ggplot(subset(df_tot_gr, pheno_type != 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~.,  nrow = 1, strip.position="right", scales = 'free_x')+
      theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), axis.title.y = element_text(size = 7),
            plot.title = element_text(size=9), axis.title.x = element_blank(),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7),
            strip.text = element_text(size=8))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+
      scale_y_continuous(trans='log')+guides(shape=FALSE)
    # coord_flip()
    
    ratio_OR_beta <- sum(df_red$pheno_type == 'CONTINUOUS')/sum(df_red$pheno_type != 'CONTINUOUS')
  }
  
  if(any(df_red$pheno_type == 'CONTINUOUS') & any(df_red$pheno_type != 'CONTINUOUS')){
    tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'h', nrow = 1,  common.legend = T, widths = c(1, ratio_OR_beta*0.6))
  }else{
    if(any(df_red$pheno_type == 'CONTINUOUS')){
      tot_pl <- pl_beta
    }else{
      tot_pl <- pl_OR
    }
  }
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise.png', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise.pdf', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}

########################################################################################
###### Treatment response plot: specific for statin and blood pressure medicine ########
########################################################################################

if(file.exists(treatmentResponseFile) & pheno_name == 'Asthma'){
  
  res_treat <- fread(treatmentResponseFile, h=T, stringsAsFactors = F, data.table = F)
  # consider only medicine specific for CAD
  p1 <- res_treat[res_treat$pheno_Field %in% c('Forced expiratory volume in 1-second (FEV1), predicted', 
                                               'C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage', 'Basophill percentage', 'Eosinophill percentage', 
                                               'Monocyte percentage', 'Neutrophill percentage') & 
                    res_treat$treat_meaning %in% c('H1 antihistamine', 'Muscarinic antagonists (SAMA)', 'Muscarinic antagonists (LAMA)', 
                    'Methylxanthines', 'Corticosteroids', 'Degranulation inhibitors', 'Leukotrienes', 'non-selective beta-2-agonists', 
                     'Selective beta-2-agonists (SABA)', 'Selective beta-2-agonists (LABA)'),]
  name_p1 <- unique(p1$treat_meaning)
  
  p2 <- res_treat[res_treat$pheno_Field %in% c('C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage',
                                               'Platelet count', 'Platelet crit', 'Platelet distribution width') & 
                    res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin'),]
  res_treat <- rbind(p1, p2)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 1])
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
                                                                  name_p1))
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_all')[[1]][1])
  df_red$gr <- factor(df_red$gr, levels = unique(df_red$gr))
  
  len_w <- 10
  len_h <- 12 + length(unique(df_red$gr))*0.1
  
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_red$gr)))
  
  pl_beta_p1 <-  ggplot(subset(df_red, treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
      theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
            plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
            axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
    
    pl_beta_p2 <-  ggplot(subset(df_red, !treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin')), 
                          aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
      geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
      theme_bw()+ 
      ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
      facet_wrap(treat_meaning~., nrow = 3, strip.position="top", scales = 'free_x')+
      theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
            plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
            axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
 
  tot_pl <- ggarrange(plotlist = list(pl_beta_p2, pl_beta_p1), align = 'hv', nrow = 2, common.legend = T, heights=c(1, 0.35))

  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_%smedAll.png', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_%smedAll.pdf', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}

if(file.exists(treatmentResponseFile) & pheno_name == 'CAD'){
  
  res_treat <- fread(treatmentResponseFile, h=T, stringsAsFactors = F, data.table = F)
  # consider only medicine specific for CAD
  p1 <- res_treat[res_treat$pheno_Field %in% c('Apolipoprotein A', 'Apolipoprotein B', 'Body fat percentage', 
                                               'Body mass index (BMI)', 'Cholesterol',  'Diastolic blood pressure, manual reading', 
                                               'Systolic blood pressure, manual reading', 'HDL cholesterol', 'LDL direct',  
                                               'Triglycerides', 'Pulse wave Arterial Stiffness index', 'Pulse rate', 'Lymphocyte percentage','Lymphocyte count', 'C-reactive protein') & 
                    res_treat$treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication'),]
  
  p2 <- res_treat[res_treat$pheno_Field %in% c('C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage','Lymphocyte count',
                                               'Platelet count', 'Platelet crit', 'Platelet distribution width') & 
                    res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin'),]
  res_treat <- rbind(p1, p2)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 1])
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
                                                                  'Blood pressure medication', 'Cholesterol lowering medication'))
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_all')[[1]][1])
  df_red$gr <- factor(df_red$gr, levels = unique(df_red$gr))
  
  len_w <- 7
  len_h <- 9 + length(unique(df_red$gr))*0.1
  
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_red$gr)))
  
  pl_beta_p1 <-  ggplot(subset(df_red, treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 1, strip.position="top")+
    theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  pl_beta_p2 <-  ggplot(subset(df_red, !treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 1, strip.position="top")+
    theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  tot_pl <- ggarrange(plotlist = list(pl_beta_p1, pl_beta_p2), align = 'hv', nrow = 2, common.legend = T, heights = c(1, 0.7))
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_%smedAll.png', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_%smedAll.pdf', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}

if(file.exists(treatmentResponsePairwiseFile) & pheno_name == 'Asthma'){
  
  res_treat <- fread(treatmentResponsePairwiseFile, h=T, stringsAsFactors = F, data.table = F)
  # consider only medicine specific for CAD
  p1 <- res_treat[res_treat$pheno_Field %in% c('Forced expiratory volume in 1-second (FEV1), predicted', 
                                               'C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage', 'Basophill percentage', 'Eosinophill percentage', 
                                               'Monocyte percentage', 'Neutrophill percentage') & 
                    res_treat$treat_meaning %in% c('H1 antihistamine', 'Muscarinic antagonists (SAMA)', 'Muscarinic antagonists (LAMA)', 
                                                   'Methylxanthines', 'Corticosteroids', 'Degranulation inhibitors', 'Leukotrienes', 'non-selective beta-2-agonists', 
                                                   'Selective beta-2-agonists (SABA)', 'Selective beta-2-agonists (LABA)'),]
  name_p1 <- unique(p1$treat_meaning)
  
  p2 <- res_treat[res_treat$pheno_Field %in% c('C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage',
                                               'Platelet count', 'Platelet crit', 'Platelet distribution width') & 
                    res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin'),]
  res_treat <- rbind(p1, p2)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 1])
  # id_keep_pheno <- unique(res_treat$pheno_id[res_treat$pvalue_diff <= 0.001])
  # id_keep_treat <- unique(res_treat$treat_id[res_treat$pvalue_diff <= 0.001])
  # df_red <- res_treat[res_treat$pheno_id %in% id_keep_pheno & res_treat$treat_id %in% id_keep_treat, ]
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
                                                                  name_p1))
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr1 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  # df_red$gr1 <- factor(df_red$gr1, levels = unique(df_red$gr1))
  df_red$gr2 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  # df_red$gr2 <- factor(df_red$gr2, levels = unique(df_red$gr2))
  
  tot_gr <- sort(unique(c(as.character(df_red$gr1),  as.character(df_red$gr2))))
  
  df_tot_gr <- list()
  for(i in 1:length(tot_gr)){
    
    tmp <- df_red[df_red$gr1 == tot_gr[i] | df_red$gr2 == tot_gr[i], ]
    feat_name <- unique(tmp$comb_name)
    df_tot_gr[[i]] <- data.frame(comb_name =c(),gr =c(), sign =c(),pheno_class =c(),new_id =c(),treat_meaning  =c(),
                                 gr_ORorBeta =c(),gr_CI_low =c(), gr_CI_up =c())
    for(j in 1:length(feat_name)){
      df_new <- data.frame(comb_name = feat_name[j], gr = tot_gr[i], sign = any(tmp$sign[tmp$comb_name == feat_name[j]] == 'yes'), 
                           pheno_class = tmp$pheno_class[tmp$comb_name == feat_name[j]][1], 
                           new_id = tmp$new_id[tmp$comb_name == feat_name[j]][1], 
                           treat_meaning = tmp$treat_meaning[tmp$comb_name == feat_name[j]][1], 
                           pheno_type = tmp$pheno_type[tmp$comb_name == feat_name[j]][1])
      name_gr <- ifelse(any(tot_gr[i] %in% tmp$gr1[tmp$comb_name == feat_name[j]]), 'gr1', 'gr2')
      df_new$gr_ORorBeta <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_ORorBeta')][1]
      df_new$gr_CI_low <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_low')][1]
      df_new$gr_CI_up <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_up')][1]
      df_tot_gr[[i]] <- rbind(df_tot_gr[[i]], df_new)
    }
  }
  df_tot_gr <- do.call(rbind, df_tot_gr)
  df_tot_gr$gr <- factor(df_tot_gr$gr, levels = unique(df_tot_gr$gr))
  
  len_w <- 10
  len_h <- 12 + length(unique(df_red$gr))*0.1
  
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_tot_gr$gr)))
  
  pl_beta_p1 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 1, strip.position="top", scales = 'free_x')+
    theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  pl_beta_p2 <-  ggplot(subset(df_tot_gr, !treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 3, strip.position="top", scales = 'free_x')+
    theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  tot_pl <- ggarrange(plotlist = list(pl_beta_p2, pl_beta_p1), align = 'hv', nrow = 2, common.legend = T, heights=c(1, 0.35))
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll.png', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll.pdf', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}




if(file.exists(treatmentResponsePairwiseFile) & pheno_name == 'CAD'){
  
  res_treat <- fread(treatmentResponsePairwiseFile, h=T, stringsAsFactors = F, data.table = F)
  # consider only medicine specific for CAD
  p1 <- res_treat[res_treat$pheno_Field %in% c('Apolipoprotein A', 'Apolipoprotein B', 'Body fat percentage', 
                                               'Body mass index (BMI)', 'Cholesterol',  'Diastolic blood pressure, manual reading', 
                                               'Systolic blood pressure, manual reading', 'HDL cholesterol', 'LDL direct',  
                                               'Triglycerides', 'Pulse wave Arterial Stiffness index', 'Pulse rate', 'Lymphocyte percentage','Lymphocyte count', 'C-reactive protein') & 
                    res_treat$treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication'),]
  
  p2 <- res_treat[res_treat$pheno_Field %in% c('C-reactive protein', 'White blood cell (leukocyte) count','Lymphocyte percentage','Lymphocyte count',
                                               'Platelet count', 'Platelet crit', 'Platelet distribution width') & 
                    res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin'),]
  res_treat <- rbind(p1, p2)
  res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
  id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 1])
  # id_keep_pheno <- unique(res_treat$pheno_id[res_treat$pvalue_diff <= 0.001])
  # id_keep_treat <- unique(res_treat$treat_id[res_treat$pvalue_diff <= 0.001])
  # df_red <- res_treat[res_treat$pheno_id %in% id_keep_pheno & res_treat$treat_id %in% id_keep_treat, ]
  df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
  
  df_red$new_id <- df_red$pheno_Field
  df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
  # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
  df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
  df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
  df_red$type_res <- 'beta'
  df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
  df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
  df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
                                                                  'Blood pressure medication', 'Cholesterol lowering medication'))
  df_red$sign <- 'no'
  df_red$sign[df_red$pvalue_corr_diff <= 0.05] <- 'yes'
  df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
  df_red$gr1 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
  # df_red$gr1 <- factor(df_red$gr1, levels = unique(df_red$gr1))
  df_red$gr2 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
  # df_red$gr2 <- factor(df_red$gr2, levels = unique(df_red$gr2))
  
  tot_gr <- sort(unique(c(as.character(df_red$gr1),  as.character(df_red$gr2))))
  
  df_tot_gr <- list()
  for(i in 1:length(tot_gr)){
    
    tmp <- df_red[df_red$gr1 == tot_gr[i] | df_red$gr2 == tot_gr[i], ]
    feat_name <- unique(tmp$comb_name)
    df_tot_gr[[i]] <- data.frame(comb_name =c(),gr =c(), sign =c(),pheno_class =c(),new_id =c(),treat_meaning  =c(),
                                 gr_ORorBeta =c(),gr_CI_low =c(), gr_CI_up =c())
    for(j in 1:length(feat_name)){
      df_new <- data.frame(comb_name = feat_name[j], gr = tot_gr[i], sign = any(tmp$sign[tmp$comb_name == feat_name[j]] == 'yes'), 
                           pheno_class = tmp$pheno_class[tmp$comb_name == feat_name[j]][1], 
                           new_id = tmp$new_id[tmp$comb_name == feat_name[j]][1], 
                           treat_meaning = tmp$treat_meaning[tmp$comb_name == feat_name[j]][1], 
                           pheno_type = tmp$pheno_type[tmp$comb_name == feat_name[j]][1])
      name_gr <- ifelse(any(tot_gr[i] %in% tmp$gr1[tmp$comb_name == feat_name[j]]), 'gr1', 'gr2')
      df_new$gr_ORorBeta <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_ORorBeta')][1]
      df_new$gr_CI_low <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_low')][1]
      df_new$gr_CI_up <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_up')][1]
      df_tot_gr[[i]] <- rbind(df_tot_gr[[i]], df_new)
    }
  }
  df_tot_gr <- do.call(rbind, df_tot_gr)
  df_tot_gr$gr <- factor(df_tot_gr$gr, levels = unique(df_tot_gr$gr))
  
  len_w <- 7
  len_h <- 9 + length(unique(df_tot_gr$gr))*0.1
  
  gr_color <- pal_d3(palette = 'category20')(length(unique(df_tot_gr$gr)))
  
  pl_beta_p1 <-  ggplot(subset(df_tot_gr, treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 1, strip.position="top")+
    theme(legend.position = 'none', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  pl_beta_p2 <-  ggplot(subset(df_tot_gr, !treat_meaning %in% c('Blood pressure medication', 'Cholesterol lowering medication')), 
                        aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
    geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
    theme_bw()+ 
    ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
    facet_wrap(treat_meaning~., nrow = 1, strip.position="top")+
    theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 9), 
          plot.title = element_text(size=9), axis.title.y = element_blank(),  axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9))+
    scale_shape_manual(values=c(1, 19))+
    scale_color_manual(values=gr_color)+guides(shape=FALSE)+
    coord_flip()
  
  tot_pl <- ggarrange(plotlist = list(pl_beta_p1, pl_beta_p2), align = 'hv', nrow = 2, common.legend = T, heights = c(1, 0.7))
  
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll.png', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'png')
  ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_pairwise_%smedAll.pdf', outFold, type_cluster, pheno_name), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
  
}



# if(file.exists(treatmentResponseFile)){
#   
#   pval_treat_CAD <- 0.001
#   
#   pheno_ann <- read.delim(color_file, header = T, stringsAsFactors = F)
#   pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4'), pheno_type = c('ICD9-10_OPCS4', 'Medications')))
#   pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'
#   pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
#   
#   res_treat <- fread(treatmentResponseFile, h=T, stringsAsFactors = F, data.table = F)
#   # consider only medicine specific for CAD
#   res_treat <- res_treat[res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
#                                                         'Blood pressure medication', 'Cholesterol lowering medication'),]
#   res_treat <- res_treat[res_treat$pheno_type == 'CONTINUOUS', ]
#   res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
#   id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 0.01])
#   # id_keep_pheno <- unique(res_treat$pheno_id[res_treat$pvalue_diff <= 0.001])
#   # id_keep_treat <- unique(res_treat$treat_id[res_treat$pvalue_diff <= 0.001])
#   # df_red <- res_treat[res_treat$pheno_id %in% id_keep_pheno & res_treat$treat_id %in% id_keep_treat, ]
#   df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
#   
#   df_red$new_id <- df_red$pheno_Field
#   df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
#   # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
#   df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
#   df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
#   df_red$type_res <- 'beta'
#   df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
#   df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
#   df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
#                                                                   'Blood pressure medication', 'Cholesterol lowering medication'))
#   df_red$sign <- 'no'
#   df_red$sign[df_red$pvalue_diff <= pval_treat_CAD] <- 'yes'
#   df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
#   df_red$gr <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_all')[[1]][1])
#   df_red$gr <- factor(df_red$gr, levels = unique(df_red$gr))
#   
#   len_w <- round(max(table(df_red$treat_meaning))/length(unique(df_red$gr))*1.3)
#   len_h <- 6
#   
#   gr_color <- pal_d3(palette = 'category20')(length(unique(df_red$gr)))
#   
#   if(any(df_red$pheno_type == 'CONTINUOUS')){
#     
#     pl_beta <-  ggplot(subset(df_red, pheno_type == 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
#       geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
#       theme_bw()+ 
#       ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
#       facet_wrap(treat_meaning~., nrow = 2, strip.position="top", scales = 'free_x')+
#       theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), 
#             plot.title = element_text(size=9), axis.title.x = element_blank(),  axis.title.y = element_text(size = 7),
#             axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=8))+
#       scale_shape_manual(values=c(1, 19))+
#       scale_color_manual(values=gr_color)+guides(shape=FALSE)
#     # coord_flip()
#   }
#   
#   if(any(df_red$pheno_type!= 'CONTINUOUS')){
#     pl_OR <-  ggplot(subset(df_red, pheno_type != 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
#       geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
#       theme_bw()+ 
#       ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
#       facet_wrap(treat_meaning~.,  nrow = 2, strip.position="top", scales = 'free_x')+
#       theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), axis.title.y = element_text(size = 7),
#             plot.title = element_text(size=9), axis.title.x = element_blank(),
#             axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7),
#             strip.text = element_text(size=8))+
#       scale_shape_manual(values=c(1, 19))+
#       scale_color_manual(values=gr_color)+
#       scale_y_continuous(trans='log')+guides(shape=FALSE)
#     # coord_flip()
#     ratio_OR_beta <- sum(df_red$type_pheno == 'CONTINUOUS')/sum(df_red$type_pheno != 'CONTINUOUS')
#   }
#   
#   if(any(df_red$pheno_type == 'CONTINUOUS') & any(df_red$pheno_type != 'CONTINUOUS')){
#     tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1,  common.legend = T)
#     len_h <- len_h+2 
#   }else{
#     if(any(df_red$pheno_type == 'CONTINUOUS')){
#       tot_pl <- pl_beta
#     }else{
#       tot_pl <- pl_OR
#     }
#   }
#   
#   
#   ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_CADmed.png', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'png')
#   ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_CADmed.pdf', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
#   
# }
# 
# 
# 
# if(file.exists(treatmentResponsePairwiseFile)){
#   
#   pval_treat_CAD <- 0.001
#   
#   pheno_ann <- read.delim(color_file, header = T, stringsAsFactors = F)
#   pheno_ann <- rbind(pheno_ann, data.frame(color = c('grey40', 'chocolate4'), pheno_type = c('ICD9-10_OPCS4', 'Medications')))
#   pheno_ann$color[pheno_ann$pheno_type == 'Family_history'] <- 'orange3'
#   pheno_ann$color[pheno_ann$pheno_type == 'Smoking'] <- 'darkgreen'
#   
#   res_treat <- fread(treatmentResponsePairwiseFile, h=T, stringsAsFactors = F, data.table = F)
#   res_treat <- res_treat[res_treat$treat_meaning %in% c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
#                                                         'Blood pressure medication', 'Cholesterol lowering medication'),]
#   res_treat <- res_treat[res_treat$pheno_type == 'CONTINUOUS', ]
#   res_treat$comb_name <- paste0(res_treat$pheno_id, '_and_', res_treat$treat_id)
#   id_keep <-  unique(res_treat$comb_name[res_treat$pvalue_diff <= 0.01])
#   # id_keep_pheno <- unique(res_treat$pheno_id[res_treat$pvalue_diff <= 0.001])
#   # id_keep_treat <- unique(res_treat$treat_id[res_treat$pvalue_diff <= 0.001])
#   # df_red <- res_treat[res_treat$pheno_id %in% id_keep_pheno & res_treat$treat_id %in% id_keep_treat, ]
#   df_red <- res_treat[res_treat$comb_name %in% id_keep, ]
#   
#   df_red$new_id <- df_red$pheno_Field
#   df_red$new_id[!is.na(df_red$pheno_meaning)] <- paste(df_red$pheno_meaning[!is.na(df_red$pheno_meaning)], df_red$pheno_Field[!is.na(df_red$pheno_meaning)], sep = '\n')
#   # df_red$comp <- factor(df_red$comp, levels = unique(df_red$comp))
#   df_red$new_id <- factor(df_red$new_id, levels = unique(df_red$new_id))
#   df_red$pheno_type <- factor(df_red$pheno_type, levels = unique(df_red$pheno_type))
#   df_red$type_res <- 'beta'
#   df_red$type_res[df_red$pheno_type!= 'CONTINUOUS'] <- 'OR'
#   df_red$type_res <- factor(df_red$type_res, levels = c('OR', 'beta'))
#   df_red$treat_meaning <- factor(df_red$treat_meaning, levels = c('Ibuprofen (e.g. Nurofen)', 'Paracetamol', 'Aspirin', 
#                                                                   'Blood pressure medication', 'Cholesterol lowering medication'))
#   df_red$sign <- 'no'
#   df_red$sign[df_red$pvalue_diff <= pval_treat_CAD] <- 'yes'
#   df_red$sign <- factor(df_red$sign, levels = c('no', 'yes'))
#   df_red$gr1 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][1])
#   # df_red$gr1 <- factor(df_red$gr1, levels = unique(df_red$gr1))
#   df_red$gr2 <- sapply(df_red$comp, function(x) strsplit(x, split = '_vs_')[[1]][2])
#   # df_red$gr2 <- factor(df_red$gr2, levels = unique(df_red$gr2))
#   
#   tot_gr <- sort(unique(c(as.character(df_red$gr1),  as.character(df_red$gr2))))
#   
#   df_tot_gr <- list()
#   for(i in 1:length(tot_gr)){
#     
#     tmp <- df_red[df_red$gr1 == tot_gr[i] | df_red$gr2 == tot_gr[i], ]
#     feat_name <- unique(tmp$comb_name)
#     df_tot_gr[[i]] <- data.frame(comb_name =c(),gr =c(), sign =c(),pheno_class =c(),new_id =c(),treat_meaning  =c(),
#                                  gr_ORorBeta =c(),gr_CI_low =c(), gr_CI_up =c())
#     for(j in 1:length(feat_name)){
#       df_new <- data.frame(comb_name = feat_name[j], gr = tot_gr[i], sign = any(tmp$sign[tmp$comb_name == feat_name[j]] == 'yes'), 
#                            pheno_class = tmp$pheno_class[tmp$comb_name == feat_name[j]][1], 
#                            new_id = tmp$new_id[tmp$comb_name == feat_name[j]][1], 
#                            treat_meaning = tmp$treat_meaning[tmp$comb_name == feat_name[j]][1], 
#                            pheno_type = tmp$pheno_type[tmp$comb_name == feat_name[j]][1])
#       name_gr <- ifelse(any(tot_gr[i] %in% tmp$gr1), 'gr1', 'gr2')
#       df_new$gr_ORorBeta <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_ORorBeta')][1]
#       df_new$gr_CI_low <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_low')][1]
#       df_new$gr_CI_up <- tmp[tmp$comb_name == feat_name[j], paste0(name_gr, '_CI_up')][1]
#       df_tot_gr[[i]] <- rbind(df_tot_gr[[i]], df_new)
#     }
#   }
#   df_tot_gr <- do.call(rbind, df_tot_gr)
#   df_tot_gr$gr <- factor(df_tot_gr$gr, levels = unique(df_tot_gr$gr))
#   
#   len_w <- round(max(table(df_tot_gr$treat_meaning))/length(unique(df_tot_gr$gr))*1.3)
#   len_h <- 6
#   
#   gr_color <- pal_d3(palette = 'category20')(length(unique(df_tot_gr$gr)))
#   
#   if(any(df_tot_gr$pheno_type == 'CONTINUOUS')){
#     
#     pl_beta <-  ggplot(subset(df_tot_gr, pheno_type == 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
#       geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
#       theme_bw()+ 
#       ylab('Adjusted Beta (95% CI)')+geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
#       facet_wrap(treat_meaning~., nrow = 2, strip.position="top", scales = 'free_x')+
#       theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), 
#             plot.title = element_text(size=9), axis.title.x = element_blank(),  axis.title.y = element_text(size = 7),
#             axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=8))+
#       scale_shape_manual(values=c(1, 19))+
#       scale_color_manual(values=gr_color)+guides(shape=FALSE)
#     # coord_flip()
#   }
#   
#   if(any(df_tot_gr$pheno_type!= 'CONTINUOUS')){
#     pl_OR <-  ggplot(subset(df_tot_gr, pheno_type != 'CONTINUOUS'), aes(x = new_id, y = gr_ORorBeta, group = gr, color = gr, shape = sign))+
#       geom_point(position=position_dodge(0.5))+geom_errorbar(aes(ymin=gr_CI_low, ymax=gr_CI_up), width=.2, position=position_dodge(0.5))+
#       theme_bw()+ 
#       ylab('Adjusted OR (95% CI)')+ geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40')+
#       facet_wrap(treat_meaning~.,  nrow = 2, strip.position="top", scales = 'free_x')+
#       theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 7), axis.title.y = element_text(size = 7),
#             plot.title = element_text(size=9), axis.title.x = element_blank(),
#             axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 7),
#             strip.text = element_text(size=8))+
#       scale_shape_manual(values=c(1, 19))+
#       scale_color_manual(values=gr_color)+
#       scale_y_continuous(trans='log')+guides(shape=FALSE)
#     # coord_flip()
#     ratio_OR_beta <- sum(df_tot_gr$type_pheno == 'CONTINUOUS')/sum(df_tot_gr$type_pheno != 'CONTINUOUS')
#   }
#   
#   if(any(df_tot_gr$pheno_type == 'CONTINUOUS') & any(df_tot_gr$pheno_type != 'CONTINUOUS')){
#     tot_pl <- ggarrange(plotlist = list(pl_OR, pl_beta), align = 'v', ncol = 1,  common.legend = T)
#     len_h <- len_h+2 
#   }else{
#     if(any(df_tot_gr$pheno_type == 'CONTINUOUS')){
#       tot_pl <- pl_beta
#     }else{
#       tot_pl <- pl_OR
#     }
#   }
#   
#   
#   ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_CADmed_pairwise.png', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'png')
#   ggsave(filename = sprintf('%stscore_zscaled_cluster%s_PGmethod_HKmetric_treatmentResponse_CADmed_pairwise.pdf', outFold, type_cluster), width = len_w, height = len_h, plot = tot_pl, device = 'pdf')
#   
# }
# 






