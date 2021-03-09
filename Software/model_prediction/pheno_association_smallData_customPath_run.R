#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# save beta, se, zcores, phenotype already processed, use custom pathway structure

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))


parser <- ArgumentParser(description="Gene and Pathway (custom) association analysis")
parser$add_argument("--pathwayStructure_file", type = "character", help = "custom pathway anntation")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--thr_reliableGenes", type = "double", nargs = '*', default = c(0.01, 0), help = "threshold for reliable genes: dev_geno_tot and test_dev_geno")
parser$add_argument("--inputFold", type = "character", help = "Folde with results from pathway analysis")
parser$add_argument("--covDat_file", type = "character", nargs = '*', help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases), associated to phenoDat_file")
parser$add_argument("--phenoDat_file", type = "character", nargs = '*', help = "file(s) with individual_ID to match and phenotype to test association, associated to covDat_file")
parser$add_argument("--names_file", type = "character", nargs = '*', help = "for each couple of covDat/phenoDat file, associated name")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical",default = F,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--ncores", type = "integer", default = 0, help = "number of cores for the parallelization (over genes, if zero not parallelized)")
parser$add_argument("--geneAnn_file", type = "character", help = "file with gene info from train, to be filtered")
parser$add_argument("--functR", type = "character", help = "Rscript with functions to be used")
parser$add_argument("--geneSetName", type = "character", help = "name pathway custom")
parser$add_argument("--abs_tscore", type = "logical", default = F, help = "if true also the pathway using absolute values of tscore is computed")
parser$add_argument("--not_rm_samepath", type = "logical", default = F, help = "if true do not remove common pathways")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFold <- args$inputFold
pathwayStructure_file <- args$pathwayStructure_file
sampleAnn_file <- args$sampleAnn_file
covDat_file <- args$covDat_file
thr_reliableGenes <- args$thr_reliableGenes
phenoDat_file <- args$phenoDat_file
names_file <- args$names_file
phenoAnn_file <- args$phenoAnn_file
cov_corr <- args$cov_corr
ncores <- args$ncores
geneAnn_file <- args$geneAnn_file
functR <- args$functR
geneSetName <- args$geneSetName
abs_tscore <- args$abs_tscore
not_rm_samepath <- args$not_rm_samepath
outFold <- args$outFold

# ##########################################################################################
# covDat_file <- 'INPUT_DATA/Covariates/scz_ersw_eur.covariateMatrix_old.txt'
# outFold <- 'OUTPUT_CMC/predict_PGC/200kb/scz_ersw_eur/devgeno0.01_testdevgeno0/'
# pathwayStructure_file <- '/home/luciat/priler_project/refData/SCZ_LoF_GeneSets_ordered.RData'
# thr_reliableGenes <- c(0.01,0)
# sampleAnn_file <- 'INPUT_DATA/Covariates/scz_ersw_eur.covariateMatrix_old.txt'
# inputFold <- 'OUTPUT_CMC/predict_PGC/200kb/scz_ersw_eur/devgeno0.01_testdevgeno0/'
# phenoDat_file <- 'INPUT_DATA/Covariates/scz_ersw_eur.phenoMatrix.txt'
# names_file <- 'SCZ_pheno'
# phenoAnn_file <- 'INPUT_DATA/Covariates/phenotypeDescription_PGCcohorts.csv'
# cov_corr <- T
# ncores <- 10
# geneSetName <- 'SCZ_LoF_GeneSets'
# geneAnn_file <- 'OUTPUT_CMC/train_CMC/200kb/resPrior_regEval_allchr.txt'
# functR <- '/home/luciat/priler_project/Software/model_prediction/pheno_association_functions.R'
# # ##########################################################################################


source(functR)

# load sample annotation
sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors=F)

# load gene annotation
geneAnn <- read.table(geneAnn_file, h=T, stringsAsFactors = F, sep = '\t')
geneAnn <- geneAnn[!(is.na(geneAnn$dev_geno) | is.na(geneAnn$test_dev_geno)), ]
geneAnn <- geneAnn[geneAnn$dev_geno >= thr_reliableGenes[1],]
geneAnn <- geneAnn[geneAnn$test_dev_geno > thr_reliableGenes[2],]

######################
#### load Tscore #####
######################

tscoreMat <- fread(sprintf('%spredictedTscores.txt', inputFold), header = T, stringsAsFactors = F, sep = '\t', check.names = F, data.table = F)
genesID <- tscoreMat[,1]
tscoreMat <- tscoreMat[,-1]
samplesID <- sapply(colnames(tscoreMat), function(x) strsplit(x, split = '.vs')[[1]][1])
samplesID <- unname(samplesID)
colnames(tscoreMat) <- samplesID
print(identical(samplesID, sampleAnn$Individual_ID))
tscoreMat <- t(tscoreMat) # samples on the rows

# remove sample that have NAs
id_s <- rowSums(is.na(tscoreMat)) == 0
sampleAnn <- sampleAnn[id_s,]
samplesID_new <- sampleAnn$Individual_ID
tscoreMat <- tscoreMat[id_s, ]

# filter geneAnn
if(!identical(genesID, geneAnn$external_gene_name)){
  print('adjust genes annotation')
  id <- sapply(genesID, function(x) which(x == geneAnn$external_gene_name))
  geneAnn <- geneAnn[id,]
}

print(paste('same gene annotation:', identical(geneAnn$external_gene_name, genesID)))

#if(!identical(geneAnn$external_gene_name, genesID)){
#  rep_gene <- names(which(table(geneAnn$external_gene_name)>1))
#  for(i in 1:length(rep_gene)){
#    tmp <- geneAnn[geneAnn$external_gene_name %in% rep_gene[i], ]
#    new_gene = data.frame(ensembl_gene_id = paste0(tmp$ensembl_gene_id, collapse = '_'), external_gene_name = rep_gene[i], dev_geno = mean(tmp$dev_geno), test_dev_geno = mean(tmp$test_dev_geno))
#    geneAnn <- geneAnn[! geneAnn$external_gene_name %in%  rep_gene[i], ]
#    geneAnn <- rbind(geneAnn, new_gene)
#  }
#  geneAnn <- geneAnn[match(genesID, geneAnn$external_gene_name),]
#}

print('Tscore mat loaded')

################################
#### load custom pathScore #####
################################
tmp <-  fread(sprintf('%sPathway_%s_scores.txt', inputFold, geneSetName), header = T, stringsAsFactors = F, sep = '\t', check.names = F, data.table = F)
pathScoreID <- tmp[,1]
tmp <- tmp[,-1]
identical(colnames(tmp), samplesID_new) # same order samples
pathScore <- as.matrix(t(tmp))
rm(tmp)

# consider only pathaways that do not have gene repetition, on pathwayScoreID add gene info
custom_pathway <- get(load(pathwayStructure_file))

path_name <- sapply(custom_pathway, function(x) x$name)
custom_pathway <- custom_pathway[which(path_name %in% pathScoreID)]
path_name <- sapply(custom_pathway, function(x) x$name)
identical(path_name, pathScoreID)

genes_path <- lapply(custom_pathway, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
ngenes_path <-  sapply(custom_pathway, function(x) length(unique(x$geneIds[x$geneIds != ""])))
ngenes_tscore_path <- sapply(genes_path, length)

# remove pathway with the same genes, recompute qvalue
if(!not_rm_samepath){
  rm_path <- c()
  len <- c()
  for(i in 1:length(genes_path)){
    # print(i)
    id <- which(sapply(genes_path, function(x) all(genes_path[[i]] %in% x) & all(x %in% genes_path[[i]])))
    len[i] <- length(id)
    ngenes_tmp <- ngenes_path[id]
    # take the one woth the lower amount of genes
    rm_path <- c(rm_path, path_name[id][-which.min(ngenes_tmp)])
  }
  
  rm_path <- unique(rm_path)
  id_rm <- which(pathScoreID %in% rm_path)
  if(length(id_rm)>0){
    pathScore <- pathScore[,-id_rm]
    pathScoreID <- pathScoreID[-id_rm]
    custom_pathway <- custom_pathway[-id_rm]
    path_name <- sapply(custom_pathway, function(x) x$name)
    identical(path_name, pathScoreID)
    genes_path <- lapply(custom_pathway, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
    ngenes_path <- ngenes_path[-id_rm]
    ngenes_tscore_path <- ngenes_tscore_path[-id_rm]
  }
}

# add number of genes info and mean test_dev_geno and dev_geno
df_path_info <- data.frame(path = pathScoreID, ngenes_tscore = unname(ngenes_tscore_path), ngenes_path = unname(ngenes_path), stringsAsFactors = F)
df_path_info$mean_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
df_path_info$sd_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
df_path_info$mean_test_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))
df_path_info$sd_test_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))

# compute mean/sd correlation for the genes belonging to the pathway (based on tscore)
df_path_info$mean_gene_corr <- NA
df_path_info$sd_gene_corr <- NA

for(i in 1:length(genes_path)){
  
  id <- which(genesID %in% genes_path[[i]])
  
  if(length(id)>1){
    # print(i)
    tmp <- cor(tscoreMat[,id])
    tmp <- tmp[lower.tri(tmp, diag = F)]
    df_path_info$mean_gene_corr[i] <- mean(tmp)
    df_path_info$sd_gene_corr[i] <- sd(tmp)  
  }
  
}

print('pathScore mat loaded')

if(abs_tscore){
  
  ##########################################
  #### load custom pathScore abstscore #####
  ##########################################
  
  tmp <-  fread(sprintf('%sPathway_%s_scores_abstscore.txt', inputFold, geneSetName), header = T, stringsAsFactors = F, sep = '\t', check.names = F, data.table = F)
  pathScoreID_abstscore <- tmp[,1]
  tmp <- tmp[,-1]
  identical(colnames(tmp), samplesID_new) # same order samples
  pathScore_abstscore <- as.matrix(t(tmp))
  rm(tmp)
  
  
  # consider only pathaways that do not have gene repetition, on pathwayScoreID add gene info
  custom_pathway <- get(load(pathwayStructure_file))
  
  path_name <- sapply(custom_pathway, function(x) x$name) 
  custom_pathway <- custom_pathway[which(path_name %in% pathScoreID_abstscore)]
  path_name <- sapply(custom_pathway, function(x) x$name)
  identical(path_name, pathScoreID_abstscore)
  
  genes_path <- lapply(custom_pathway, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
  ngenes_path <-  sapply(custom_pathway, function(x) length(unique(x$geneIds[x$geneIds != ""])))
  ngenes_tscore_path <- sapply(genes_path, length)
  
  # remove pathway with the same genes, recompute qvalue
  rm_path <- c()
  len <- c()
  for(i in 1:length(genes_path)){
    # print(i)
    id <- which(sapply(genes_path, function(x) all(genes_path[[i]] %in% x) & all(x %in% genes_path[[i]])))
    len[i] <- length(id)
    ngenes_tmp <- ngenes_path[id]
    # take the one woth the lower amount of genes
    rm_path <- c(rm_path, path_name[id][-which.min(ngenes_tmp)])
  }
  
  rm_path <- unique(rm_path)
  id_rm <- which(pathScoreID_abstscore %in% rm_path)
  if(length(id_rm)>0){
    pathScore_abstscore <- pathScore_abstscore[,-id_rm]
    pathScoreID_abstscore <- pathScoreID_abstscore[-id_rm]
    custom_pathway <- custom_pathway[-id_rm]
    path_name <- sapply(custom_pathway, function(x) x$name)
    identical(path_name, pathScoreID_abstscore)
    genes_path <- lapply(custom_pathway, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
    ngenes_path <- ngenes_path[-id_rm]
    ngenes_tscore_path <- ngenes_tscore_path[-id_rm]
  }
  
  # add number of genes info and mean test_dev_geno and dev_geno
  df_path_info_abstscore <- data.frame(path = pathScoreID_abstscore, ngenes_tscore = unname(ngenes_tscore_path), ngenes_path = unname(ngenes_path), stringsAsFactors = F)
  df_path_info_abstscore$mean_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_path_info_abstscore$sd_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_path_info_abstscore$mean_test_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))
  df_path_info_abstscore$sd_test_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))
  
  # compute mean/sd correlation for the genes belonging to the pathway (based on tscore)
  df_path_info_abstscore$mean_gene_corr <- NA
  df_path_info_abstscore$sd_gene_corr <- NA
  
  for(i in 1:length(genes_path)){
    
    id <- which(genesID %in% genes_path[[i]])
    
    if(length(id)>1){
      # print(i)
      tmp <- cor(tscoreMat[,id])
      tmp <- tmp[lower.tri(tmp, diag = F)]
      df_path_info_abstscore$mean_gene_corr[i] <- mean(tmp)
      df_path_info_abstscore$sd_gene_corr[i] <- sd(tmp)  
    }
    
  }
  
  print('pathScore mat abstscore loaded')
  
  
}


#### load phenotype annotation ####
phenoAnn <- fread(file = phenoAnn_file, header = T, stringsAsFactors = F, data.table = F)

###############################################################################################################################
for(n in 1:length(phenoDat_file)){
  
  phenoDat <- fread(file = phenoDat_file[n], header = T, stringsAsFactors = F, data.table = F, check.names=F)
  covDat <- fread(file = covDat_file[n], header = T, stringsAsFactors = F, data.table = F, check.names=F)
  
  print(paste0('############ phenoFile/covFile ', names_file[n], ' ############'))
  pheno_names <- colnames(phenoDat)[!colnames(phenoDat) %in% 'Individual_ID']
  print(pheno_names)
  if(cov_corr){
    tot_var <- merge(x = covDat[, ! colnames(covDat) %in% c('genoSample_ID', 'RNASample_ID', 'Dx')], y = phenoDat, by.x = 'Individual_ID', by.y = 'Individual_ID', sort = F)
    cov_names <- colnames(covDat[, ! colnames(covDat) %in% c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx')])
  }else{
    tot_var <- phenoDat
  }
  
  # filter samples
  tot_var <- tot_var[tot_var$Individual_ID %in% samplesID_new,]
  id_samples <- match(tot_var$Individual_ID, samplesID_new)
  samples_tmp <- samplesID_new[id_samples]
  print(identical(samples_tmp, tot_var$Individual_ID))
  tot_var <- tot_var[, ! colnames(tot_var) %in% 'Individual_ID']
  colnames(tot_var)[colnames(tot_var) %in% pheno_names] <- paste0('p', pheno_names)
  if(cov_corr){colnames(tot_var)[colnames(tot_var) %in% cov_names] <- paste0('c', cov_names)}
  print(colnames(tot_var))  
  
  # prepare phenotype
  phenoAnn_tmp <- phenoAnn[phenoAnn$pheno_id %in% pheno_names, c('pheno_id', 'Field', 'transformed_type')]
  print(str(phenoAnn_tmp))
  names_df <- t(sapply(phenoAnn_tmp$pheno_id, function(x) paste0(x, c('_beta', '_se_beta','_z_t','_pval'))))
  names_df_qval <- sapply(phenoAnn_tmp$pheno_id, function(x) paste0(x, c('_qval')))
  
  ############################
  #### tscore assocaition ####
  ############################
  tscoreMat_association <- cbind(tscoreMat[id_samples, ], tot_var)
  colnames(tscoreMat_association)[1:length(genesID)] <- paste0('X',1:length(genesID))
  df_tscore_info <- data.frame(ensembl_gene_id =  geneAnn$ensembl_gene_id, external_gene_name = geneAnn$external_gene_name,
                               dev_geno = geneAnn$dev_geno, test_dev_geno = geneAnn$test_dev_geno, stringsAsFactors = F)
  
  if(ncores >0){registerDoParallel(cores = ncores)}
  df_corr_tscore <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
  df_pi1 <- data.frame(pheno_id = phenoAnn_tmp$pheno_id,  tscore = rep(NA, nrow(phenoAnn_tmp)), pathScore = rep(NA, nrow(phenoAnn_tmp)), stringsAsFactors = F)
  if(abs_tscore){df_pi1$pathScore_abstscore <- rep(NA, nrow(phenoAnn_tmp))}  
  
  
  for(j in 1:nrow(phenoAnn_tmp)){
    
    print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))
    
    if(ncores>0){
      # parallelize over genes
      output <- foreach(x=1:nrow(df_tscore_info), .combine = rbind)%dopar%{
        # print(x)
        compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = tscoreMat_association)
      }   
    }else{
      output <- matrix(nrow = nrow(df_tscore_info), ncol = 4)
      for(x in 1:nrow(df_tscore_info)){
        # print(x)
        output[x,] <- compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = tscoreMat_association)
      }
    }
    
    # if(any(output == 'wrong pheno type annotation')){ stop("wrong pheno type annotation") }
    
    colnames(output) <- names_df[j,]
    rownames(output) <- NULL
    df_corr_tscore[[j]] <- cbind(df_tscore_info, output)
    
    # qvalue:
    qval <- tryCatch(qvalue(as.vector(df_corr_tscore[[j]][, names_df[j,4]])))
    df_corr_tscore[[j]] <- cbind(df_corr_tscore[[j]], qval$qvalue)
    colnames(df_corr_tscore[[j]])[ncol(df_corr_tscore[[j]])] <-  names_df_qval[j]
    df_pi1$tscore[j] <- 1 - qval$pi0
    
  }
  
  rm(tscoreMat_association)
  print('tscore completed')
  
  ##########################
  #### pathScore custom ####
  ##########################
  
  pathScore_association <- cbind(pathScore[id_samples, ], tot_var)
  colnames(pathScore_association)[1:length(pathScoreID)] <- paste0('X',1:length(pathScoreID))
  
  df_corr_path <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
  info_path <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
  
  for(j in 1:nrow(phenoAnn_tmp)){
    
    print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))
    
    if(ncores>0){
      # parallelize over genes
      output <- foreach(x=1:nrow(df_path_info), .combine = rbind)%dopar%{
        # print(x)
        compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_association)
      }   
    }else{
      output <- matrix(nrow = nrow(df_path_info), ncol = 4)
      for(x in 1:nrow(df_path_info)){
        print(x)
        output[x,] <- compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_association)
      }
    }
    
    colnames(output) <- names_df[j,]
    rownames(output) <- NULL
    df_corr_path[[j]] <- cbind(df_path_info, output)
    
    # qvalue:
    qval <- tryCatch(qvalue(as.vector(df_corr_path[[j]][, names_df[j,4]])), error=function(...) NULL)
    if(is.null(qval)){
      qval <- list(qvalue = rep(NA, nrow(df_corr_path[[j]])), pi0 = NA)  
    }
    df_corr_path[[j]] <- cbind(df_corr_path[[j]], qval$qvalue)
    colnames(df_corr_path[[j]])[ncol(df_corr_path[[j]])] <-  names_df_qval[j]
    df_pi1$pathScore[j] <- 1 - qval$pi0
    
    # create list object for each pathway with tscore association results
    info_path[[j]] <- vector(mode = 'list', length = nrow(df_path_info))
    for(x in 1:nrow(df_path_info)){
      
      tmp <- custom_pathway[which(sapply(custom_pathway, function(x) x$name) == df_path_info$path[x])]
      info_path[[j]][[x]] <- list(path = df_corr_path[[j]][x,], genes_path = tmp[[1]]$geneIds, tscore = df_corr_tscore[[j]][df_corr_tscore[[j]]$external_gene_name %in% tmp[[1]]$geneIds,])
      
    }
    
  }
  
  rm(pathScore_association)
  print('custom pathScore completed')
  
  if(abs_tscore){
    ####################################
    #### pathScore custom abstscore ####
    ####################################
    
    pathScore_association <- cbind(pathScore_abstscore[id_samples, ], tot_var)
    colnames(pathScore_association)[1:length(pathScoreID_abstscore)] <- paste0('X',1:length(pathScoreID_abstscore))
    
    df_corr_path_abstscore <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
    info_path_abstscore <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
    
    for(j in 1:nrow(phenoAnn_tmp)){
      
      print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))
      
      if(ncores>0){
        # parallelize over genes
        output <- foreach(x=1:nrow(df_path_info_abstscore), .combine = rbind)%dopar%{
          # print(x)
          compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_association)
        }   
      }else{
        output <- matrix(nrow = nrow(df_path_info_abstscore), ncol = 4)
        for(x in 1:nrow(df_path_info_abstscore)){
          print(x)
          output[x,] <- compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_association)
        }
      }
      
      colnames(output) <- names_df[j,]
      rownames(output) <- NULL
      df_corr_path_abstscore[[j]] <- cbind(df_path_info_abstscore, output)
      
      # qvalue:
      qval <- qvalue(as.vector(df_corr_path_abstscore[[j]][, names_df[j,4]]))
      df_corr_path_abstscore[[j]] <- cbind(df_corr_path_abstscore[[j]], qval$qvalue)
      colnames(df_corr_path_abstscore[[j]])[ncol(df_corr_path_abstscore[[j]])] <-  names_df_qval[j]
      df_pi1$pathScore_abstscore[j] <- 1 - qval$pi0
      
      # create list object for each pathway with tscore association results
      info_path_abstscore[[j]] <- vector(mode = 'list', length = nrow(df_path_info_abstscore))
      for(x in 1:nrow(df_path_info_abstscore)){
        
        tmp <- custom_pathway[which(sapply(custom_pathway, function(x) x$name) == df_path_info_abstscore$path[x])]
        info_path_abstscore[[j]][[x]] <- list(path = df_corr_path_abstscore[[j]][x,], genes_path = tmp[[1]]$geneIds, tscore = df_corr_tscore[[j]][df_corr_tscore[[j]]$external_gene_name %in% tmp[[1]]$geneIds,])
        
      }
      
    }
    
    rm(pathScore_association)
    print('custom pathScore abstscore completed')
  }
  ################################################################################################################################################################################
  
  final <- list(pheno = phenoAnn[phenoAnn$pheno_id %in% pheno_names, ], tscore = df_corr_tscore, pathScore = df_corr_path, pi1 = df_pi1, info_pathScore = info_path)
  if(abs_tscore){
    final$pathScore_abstscore <- df_corr_path_abstscore
    final$info_pathScore_abstscore <- info_path_abstscore
  }  
  
  # save results
  filename <- ifelse(cov_corr, sprintf('%spval_%s_covCorr_customPath_%s.RData', outFold, names_file[n], geneSetName), sprintf('%spval_%s_customPath_%s.RData', outFold, names_file[n], geneSetName))
  print(filename)
  save(final, file = filename)
  
  print('final result saved')
  
}

