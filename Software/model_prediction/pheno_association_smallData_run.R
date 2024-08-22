#!/usr/bin/env Rscript
# v3: save also beta, se, zcores, phenotype already processed

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))


parser <- ArgumentParser(description="Gene and Pathwasy association analysis")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data")
parser$add_argument("--thr_reliableGenes", type = "double", nargs = '*', default = c(0.01, 0), help = "threshold for reliable genes: dev_geno_tot and test_dev_geno  [default %(default)s]")
parser$add_argument("--inputFold", type = "character", help = "Folde with results from pathway analysis")
parser$add_argument("--covDat_file", type = "character", nargs = '*', help = "file with sample info for new genotype data, associated to phenoDat_file")
parser$add_argument("--phenoDat_file", type = "character", nargs = '*', help = "file(s) with individual_ID to match and phenotype to test association, associated to covDat_file")
parser$add_argument("--names_file", type = "character", nargs = '*', help = "for each couple of covDat/phenoDat file, associated name")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical", default = T,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)  [default %(default)s]")
parser$add_argument("--ncores", type = "integer", default = 0, help = "number of cores for the parallelization (over genes, if zero not parallelized)  [default %(default)s]")
parser$add_argument("--geneAnn_file", type = "character", help = "file with gene info from train, to be filtered")
parser$add_argument("--functR", type = "character", help = "Rscript with functions to be used")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFold <- args$inputFold
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
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
outFold <- args$outFold

##########################################################################################
# inputFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/CG/devgeno0.01_testdevgeno0/'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/CG/devgeno0.01_testdevgeno0/'
# GOterms_file <- '/psycl/g/mpsziller/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/psycl/g/mpsziller/lucia/refData/ReactomePathways.gmt'
# thr_reliableGenes <- c(0.01,0)
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/CG/covariateMatrix.txt'
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/CG/covariateMatrix.txt'
# phenoDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/CG/phenoMatrix.txt'
# names_file <- 'CAD_pheno'
# phenoAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/phenotypeDescription_SchunkertCohorts.csv'
# cov_corr <- T
# ncores <- 10
# geneAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/resPrior_regEval_allchr.txt'
# functR <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/RSCRIPTS/SCRIPTS_v2/AssociationAnalysis_functions_run.R'
##########################################################################################

source(functR)

# load sample annotation
sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors=F)
# change Individual_ID name if start with a number or '-' or '*' present
#sampleAnn$Temp_ID <- sampleAnn$Individual_ID
#id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) >= 2)
#id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
#id_a <-  sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '[*]')[[1]]) >= 2)
#if(any(id_h)){
#  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
#}
#if(any(id_n)){
#  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Temp_ID[id_n], function(x) paste0('X', x))   
#}
#if(any(id_a)){
#  sampleAnn$Temp_ID[id_a] <- sapply(sampleAnn$Temp_ID[id_a], function(x) paste0(strsplit(x, split = '[*]')[[1]], collapse = '_'))   
#}

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
#samplesID <- sapply(colnames(tscoreMat), function(x) strsplit(x, split = '.vs')[[1]][1])
#samplesID <- unname(samplesID)
#colnames(tscoreMat) <- samplesID
print(identical(colnames(tscoreMat), sampleAnn$Individual_ID))
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

##################################
#### load pathScore Reactome #####
##################################
if (!is.null(reactome_file)) {
  tmp <-  fread(sprintf('%sPathway_Reactome_scores.txt', inputFold), header = T, stringsAsFactors = F, sep = '\t', check.names = F, data.table = F)
  pathScoreID_reactome <- tmp[,1]
  tmp <- tmp[,-1]
  identical(colnames(tmp), samplesID_new) # same order samples
  pathScore_reactome <- as.matrix(t(tmp))
  rm(tmp)


  # consider only pathaways that do not have gene repetition, on pathwayScoreID add gene info
  gs <- readGmt(reactome_file)

  gs=lapply(gs,function(X){
    X@ids=gsub(",1.0","",X@ids)
    return(X)
  })
  gs_name <- sapply(gs, function(x) x@reference)
  gs <- gs[which(gs_name %in% pathScoreID_reactome)]
  gs_name <- unname(sapply(gs, function(x) x@reference))
  identical(gs_name, pathScoreID_reactome)

  genes_path <- lapply(gs, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x@ids])
  ngenes_path <-  sapply(gs, function(x) length(unique(x@ids[x@ids != ""])))
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
    rm_path <- c(rm_path, gs_name[id][-which.min(ngenes_tmp)])
  }

  rm_path <- unique(rm_path)
  id_rm <- which(pathScoreID_reactome %in% rm_path)
  if(length(id_rm)>0){
    pathScore_reactome <- pathScore_reactome[,-id_rm]
    pathScoreID_reactome <- pathScoreID_reactome[-id_rm]
    gs <- gs[-id_rm]
    gs_name <- unname(sapply(gs, function(x) x@reference))
    identical(gs_name, pathScoreID_reactome)
    genes_path <- lapply(gs, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x@ids])
    ngenes_path <- ngenes_path[-id_rm]
    ngenes_tscore_path <- ngenes_tscore_path[-id_rm]
  }

  # add number of genes info and mean test_dev_geno and dev_geno
  df_pathR_info <- data.frame(path = pathScoreID_reactome, ngenes_tscore = unname(ngenes_tscore_path), ngenes_path = unname(ngenes_path), stringsAsFactors = F)
  df_pathR_info$mean_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathR_info$sd_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathR_info$mean_test_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathR_info$sd_test_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))

  # compute mean/sd correlation for the genes belonging to the pathway (based on tscore)
  df_pathR_info$mean_gene_corr <- NA
  df_pathR_info$sd_gene_corr <- NA

  for(i in 1:length(genes_path)){

    id <- which(genesID %in% genes_path[[i]])

    if(length(id)>1){
      # print(i)
      tmp <- cor(tscoreMat[,id])
      tmp <- tmp[lower.tri(tmp, diag = F)]
      df_pathR_info$mean_gene_corr[i] <- mean(tmp)
      df_pathR_info$sd_gene_corr[i] <- sd(tmp)
    }

  }

  print('pathScore reactome mat loaded')
}

############################
#### load pathScore GO #####
############################
if (!is.null(GOterms_file)) {
  tmp <-  fread(sprintf('%sPathway_GO_scores.txt', inputFold), header = T, stringsAsFactors = F, sep = '\t', check.names = F, data.table = F)
  pathScoreID_GO <- tmp[,1]
  tmp <- tmp[,-1]
  identical(colnames(tmp), samplesID_new) # same order samples
  pathScore_GO <- as.matrix(t(tmp))
  rm(tmp)

  # consider only pathaways that do not heave gene repetition, on pathwayScoreID add gene info
  go <- get(load(GOterms_file))

  go_name <- sapply(go, function(x) x$GOID)
  go <- go[which(go_name %in% pathScoreID_GO)]
  go_name <- sapply(go, function(x) x$GOID)
  go_path_name <- sapply(go, function(x) x$Term)
  go_ont_name <- sapply(go, function(x) x$Ontology)
  identical(go_name, pathScoreID_GO)

  genes_path <- lapply(go, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
  ngenes_path <-  sapply(go, function(x) length(unique(x$geneIds[x$geneIds != ""])))
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
    rm_path <- c(rm_path, go_name[id][-which.min(ngenes_tmp)])
  }

  rm_path <- unique(rm_path)
  id_rm <- which(pathScoreID_GO %in% rm_path)
  if(length(id_rm)>0){
    pathScore_GO <- pathScore_GO[,-id_rm]
    pathScoreID_GO <- pathScoreID_GO[-id_rm]
    go <- go[-id_rm]
    go_name <- sapply(go, function(x) x$GOID)
    go_path_name <- sapply(go, function(x) x$Term)
    go_ont_name <- sapply(go, function(x) x$Ontology)
    identical(go_name, pathScoreID_GO)
    genes_path <- lapply(go, function(x) geneAnn$external_gene_name[geneAnn$external_gene_name %in% x$geneIds])
    ngenes_path <- ngenes_path[-id_rm]
    ngenes_tscore_path <- ngenes_tscore_path[-id_rm]
  }

  # add number of genes info and mean test_dev_geno and dev_geno
  df_pathGO_info <- data.frame(path_id = pathScoreID_GO, path = go_path_name, path_ont = go_ont_name, ngenes_tscore = unname(ngenes_tscore_path), ngenes_path = unname(ngenes_path), stringsAsFactors = F)
  df_pathGO_info$mean_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathGO_info$sd_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathGO_info$mean_test_dev_geno <- unname(sapply(genes_path, function(x) mean(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))
  df_pathGO_info$sd_test_dev_geno <- unname(sapply(genes_path, function(x) sd(geneAnn$test_dev_geno[geneAnn$external_gene_name %in% x])))

  # compute mean/sd correlation for the genes belonging to the pathway (based on tscore)
  df_pathGO_info$mean_gene_corr <- NA
  df_pathGO_info$sd_gene_corr <- NA

  for(i in 1:length(genes_path)){

    id <- which(genesID %in% genes_path[[i]])

    if(length(id)>1){
      # print(i)
      tmp <- cor(tscoreMat[,id])
      tmp <- tmp[lower.tri(tmp, diag = F)]
      df_pathGO_info$mean_gene_corr[i] <- mean(tmp)
      df_pathGO_info$sd_gene_corr[i] <- sd(tmp)
    }

  }

  print('pathScore GO mat loaded')
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
    cov_names <- colnames(covDat[, ! colnames(covDat) %in% c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx'), drop = FALSE])
  }else{
    tot_var <- phenoDat
  }

  ## change Individual_ID name if start with a number or '-' or '*' present
  #tot_var$Temp_ID <- tot_var$Individual_ID
  #id_h <- sapply(tot_var$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) >= 2)
  #id_n <- sapply(tot_var$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
  #id_a <-  sapply(tot_var$Individual_ID, function(x) length(strsplit(x, split = '[*]')[[1]]) >= 2)
  #if(any(id_h)){
  #  tot_var$Temp_ID[id_h] <- sapply(tot_var$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
  #}
  #if(any(id_n)){
  #  tot_var$Temp_ID[id_n] <- sapply(tot_var$Temp_ID[id_n], function(x) paste0('X', x))   
  #}
  #if(any(id_a)){
  #  tot_var$Temp_ID[id_a] <- sapply(tot_var$Temp_ID[id_a], function(x) paste0(strsplit(x, split = '[*]')[[1]], collapse = '_'))   
  #}
  
  # filter samples
  tot_var <- tot_var[tot_var$Individual_ID %in% samplesID_new,]
  id_samples <- match(tot_var$Individual_ID, samplesID_new)
  samples_tmp <- samplesID_new[id_samples]
  print(identical(samples_tmp, tot_var$Individual_ID))
  tot_var <- tot_var[, ! colnames(tot_var) %in% 'Individual_ID', drop = FALSE]
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
  df_pi1 <- data.frame(pheno_id = phenoAnn_tmp$pheno_id,  tscore = rep(NA, nrow(phenoAnn_tmp)), pathScore_reactome = rep(NA, nrow(phenoAnn_tmp)), pathScore_GO = rep(NA, nrow(phenoAnn_tmp)), stringsAsFactors = F)
  
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
    qval <- qvalue(as.vector(df_corr_tscore[[j]][, names_df[j,4]]))
    df_corr_tscore[[j]] <- cbind(df_corr_tscore[[j]], qval$qvalue)
    colnames(df_corr_tscore[[j]])[ncol(df_corr_tscore[[j]])] <-  names_df_qval[j]
    df_pi1$tscore[j] <- 1 - qval$pi0
    
  }
  
  rm(tscoreMat_association)
  print('tscore completed')
  
  ########################################
  #### pathScore Reactome assocaition ####
  ########################################
  if (!is.null(reactome_file)) {
    pathScore_r_association <- cbind(pathScore_reactome[id_samples, ], tot_var)
    colnames(pathScore_r_association)[1:length(pathScoreID_reactome)] <- paste0('X',1:length(pathScoreID_reactome))

    df_corr_pathR <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
    info_pathR <- vector(mode = 'list', length = nrow(phenoAnn_tmp))

    for(j in 1:nrow(phenoAnn_tmp)){

      print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))

      if(ncores>0){
        # parallelize over genes
        output <- foreach(x=1:nrow(df_pathR_info), .combine = rbind)%dopar%{
          # print(x)
          compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_r_association)
        }
      }else{
        output <- matrix(nrow = nrow(df_pathR_info), ncol = 4)
        for(x in 1:nrow(df_pathR_info)){
          # print(x)
          output[x,] <- compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_r_association)
        }
      }

      colnames(output) <- names_df[j,]
      rownames(output) <- NULL
      df_corr_pathR[[j]] <- cbind(df_pathR_info, output)

      # qvalue:
      qval <- qvalue(as.vector(df_corr_pathR[[j]][, names_df[j,4]]))
      df_corr_pathR[[j]] <- cbind(df_corr_pathR[[j]], qval$qvalue)
      colnames(df_corr_pathR[[j]])[ncol(df_corr_pathR[[j]])] <-  names_df_qval[j]
      df_pi1$pathScore_reactome[j] <- 1 - qval$pi0

      # create list object for each pathway with tscore association results
      info_pathR[[j]] <- vector(mode = 'list', length = nrow(df_pathR_info))
      for(x in 1:nrow(df_pathR_info)){

        tmp <- gs[which(sapply(gs, function(x) x@reference) == df_pathR_info$path[x])]
        info_pathR[[j]][[x]] <- list(path = df_corr_pathR[[j]][x,], genes_path = tmp[[1]]@ids, tscore = df_corr_tscore[[j]][df_corr_tscore[[j]]$external_gene_name %in% tmp[[1]]@ids,])

      }

    }

    rm(pathScore_r_association)
    print('pathScore reactome completed')
  } else {
    df_corr_pathR <- data.frame()
    info_pathR <- data.frame()
  }
  
  ##################################
  #### pathScore GO assocaition ####
  ##################################
  if (!is.null(reactome_file)) {
    pathScore_GO_association <- cbind(pathScore_GO[id_samples, ], tot_var)
    colnames(pathScore_GO_association)[1:length(pathScoreID_GO)] <- paste0('X',1:length(pathScoreID_GO))

    df_corr_pathGO <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
    info_pathGO <- vector(mode = 'list', length = nrow(phenoAnn_tmp))

    for(j in 1:nrow(phenoAnn_tmp)){

      print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))

      if(ncores>0){
        # parallelize over genes
        output <- foreach(x=1:nrow(df_pathGO_info), .combine = rbind)%dopar%{
          # print(x)
          compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_GO_association)
        }
      }else{
        output <- matrix(nrow = nrow(df_pathGO_info), ncol = 4)
        for(x in 1:nrow(df_pathGO_info)){
          # print(x)
          output[x,] <- compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScore_GO_association)
        }
      }

      colnames(output) <- names_df[j,]
      rownames(output) <- NULL
      df_corr_pathGO[[j]] <- cbind(df_pathGO_info, output)

      # qvalue:
      qval <- qvalue(as.vector(df_corr_pathGO[[j]][, names_df[j,4]]))
      df_corr_pathGO[[j]] <- cbind(df_corr_pathGO[[j]], qval$qvalue)
      colnames(df_corr_pathGO[[j]])[ncol(df_corr_pathGO[[j]])] <-  names_df_qval[j]
      df_pi1$pathScore_GO[j] <- 1 - qval$pi0

      # create list object for each pathway with tscore association results
      # create list object for each pathway with tscore association results
      info_pathGO[[j]] <- vector(mode = 'list', length = nrow(df_pathGO_info))
      for(x in 1:nrow(df_pathGO_info)){

        tmp <- go[which(sapply(go, function(x) x$GOID) == df_pathGO_info$path_id[x])]
        info_pathGO[[j]][[x]] <- list(path = df_corr_pathGO[[j]][x,], genes_path = tmp[[1]]$geneIds, tscore = df_corr_tscore[[j]][df_corr_tscore[[j]]$external_gene_name %in% tmp[[1]]$geneIds,])

      }
    }

    rm(pathScore_GO_association)
    print('pathScore go completed')
  } else {
    df_corr_pathGO <- data.frame()
    info_pathGO <- data.frame()
  }

  ################################################################################################################################################################################
  
  final <- list(pheno = phenoAnn[phenoAnn$pheno_id %in% pheno_names, ], tscore = df_corr_tscore, pathScore_reactome = df_corr_pathR, pathScore_GO = df_corr_pathGO, pi1 = df_pi1, 
                info_pathScore_reactome = info_pathR, info_pathScore_GO = info_pathGO)
  
  # save results
  filename <- ifelse(cov_corr, sprintf('%spval_%s_covCorr.RData', outFold, names_file[n]), sprintf('%spval_%s.RData', outFold, names_file[n]))
  print(filename)
  save(final, file = filename)
  
  print('final result saved')
  
}

