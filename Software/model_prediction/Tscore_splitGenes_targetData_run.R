#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# compute T-scores wrt a reference set, UKBB sample size

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))


parser <- ArgumentParser(description="Differential pathway analysis for big sample size, use t-statistic instead of moderate one, split per genes")

parser$add_argument("--input_file", type = "character", nargs = '*', help = "output path prediction folder (or cov corrected data), splitted files")
parser$add_argument("--covDat_notref_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--covDat_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases), must include covDat_notref_file samples")
parser$add_argument("--nFolds", type="integer", default = 20, help = "number of fold for control partition")
parser$add_argument("--perc_comp", type="double", default = 0.5, help = "percentage of samples used for comparison")
parser$add_argument("--split_gene_id", type="integer", help = "split partiton for genes")
parser$add_argument("--split_tot", type="integer", default = 100, help = "total number of partitions")
parser$add_argument("--ncores", type="integer", default = 10, help = "n cores for the parallelization")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
input_file <- args$input_file
covDat_notref_file <- args$covDat_notref_file
covDat_file <- args$covDat_file
perc_comp <- args$perc_comp
nFolds <- args$nFolds
split_gene_id <- args$split_gene_id
split_tot <- args$split_tot
ncores <- args$ncores
outFold <- args$outFold

# ####################################################################
# input_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/split', 1:100, '_predictedExpression_filt.txt')
# nFolds <- 10
# perc_comp <- 0.7
# ncores <- 10
# covDat_file <- '/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/covariateMatrix_LipoProteinDisorder.txt'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# split_gene_id <- 1
# split_tot <- 100
# ###################################################################

#########################
### function needed #####
#########################
# function: extract correct ref matrix 
extract_values <- function(inputData, colnames_exp_list, sample_list){
  
  if(length(sample_list) == 1){
    
    id_list = which(sapply(colnames_exp_list, function(x) sample_list %in% x))
    new_mat <- inputData[[id_list]][, colnames_exp_list[[id_list]] %in% sample_list]
    
  }else{
    
    tmp <- list()
    for(i in 1:length(inputData)){
      
      if(any(colnames_exp_list[[i]] %in% sample_list)){
        id <- which(colnames_exp_list[[i]] %in% sample_list)
        tmp[[i]] <- matrix(inputData[[i]][, id], ncol = length(id), nrow = nrow(inputData[[i]]))
        colnames(tmp[[i]]) <- colnames_exp_list[[i]][id]
      }
    }
    
    tmp <- do.call(cbind, tmp)
    tmp <- tmp[, match(sample_list, colnames(tmp))]
    print(identical(colnames(tmp), sample_list))
    
    new_mat <- big.matrix(ncol=length(sample_list) , nrow = nrow(inputData[[i]]), type = "double", init = 0, dimnames = NULL, shared = F)
    new_mat[,] <-  tmp
    
  } 
  return(new_mat)
  
}


# function to compute t statistic
t_stat <- function(x){
  n <- ncol(x)
  mu <- rowMeans(x)
  sigma <- rowSds(x)
  return(mu/(sigma/sqrt(n)))
}

### load data ###
### NOTE: if splitted the genes should be already filtered for gene performance
sampleAnn_tot <- read.table(covDat_file, header = T, stringsAsFactors = F)
sampleAnn_tot$ref <- T
sampleAnn_notref <- read.table(covDat_notref_file, h=T,stringsAsFactors = F)
sampleAnn_tot$ref[sampleAnn_tot$Individual_ID %in% sampleAnn_notref$Individual_ID] <- F

# load block by block, store as big matrix
expDat <-  fread(input_file[1],  header = T, stringsAsFactors = F, data.table = F) 
common_samples <- intersect(colnames(expDat), sampleAnn_tot$Individual_ID)
id_samples <- match(common_samples, colnames(expDat))

colnames_exp <- c()
colnames_exp_list <- vector(mode = 'list', length = length(input_file))
eMat <- vector(mode = 'list', length = length(input_file))

# consider only some genes:
step_genes <- round(nrow(expDat)/split_tot)
split_id = as.vector(sapply(1:split_tot, function(x) rep(x, step_genes)))
if(length(split_id)<nrow(expDat)){
  split_id <- c(split_id, rep(split_tot, nrow(expDat) - length(split_id)))  
}else{
  if(length(split_id)>nrow(expDat)){
    split_id <- split_id[-((nrow(expDat)+1):length(split_id))]
  }
}

split_df <- data.frame(id = 1:nrow(expDat), split = split_id)
print(tail(split_df))
gene_id <- split_df$id[split_df$split == split_gene_id]
geneInfo <- expDat[gene_id, -id_samples]

colnames_exp <- c(colnames_exp, colnames(expDat[, id_samples]))
colnames_exp_list[[1]] <- colnames(expDat[, id_samples])

print(str(id_samples))
print(str(geneInfo))

eMat[[1]] <- big.matrix(ncol=length(id_samples) , nrow = nrow(geneInfo), type = "double", init = 0, dimnames = NULL, shared = F)
eMat[[1]][,] <- as.matrix(expDat[gene_id,id_samples])

for(i in 2:length(input_file)){
  
  print(i)
  
  expDat <-  fread(input_file[i],  header = T, stringsAsFactors = F, data.table = F) 
  common_samples <- intersect(colnames(expDat), sampleAnn_tot$Individual_ID)
  id_samples <- match(common_samples, colnames(expDat))
  
  if(!identical(geneInfo$ensembl_gene_id, expDat[gene_id, 'ensembl_gene_id'])){print("ERROR: different gene annotation")}
  colnames_exp_list[[i]] <- colnames(expDat[, id_samples])
  colnames_exp <- c(colnames_exp, colnames(expDat[, id_samples]))
  
  eMat[[i]] <- big.matrix(ncol=length(id_samples) , nrow = nrow(geneInfo), type = "double", init = 0, dimnames = NULL, shared = F)
  eMat[[i]][,] <- as.matrix(expDat[gene_id,id_samples])
  
  if(i%%20 == 0){print(mem_used())}
  
}

# reorder total annotation file based on the big.matrix
if(!identical(sampleAnn_tot$Individual_ID, colnames_exp)){
  id_samples <- match(colnames_exp, sampleAnn_tot$Individual_ID)
  sampleAnn_tot <- sampleAnn_tot[id_samples,]
  print(identical(sampleAnn_tot$Individual_ID, colnames_exp))
}


# change Individual_ID name if start with a number or '-' present
sampleAnn <- sampleAnn_tot
sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) == 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Individual_ID[id_n], function(x) paste0('X', x))   
}
colnames_exp <- sampleAnn$Temp_ID

# add Dx info if not present
if(!'Dx' %in% colnames(sampleAnn)){sampleAnn$Dx <- 0}
sampleAnn$split_id = unlist(lapply(1:length(colnames_exp_list), function(x) rep(x, length(colnames_exp_list[[x]]))))


# STEP 1: LOAD AND CONFIGURE INPUT DATA FOR SCORECARD CALCULATION
sampleIds = colnames_exp
curType = "Gene_expression"
geneNames = geneInfo$external_gene_name
ind <- geneNames==""
geneNames <- geneNames[!ind]
# remove genes for each eMat file
if(length(geneNames)==1){
  inputData <- vector(length = length(eMat), mode = 'list')
  for(i in 1:length(inputData)){
    inputData[[i]] <- matrix(eMat[[i]][!ind,], nrow = length(geneNames), ncol=ncol(eMat[[i]]))
  }
}else{
  inputData <- lapply(eMat, function(x) as.big.matrix(x[!ind,], type = 'double', shared = F))
}

rm(eMat)

sampleGroups = list()
comparisons = list()
sampleGroups[["all"]] = sampleIds
sampleGroups[["Case"]] = sampleAnn$Temp_ID[!sampleAnn$ref] 
sampleGroups[["Ref"]] = sampleAnn$Temp_ID[sampleAnn$Dx==0 & sampleAnn$ref] 
comparisons[[1]] = c("pairwise","all","Ref")

# produce a genes*samples table of t-scores for each comparison
for (i in 1:length(sampleGroups)){sampleGroups[[i]] = unique(sampleGroups[[i]])}

curComparison = comparisons[[1]]
print(paste("Processing: ",paste(curComparison,collapse="::"),sep=""))
#allSamp = sampleGroups[[curComparison[[2]]]] # all
#reference = sampleGroups[[curComparison[[3]]]] # ref (controls)
#overlap = intersect(allSamp,reference) # ref
overlap = sampleGroups[[curComparison[[3]]]]
nSel=floor(length(overlap)*perc_comp) #floor(length(overlap)/nFolds)

# not balanced 
ovMat=sapply(seq(1,nFolds),function(X){
  
  vec=rep(F,length(overlap))
  set.seed(42+X)
  vec[sample.int(length(overlap),nSel)]=T
  return(vec)
  
})

rownames(ovMat)=overlap

# res must be a bigmatrix
# parallelize wrt nFolds

# extract outside ref and comparison
excluded <- lapply(1:ncol(ovMat), function(x) rownames(ovMat)[ovMat[,x]])
# curCases <- lapply(excluded, function(x) union(x,setdiff(allSamp,overlap))) # union part of controls and cases
curCases <- sampleGroups[["Case"]] 
# curReference <- lapply(excluded, function(x) setdiff(reference,x))
curReference <- lapply(excluded, function(x) setdiff(sampleGroups[["Ref"]],x))

print(mem_used())

res_mat <- vector(mode = 'list', length = nFolds)
# res_mat <- foreach(i=1:nFolds)%dopar%{
for (i in 1:nFolds){
  
  # excluded = rownames(ovMat)[ovMat[,i]]
  comparisonLabel = paste("Reference_excluding_fold",i,sep="_")
  print(comparisonLabel)
  
  ref <- extract_values(inputData, colnames_exp_list, curReference[[i]])
  comp <- extract_values(inputData, colnames_exp_list, curCases)
  
  if(length(geneNames)==1){
    ref <- matrix(ref[,], nrow = 1)
    comp <- matrix(comp[,], nrow = 1)
  }else{
    ref <- ref[,]
    comp <- comp[,]
  }
  
  print(dim(comp))
  print(dim(ref))
  
  cl <- makeCluster(ncores)
  clusterExport(cl, "comp")
  clusterExport(cl, "ref")
  clusterExport(cl, "t_stat")
  clusterEvalQ(cl, library("matrixStats"))
  
  tmp <- parLapply(cl, 1:length(curCases), function(j) t_stat(comp[,j] - ref))
  stopCluster(cl)
  
  res_mat[[i]] <- do.call(cbind, tmp)
  
  colnames(res_mat[[i]]) <- curCases[[i]]
  rownames(res_mat[[i]]) <- geneNames
  
}

print('Tscore computation finished: reorganize')

tscoreTables <- data.frame(geneId = geneNames)
cl <- makeCluster(ncores)
clusterExport(cl, "res_mat")
clusterExport(cl, "curCases")
clusterExport(cl, "geneNames")

tmp <- parLapply(cl, 1:length(curCases), function(s){
  new <- sapply(res_mat, function(x) x[,s]) 
  rowMeans(matrix(new, nrow = length(geneNames), ncol = length(res_mat)))
})
stopCluster(cl)

tscoreTables <- do.call(cbind,tmp)
colnames(tscoreTables) <- curCases
tscoreTables <- data.frame(geneId = geneNames, tscoreTables)

# save results (.RData)
save(tscoreTables, file = sprintf('%spredictedTscores_splitGenes%i.RData', outFold, split_gene_id))

