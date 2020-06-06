#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser(description="split data in base and target to compute PRS, note: target should not contain samples used for reference")

parser$add_argument("--covDat_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--input_file", type = "character", nargs = '*', help = "total path predicted gene expression, splitted files")
parser$add_argument("--nFolds", type="integer", default = 20, help = "number of fold for control partition")
parser$add_argument("--perc_comp", type="double", default = 0.5, help = "percentage of samples used for comparison")
parser$add_argument("--perc_ext", type="double", default = 0.2, help = "percentage of samples to keep as external")
parser$add_argument("--seed_fixed", type="integer", help = "seed for sample extraction")
parser$add_argument("--thr_n", type="integer", default = 10000, help = "maximum number of elements for a class")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
covDat_file <- args$covDat_file
input_file <- args$input_file
perc_comp <- args$perc_comp
nFolds <- args$nFolds
seed_fixed <- args$seed_fixed
thr_n <- args$thr_n
perc_ext <- args$perc_ext
outFold <- args$outFold

# ################################################################################################################################
# input_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/split', 1:100, '_predictedExpression_filt.txt')
# nFolds <- 10
# perc_comp <- 0.7
# perc_ext <- 0.2
# seed_fixed <- 1234
# thr_n <- 10000
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/covariateMatrix_latestW.txt'
# ################################################################################################################################

sampleAnn_tot <- read.table(covDat_file, header = T, stringsAsFactors = F)

colnames_exp <- c()
colnames_exp_list <- vector(mode = 'list', length = length(input_file))

expDat <-  fread(input_file[1],  header = T, stringsAsFactors = F, data.table = F) 
common_samples <- intersect(colnames(expDat), sampleAnn_tot$Individual_ID)
id_samples <- match(common_samples, colnames(expDat))
colnames_exp <- c(colnames_exp, colnames(expDat[, id_samples]))
colnames_exp_list[[1]] <- colnames(expDat[, id_samples])

# for(i in 2:10){
for(i in 2:length(input_file)){
  
  print(i)
  
  expDat <-  fread(input_file[i],  header = T, stringsAsFactors = F, data.table = F) 
  common_samples <- intersect(colnames(expDat), sampleAnn_tot$Individual_ID)
  id_samples <- match(common_samples, colnames(expDat))
  colnames_exp_list[[i]] <- colnames(expDat[, id_samples])
  colnames_exp <- c(colnames_exp, colnames(expDat[, id_samples]))
}

# reorder total annotation file based on the big.matrix
if(!identical(sampleAnn_tot$Individual_ID, colnames_exp)){
  id_samples <- match(colnames_exp, sampleAnn_tot$Individual_ID)
  sampleAnn_tot <- sampleAnn_tot[id_samples,]
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

sampleGroups = list()
comparisons = list()
sampleNames = NULL
sampleGroups[["all"]] = sampleIds
sampleGroups[["Case"]] = sampleAnn$Temp_ID[sampleAnn$Dx==1] 
sampleGroups[["Ref"]] = sampleAnn$Temp_ID[sampleAnn$Dx==0] 
comparisons[[1]] = c("pairwise","all","Ref")

curComparison = comparisons[[1]]
print(paste("Processing: ",paste(curComparison,collapse="::"),sep=""))
allSamp = sampleGroups[[curComparison[[2]]]] # all
reference = sampleGroups[[curComparison[[3]]]] # ref (controls)
overlap = intersect(allSamp,reference) # ref
nSel=floor(length(overlap)*perc_comp) #floor(length(overlap)/nFolds)

# not balanced 
ovMat=sapply(seq(1,nFolds),function(X){
  
  vec=rep(F,length(overlap))
  set.seed(42+X)
  vec[sample.int(length(overlap),nSel)]=T
  return(vec)
  
})


rownames(ovMat)=overlap

# extract outside ref and comparison
excluded <- lapply(1:ncol(ovMat), function(x) rownames(ovMat)[ovMat[,x]])
curCases <- lapply(excluded, function(x) union(x,setdiff(allSamp,overlap))) # union part of controls and cases
curReference <- lapply(excluded, function(x) setdiff(reference,x))

# find controls never used as reference
controls_not_ref <- setdiff(reference, unique(unlist(curReference)))
# randomly select cases in the same size of controls

if(length(controls_not_ref)>thr_n){
  set.seed(seed_fixed - 10)
  controls_not_ref <- sample(controls_not_ref, size = thr_n, replace = F)
}

if(length(controls_not_ref)>=length(which(sampleAnn$Dx==1))){
 
  N_male <- length(which(sampleAnn$Gender[sampleAnn$Dx == 1] == 0))
  N_female <- length(which(sampleAnn$Gender[sampleAnn$Dx == 1] == 1))
  set.seed(seed_fixed*10)
  sel_controls_male <- sample(controls_not_ref[controls_not_ref %in% sampleAnn$Temp_ID[sampleAnn$Gender == 0]], size = N_male, replace = F)  
  set.seed(seed_fixed*10+20)
  sel_controls_female <- sample(controls_not_ref[controls_not_ref %in% sampleAnn$Temp_ID[sampleAnn$Gender == 1]], size = N_female, replace = F)
  sel_controls <- c(sel_controls_male, sel_controls_female)
  sel_cases <-sampleAnn$Temp_ID[sampleAnn$Dx==1]

}else{
  
  # match cases by gender
  N_male <- length(which(sampleAnn$Gender[sampleAnn$Individual_ID %in% controls_not_ref] == 0))
  N_female <- length(which(sampleAnn$Gender[sampleAnn$Individual_ID %in% controls_not_ref] == 1))

  if(length(which(sampleAnn$Dx==1 & sampleAnn$Gender == 0))<=N_male*0.9 | length(which(sampleAnn$Dx==1 & sampleAnn$Gender == 1))<=N_female*0.9){
    set.seed(seed_fixed)
    sel_cases <- sample(sampleAnn$Temp_ID[sampleAnn$Dx==1], size = length(controls_not_ref), replace = F)
  }else{
    set.seed(seed_fixed)
    sel_cases_male <- sample(sampleAnn$Temp_ID[sampleAnn$Dx==1 & sampleAnn$Gender == 0], size = N_male, replace = F)
    set.seed(seed_fixed+1)
    sel_cases_female <- sample(sampleAnn$Temp_ID[sampleAnn$Dx==1 & sampleAnn$Gender == 1], size = N_female, replace = F)
    sel_cases <- c(sel_cases_female, sel_cases_male)
  }
  sel_controls <- controls_not_ref
}

# write new covariate file (base data and target data)

if(length(controls_not_ref)<length(which(sampleAnn$Dx==1))){
  cov_base <- sampleAnn[!sampleAnn$Temp_ID %in% c(sel_controls, sel_cases), ]
  cov_base <- cov_base[, !colnames(cov_base) %in% c('Temp_ID', 'split_id')]
  # save
  write.table(sprintf('%scovariateMatrix_base_seed%i.txt', outFold, seed_fixed), x = cov_base, col.names = T, row.names = F, sep = '\t', quote = F)
}

cov_target <- sampleAnn[sampleAnn$Temp_ID %in% c(sel_controls, sel_cases), ]
cov_target <- cov_target[, !colnames(cov_target) %in% c('Temp_ID', 'split_id')]
#save
write.table(sprintf('%scovariateMatrix_target_seed%i.txt', outFold, seed_fixed), x = cov_target, col.names = T, row.names = F, sep = '\t', quote = F)

# split into external data
set.seed(9)
id_ext <- sample(1:nrow(cov_target), size = round(nrow(cov_target)*perc_ext),replace = F)
write.table(sprintf('%scovariateMatrix_target_seed%i_externalData.txt', outFold, seed_fixed), x = cov_target[id_ext,], col.names = T, row.names = F, sep = '\t', quote = F)
write.table(sprintf('%scovariateMatrix_target_seed%i_trainData.txt', outFold, seed_fixed), x = cov_target[-id_ext,],col.names = T, row.names = F, sep = '\t', quote = F)


