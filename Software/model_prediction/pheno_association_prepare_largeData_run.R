#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# consider results from big matrix, prepare for association, split pathway RData, create info files

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(PGSEA))

parser <- ArgumentParser(description="Gene and Pathway association analysis, prepare files")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--inputFold", type = "character", help = "Folder with results from Tscore/pathway scores computation")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data")
parser$add_argument("--split_tot", type = "integer",default = 100,  help = "number of split per genes")
parser$add_argument("--geneAnn_file", type = "character", help = "file with gene info from train, to be filtered")
parser$add_argument("--thr_reliableGenes", type = "double", nargs = '*', default = c(0.01, 0), help = "threshold for reliable genes: dev_geno_tot and test_dev_geno")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFold <- args$inputFold
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
sampleAnn_file <- args$sampleAnn_file
geneAnn_file <- args$geneAnn_file
thr_reliableGenes <- args$thr_reliableGenes
split_tot <- args$split_tot
outFold <- args$outFold

##########################################################################################
# inputFold <- 'OUTPUT_GTEx/predict_UKBB/Brain_Cortex/200kb/noGWAS/devgeno0.01_testdevgeno0/'
# sampleAnn_file <- 'INPUT_DATA/Covariates/covariatesMatrix.txt'
# split_tot <- 100
# geneAnn_file <- 'OUTPUT_GTEx/train_GTEx/Brain_Cortex/200kb/noGWAS/resPrior_regEval_allchr.txt'
# thr_reliableGenes <- c(0.01, 0)
# GOterms_file <- '/psycl/g/mpsziller/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/psycl/g/mpsziller/lucia/refData/ReactomePathways.gmt'
# outFold <- 'OUTPUT_GTEx/predict_UKBB/Brain_Cortex/200kb/noGWAS/devgeno0.01_testdevgeno0/'
##########################################################################################

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

sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) == 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Individual_ID[id_n], function(x) paste0('X', x))   
}


tscoreMat_list <- list()
samplesID <- list()
genesID <- NULL

for(i in 1:split_tot){
  print(i)
  if(file.exists(sprintf('%spredictedTscores_splitGenes%i.RData', inputFold, i))){
    tmp <- get(load(sprintf('%spredictedTscores_splitGenes%i.RData', inputFold, i)))
    genesID <- c(genesID, tmp$geneId)
    samplesID[[i]] <- colnames(tmp)[-1]
    tmp <- as.matrix(t(tmp[, -1]))
    tscoreMat_list[[i]] <- big.matrix(ncol=ncol(tmp), nrow = nrow(tmp), type = "double", init = 0, dimnames = NULL, shared = F)
    tscoreMat_list[[i]][,] <- tmp
  }else{
    print(paste('split',i,'does not exist'))
  }
}

split_tot <- length(tscoreMat_list)

# try(if(any(!sapply(samplesID, function(x) identical(sampleAnn$Temp_ID, x)))) stop("ERROR: different samples in the splitted genes"))
# samplesID <- samplesID[[1]]
try(if(any(!sapply(samplesID, function(x) all(sampleAnn$Temp_ID %in% x) & all(x %in% sampleAnn$Temp_ID)))) stop("ERROR: different samples in the splitted genes"))

if(!identical(sampleAnn$Temp_ID, samplesID[[1]])){
    id_samples <- match(samplesID[[1]], sampleAnn$Temp_ID)
    sampleAnn <- sampleAnn[id_samples,]
}

tscoreMat <- big.matrix(ncol=length(genesID), nrow = nrow(sampleAnn), type = "double", init = 0, dimnames = NULL, shared = F)
ngenes <- sapply(tscoreMat_list, ncol)
genes_id_start <- c(0,cumsum(ngenes[-1]))+1
genes_id_end <- cumsum(ngenes)
# correct for the last split (different dimensions than the others)
genes_id_start[split_tot] = genes_id_end[split_tot-1] +1
print(sum(genes_id_end-genes_id_start+1) == length(genesID))

for(i in 1:split_tot){
  print(i)
  tscoreMat[,genes_id_start[i]:genes_id_end[i]] <- tscoreMat_list[[i]][,]
}

rm(tscoreMat_list)

# remove sample that have NAs
id_s <- rowSums(is.na(tscoreMat[,])) == 0
sampleAnn <- sampleAnn[id_s,]
samplesID_new <- sampleAnn$Temp_ID
if(!all(id_s)){
  if(length(genesID)<=6000){
    tscoreMat <- as.big.matrix(tscoreMat[id_s, ], type = 'double', shared = F)    
  }else{
    old <- tscoreMat[id_s, ]
    tscoreMat <- big.matrix(ncol=length(genesID), nrow = length(which(id_s)), type = "double", init = 0, dimnames = NULL, shared = F)
    for(i in 1:length(genesID)){
      # print(i)
      tscoreMat[,i] <- old[,i]
    }
    rm(old)
  }
}  
  

# filter geneAnn
if(!identical(genesID, geneAnn$external_gene_name)){
  print('adjust genes annotation')
  id <- match(genesID,  geneAnn$external_gene_name)
  geneAnn <- geneAnn[id,]
}

df_tscore_info <- data.frame(ensembl_gene_id =  geneAnn$ensembl_gene_id, external_gene_name = geneAnn$external_gene_name,
                             dev_geno = geneAnn$dev_geno, test_dev_geno = geneAnn$test_dev_geno, stringsAsFactors = F)
save(df_tscore_info, file = sprintf('%stscore_info.RData', outFold))

print('Tscore mat preprocessed')

##################################
#### load pathScore Reactome #####
##################################
tmp <- get(load(sprintf('%sPathway_Reactome_scores.RData', inputFold)))
pathScoreID_reactome <- tmp[,1]
tmp <- tmp[,-1]
identical(colnames(tmp), samplesID_new) # same order samples
# assign to bigmatrix (less space)
pathScore_reactome <- big.matrix(ncol=length(pathScoreID_reactome), nrow = length(samplesID_new), type = "double", init = 0, dimnames = NULL, shared = F)
pathScore_reactome[,] <- as.matrix(t(tmp))
rm(tmp)
rm(scores)

# consider only pathaways that do not heave gene repetition, on pathwayScoreID add gene info
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
  pathScore_reactome <- as.big.matrix(pathScore_reactome[,-id_rm], type = 'double', shared = F)
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

save(df_pathR_info, file = sprintf('%spathScore_Reactome_info.RData', outFold))

# split path
step_path <- round(ncol(pathScore_reactome)/split_tot)
split_id = as.vector(sapply(1:split_tot, function(x) rep(x, step_path)))
if(length(split_id)<ncol(pathScore_reactome)){
  split_id <- c(split_id, rep(split_tot, ncol(pathScore_reactome) - length(split_id)))  
}else{
  if(length(split_id)>ncol(pathScore_reactome)){
    split_id <- split_id[-((ncol(pathScore_reactome)+1):length(split_id))]
  }
}

split_df <- data.frame(id = 1:ncol(pathScore_reactome), split = split_id)
for(i in 1:split_tot){
  print(i)
  path_id <- split_df$id[split_df$split == i]  
  scores <- as.data.frame(t(pathScore_reactome[, path_id]))
  colnames(scores) <- samplesID_new
  scores <- cbind(data.frame(path = pathScoreID_reactome[path_id], stringsAsFactors = F), scores)
  save(scores, file = sprintf('%sPathway_Reactome_scores_splitPath%i.RData', outFold, i))
  
}

print('pathScore reactome mat processed')


############################
#### load pathScore GO #####
############################

tmp <- get(load(sprintf('%sPathway_GO_scores.RData', inputFold)))
rm(scores)
pathScoreID_GO <- tmp[,1]
tmp <- tmp[,-1]
identical(colnames(tmp), samplesID_new) # same order samples
# assign to bigmatrix (less space)
pathScore_GO <- big.matrix(ncol=length(pathScoreID_GO), nrow = length(samplesID_new), type = "double", init = 0, dimnames = NULL, shared = F)
# problem with big vector
if(length(pathScoreID_GO)<=6000){
  pathScore_GO[,] <- as.matrix(t(tmp))  
}else{
  tmp <- as.matrix(t(tmp))
  for(i in 1:length(pathScoreID_GO)){
    # print(i)
    pathScore_GO[,i] <- tmp[,i]
  }
}

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
  tmp <- pathScore_GO[,-id_rm]
  pathScore_GO <- big.matrix(ncol=ncol(tmp), nrow = nrow(tmp), type = "double", init = 0, dimnames = NULL, shared = F)
  if(ncol(tmp)<=6000){
    pathScore_GO[,] <- tmp  
  }else{
    for(i in 1:ncol(tmp)){
      # print(i)
      pathScore_GO[,i] <- tmp[,i]
    }
  }
  rm(tmp)
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

save(df_pathGO_info, file = sprintf('%spathScore_GO_info.RData', outFold))

# split path
step_path <- round(ncol(pathScore_GO)/split_tot)
split_id = as.vector(sapply(1:split_tot, function(x) rep(x, step_path)))
if(length(split_id)<ncol(pathScore_GO)){
  split_id <- c(split_id, rep(split_tot, ncol(pathScore_GO) - length(split_id)))  
}else{
  if(length(split_id)>ncol(pathScore_GO)){
    split_id <- split_id[-((ncol(pathScore_GO)+1):length(split_id))]
  }
}

split_df <- data.frame(id = 1:ncol(pathScore_GO), split = split_id)
for(i in 1:split_tot){
  
  path_id <- split_df$id[split_df$split == i]  
  scores <- as.data.frame(t(pathScore_GO[, path_id]))
  colnames(scores) <- samplesID_new
  scores <- cbind(data.frame(path = pathScoreID_GO[path_id], stringsAsFactors = F), scores)
  save(scores, file = sprintf('%sPathway_GO_scores_splitPath%i.RData', outFold, i))
  
}

print('pathScore GO mat processed')



