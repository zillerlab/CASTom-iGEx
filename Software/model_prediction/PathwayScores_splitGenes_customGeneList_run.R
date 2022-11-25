#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(parallel))


parser <- ArgumentParser(description="Differential pathway analysis, to be used with a high number of samples > 10000")

parser$add_argument("--split_tot", type="integer", default = 0, help = "total number of partitions")
parser$add_argument("--input_file", type = "character", help = "output file predicted Tscore, if split_tot not 0 refers to common part")
parser$add_argument("--pathwayStruct_file", type = "character", help = "custom pathway structure .RData")
parser$add_argument("--covDat_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--geneSetName", type = "character", help = "name pathway custom")
parser$add_argument("--ncores", type="integer", default = 10, help = "n cores for the parallelization")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
input_file <- args$input_file
pathwayStruct_file <- args$pathwayStruct_file
split_tot <- args$split_tot
geneSetName <- args$geneSetName
covDat_file <- args$covDat_file
ncores <- args$ncores
outFold <- args$outFold

# ####################################################################
# input_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/predictedTscores_splitGenes'
# split_tot <- 100
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/covariateMatrix_latestW.txt'
# pathwayStruct_file <- '/psycl/g/mpsziller/lucia/castom-igex/refData/WikiPathways_2019_Human.RData'
# ncores <- 10
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# geneSetName <- 'WikiPath2019Human'
####################################################################

#########################
### function needed #####
#########################
# modified version of PGSEA that returns the mean t-score of the genes in each gene set
# consider only 1 gene set, parallelized
modPGSEA_bigmat = function(geneSet, range_gene=c(1,Inf), type) {
  
  
  # determine the rows in the t-score table (exprs) that belong to the current gene set
  if (type == "Reactome") {
    geneSetIds <- geneSet@ids
  }else{
    geneSetIds <- unique(geneSet$geneIds)
  }
  ix <- match(geneSetIds, genesID)
  # skip gene sets that are not sufficiently present in the t-score table (exprs)
  ix <- ix[!is.na(ix)]
  present <- sum(!is.na(ix))
  
  if(present < range_gene[1] | present > range_gene[2]){
    pVal = rep(NA,nrow(exprs))
    mVal = rep(NA,nrow(exprs))
    
  }else{
    
    # calculate an aggregate of the t-scores for the genes in the current gene set
    
    selected = exprs[,ix,drop=F]; if (!is.matrix(selected)) selected <- as.matrix(selected)
    #unselected = exprs[-ix,,drop=F]; if (!is.matrix(unselected)) unselected <- as.matrix(unselected) 
    if(ncol(selected) == 1){
      mVal <- as.vector(selected)
    }else{
      mVal <- rowMeans(selected) # do not scale or normalize because the t-scores have an absolute interpretation 
    }
    
    if(ncol(selected)<5){
      pVal = rep(1,nrow(selected))
    }else{
      pVal = as.numeric(apply(selected, 1, function(X){ unlist(t.test(X))[["p.value"]] }))
    }
    
  }
  
  return(list(results = mVal, p.results = pVal))
  
}


### load data ###
sampleAnn <- read.table(covDat_file, header = T, stringsAsFactors = F)
sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) == 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Individual_ID[id_n], function(x) paste0('X', x))   
}
if(!'Dx' %in% colnames(sampleAnn)){sampleAnn$Dx <- 0}

if(split_tot>0){
  
  tscoreMat_list <- list()
  samplesID <- list()
  genesID <- NULL
  
  for(i in 1:split_tot){
    print(i)
    if(file.exists(sprintf('%s%i.RData', input_file, i))){
      tmp <- get(load(sprintf('%s%i.RData', input_file, i)))
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
  
  try(if(any(!sapply(samplesID, function(x) all(sampleAnn$Temp_ID %in% x) & all(x %in% sampleAnn$Temp_ID)))) stop("ERROR: different samples in the splitted genes"))
  
  if(!identical(sampleAnn$Temp_ID, samplesID[[1]])){
    id_samples <- match(samplesID[[1]], sampleAnn$Temp_ID)
    sampleAnn <- sampleAnn[id_samples,]
  }
  
  # samplesID <- samplesID[[1]]
  
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
  
}else{
  
  tmp <- get(load(input_file))
  genesID <-tmp$geneId
  samplesID <- colnames(tmp)[-1]
  try(if(!(all(sampleAnn$Temp_ID %in% samplesID) & (samplesID %in% sampleAnn$Temp_ID))) stop("ERROR: different samples in the splitted genes"))
  
  if(!identical(sampleAnn$Temp_ID, samplesID)){
    id_samples <- match(samplesID, sampleAnn$Temp_ID)
    sampleAnn <- sampleAnn[id_samples,]
  }
  
  tmp <- as.matrix(t(tmp[, -1]))
  tscoreMat <- big.matrix(ncol=ncol(tmp), nrow = nrow(tmp), type = "double", init = 0, dimnames = NULL, shared = F)
  tscoreMat[,] <- tmp
}


# remove sample that have NAs
id_s <- rowSums(is.na(tscoreMat[,])) == 0
sampleAnn <- sampleAnn[id_s,]
samplesID_new <- sampleAnn$Temp_ID
exprs <- tscoreMat[id_s, ]

#######################################################
#################### custom pathway ###################
#######################################################

pathwayCustom <- get(load(pathwayStruct_file))

minGenesPerGeneSet = 1
maxGenesPerGeneSet = 50000

# consider only pathway that include at least 1 gene
id_path <- which(sapply(pathwayCustom, function(x) any(genesID %in% x$geneIds)))
pathwayCustom <- pathwayCustom[id_path]

geneSetNames <- sapply(pathwayCustom, function(x) x$name)
geneSetLengths <- sapply(pathwayCustom,function(X) { length(unique(X$geneIds))})
curGeneSets <- pathwayCustom

# parallelize per element in the pathway
cl <- makeCluster(ncores)
clusterExport(cl, "exprs")
clusterExport(cl, "genesID")
clusterExport(cl, "modPGSEA_bigmat")
clusterExport(cl, "curGeneSets")
clusterExport(cl, "minGenesPerGeneSet")
clusterExport(cl, "maxGenesPerGeneSet")
  
results <- parLapply(cl, 1:length(curGeneSets), function(x) 
  modPGSEA_bigmat(geneSet = curGeneSets[[x]], 
                  range_gene=c(minGenesPerGeneSet, maxGenesPerGeneSet), 
                  type = 'custom'))
stopCluster(cl)
scores <- t(sapply(results, function(x) x$results))               
pvalues <- t(sapply(results, function(x) x$p.results))                 
rm(results)
  
valid = apply(scores, 1, function(X) { sum(is.na(X)) == 0 })
scores = scores[valid,]
pvalues = pvalues[valid,]
  
scores = data.frame(pathID = geneSetNames[valid], scores)
colnames(scores)[-1] = samplesID_new
pvalues = data.frame(pathID = geneSetNames[valid], pvalues)
colnames(pvalues)[-1] = samplesID_new
  
save(scores, file = sprintf("%sPathway_%s_scores.RData", outFold, geneSetName))
save(pvalues, file = sprintf("%sPathway_%s_pval.RData", outFold, geneSetName))
  
if(any(sampleAnn$Dx == 1)){
    
  cases=is.element(colnames(scores),sampleAnn$Temp_ID[sampleAnn$Dx==1]) ## mod LT
  pvalCombined = apply(pvalues[,cases],1,function(X) { pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE) })
  # save
  df <- data.frame(chisq = sort(pvalCombined), pathway = pvalues$pathID[order(pvalCombined)])
  # save(df, file = sprintf('%s/PathwaySummaryPvalues_cases_%s.RData', outFold, geneSetName))
  write.table(x = df, file = sprintf('%sPathwaySummaryPvalues_cases_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
}
  
controls=is.element(colnames(scores),sampleAnn$Temp_ID[sampleAnn$Dx==0]) ## mod LT
pvalCombined_controls = apply(pvalues[,controls],1,function(X) {pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE)})
# save
df <- data.frame(chisq = sort(pvalCombined_controls), pathway = pvalues$pathID[order(pvalCombined_controls)])
# save
write.table(x = df, file = sprintf('%sPathwaySummaryPvalues_controls_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
  
if(any(sampleAnn$Dx == 1)){
    
    thres=0.05
    rAll=rowSums(pvalues[, -1]<=thres)
    rCase=rowSums(pvalues[,cases]<=thres)
    rFrac=rCase/rAll
    ind=!is.nan(rFrac)
    pEn=cbind(rAll,rCase,rFrac)
    
    write.table(data.frame(cbind(pathway=pvalues$pathID[ind], pEn[ind,])),
                sprintf("%sPathwayEnrichment_%s_FracCases.txt", outFold, geneSetName),sep="\t",quote=F, row.names = F)
    
    padj=apply(pvalues[, -1],2,function(X){
      return(qvalue(X)$qvalue)
    })
    
    #find gene sets that show high deviation
    thres=0.05
    rAll=rowSums(padj<=thres)
    rCase=rowSums(padj[,cases[-1]]<=thres)
    # rAll=rowSums(pvalues<=thres)
    # rCase=rowSums(pvalues[,cases]<=thres)
    rFrac=rCase/rAll
    ind=!is.nan(rFrac)
    pEn=cbind(rAll,rCase,rFrac)
    
    write.table(data.frame(cbind(pathway=pvalues$pathID[ind], pEn[ind,])),
                sprintf("%sPathwayEnrichment_%s_FracCases_pvalAdj.txt",outFold, geneSetName),sep="\t",quote=F, row.names = F)
    
    
}



