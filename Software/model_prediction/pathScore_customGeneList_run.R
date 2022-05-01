#!/usr/bin/env Rscript
# load specific pathway object and compute pathway score for them

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(qvalue))

parser <- ArgumentParser(description="Compute PathCcore for a specific pathway structure")

parser$add_argument("--pathwayStruct_file", type = "character", help = "List of genes (.RData)")
parser$add_argument("--tscore_file", type = "character", help = "already computed tscore")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--geneSetName", type = "character", help = "name pathway custom")
parser$add_argument("--abs_tscore", type = "logical", default = F, help = "if T also the pathway using absolute values of tscore is computed [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
pathwayStruct_file <- args$pathwayStruct_file
tscore_file <- args$tscore_file
sampleAnn_file <- args$sampleAnn_file
geneSetName <- args$geneSetName
abs_tscore <- args$abs_tscore
outFold <- args$outFold

# ####################################################################
# pathwayStruct_file <- '/home/luciat/eQTL_PROJECT/refData/CMC_GeneSets_Hypothesis-driven-for-Enrichement.RData'
# tscore_file <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/scz_aber_eur/devgeno0.01_testdevgeno0/predictedTscores.txt'
# sampleAnn_file <- '/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/scz_aber_eur.covariateMatrix.txt'
# outFold <- '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/scz_aber_eur/devgeno0.01_testdevgeno0/'
# geneSetName <- 'CMC_GeneSets_'
# ####################################################################

#########################
### function needed #####
#########################
# modified version of PGSEA that returns the mean t-score of the genes in each gene set
modPGSEA = function (exprs, geneSets, range=c(1,Inf), center=F, p.value=0.005, method="mean") {
  # initialize and prepare results matrix containing mean t-scores for each gene set and sample
  #if (is(exprs, "ExpressionSet")) exprs <- exprs(exprs)
  if (!is.list(geneSets)) stop("geneSets need to be a list")
  if (center) exprs <- scale(exprs, scale = FALSE)
  results <- matrix(NA, length(geneSets), ncol(exprs))
  rownames(results) <- names(geneSets)
  colnames(results) <- colnames(exprs)
  mode(results) <- "numeric"
  if (is.logical(p.value)) p.results <- results
  for (i in 1:length(geneSets)) {
    # determine the rows in the t-score table (exprs) that belong to the current gene set
    if (class(geneSets[[i]]) == "smc") {
      geneSetIds <- geneSets[[i]]@ids
    } else if (class(geneSets[[i]]) %in% c("GeneColorSet", "GeneSet")) {
      geneSetIds <- geneSets[[i]]@geneIds
    } else {
      geneSetIds <- unique(geneSets[[i]]$geneIds)
    }
    if (options()$verbose) cat("Testing region ", i, "\n")
    ix <- match(geneSetIds, rownames(exprs))
    # skip gene sets that are not sufficiently present in the t-score table (exprs)
    ix <- ix[!is.na(ix)]
    present <- sum(!is.na(ix))
    if (present < range[1]) {
      if (options()$verbose) cat("Skipping gene set ", i, " (too small): ", present, ",\n")
      next
    }
    if (present > range[2]) {
      if (options()$verbose) cat("Skipping gene set ", i, " (too large): ", present, "\n")
      next
    }
    # calculate an aggregate of the t-scores for the genes in the current gene set
    if (method == "mean") curFunction = mean
    if (method == "median") curFunction = median
    if (method == "abs_mean"){curFunction = function(x){mean(abs(x))}}
    selected = exprs[ix,,drop=F]; if (!is.matrix(selected)) selected <- as.matrix(selected)
    #unselected = exprs[-ix,,drop=F]; if (!is.matrix(unselected)) unselected <- as.matrix(unselected) 
    mVal = try(apply(selected,2,curFunction)) # do not scale or normalize because the t-scores have an absolute interpretation
    if (nrow(selected)<5) {
      pVal = rep(1,ncol(selected))
    } else {
      pVal = as.numeric(try(apply(selected,2,function(X) { unlist(t.test(X))[["p.value"]] })))
    }
    stat = list()
    for (q in 1:length(mVal)) {
      # use the absolute difference (mVal) rather than the z score as score variable
      stat[[q]] <- list(statistic = mVal[q], p.value = pVal[q])
    }
    names(stat) <- names(mVal)          
    # transfer score values into results matrix
    if (is.list(stat)) {
      ps <- unlist(lapply(stat, function(x) x$p.value))
      stat <- unlist(lapply(stat, function(x) x$statistic))
      if (!is.na(p.value)) {
        if (is.numeric(p.value)) {
          stat[ps > p.value] <- NA
        }
        else {
          p.results[i, ] <- ps
        }
      }
    }
    results[i, ] <- as.numeric(stat)
    # set all values of gene sets that do not fulfill the range requirements to NA
    for (w in 1:ncol(selected)) {
      if (sum(!is.na(selected[,w]))<range[1] | sum(!is.na(selected[,w])) > range[2]) results[i, w] <- NA
    }
  }
  if (is.logical(p.value) & !is.na(p.value)) {
    return(list(results = results, p.results = p.results))
  } else {
    return(results)
  }
}


##### load data #####
sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors = F, sep = '\t') 
pathwayCustom <- get(load(pathwayStruct_file))

tscoreMat <- read.table(tscore_file, h=T, sep = '\t', check.names = F)
geneID <- tscoreMat[,1]
tscoreMat <- tscoreMat[,-1]
sampleID <- sapply(colnames(tscoreMat), function(x) strsplit(x, split = '.vs')[[1]][1])
sampleID <- unname(sampleID)
colnames(tscoreMat) <- sampleID
# remove NA samples
sampleID_new <- colnames(tscoreMat)[colSums(is.na(tscoreMat))==0]
tscoreMat <- tscoreMat[,colSums(is.na(tscoreMat))==0]
sampleAnn <- sampleAnn[sampleAnn$Individual_ID %in% sampleID_new,]

geneSetNames <- sapply(pathwayCustom, function(x) x$name)
exprs <- as.matrix(tscoreMat)
rownames(exprs) <- geneID
results <- modPGSEA(exprs = exprs, geneSets = pathwayCustom, center=F, p.value=T) 

scores = results[["results"]]; dimnames(scores) = list(geneSetNames,gsub("\\.t","",dimnames(scores)[[2]]))
pvalues = results[["p.results"]]; dimnames(pvalues) = list(geneSetNames,paste("pval",gsub("\\.t",".pval",dimnames(pvalues)[[2]]),sep="_"))
valid = apply(scores,1,function(X) { sum(is.na(X)) == 0 })
insufficientCoverage = paste(geneSetNames[!valid])
scores = scores[valid,]
pvalues = pvalues[valid,]

output = data.frame(rownames(scores),scores)
names(output) = c("Name",colnames(scores))
output2 = data.frame(rownames(scores),pvalues)
names(output2) = c("Name",colnames(pvalues))
write.table(output,sprintf("%sPathway_%s_scores.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)     
write.table(output2,sprintf("%sPathway_%s_pval.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)

identical(colnames(scores), sampleAnn$Individual_ID)# same position due to zscore computation

sigScores = scores; #sigScores[pvalues>0.005] = NA
if(any(sampleAnn$Dx == 1)){
  
  cases=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==1,1]) ## mod LT
  pvalCombined = apply(pvalues[,cases],1,function(X) { pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE) })
  # save
  write.table(x = data.frame(chisq = sort(pvalCombined), pathway = names(sort(pvalCombined))), 
              file = sprintf('%sPathwaySummaryPvalues_cases_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
}

controls=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==0,1]) ## mod LT
pvalCombined_controls = apply(pvalues[,controls],1,function(X) {pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE)})
# save
write.table(x = data.frame(chisq = sort(pvalCombined_controls), pathway = names(sort(pvalCombined_controls))), 
            file = sprintf('%sPathwaySummaryPvalues_controls_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)


if(abs_tscore){
  # compute pathscore using absolute value of tscore
  results = modPGSEA(exprs = exprs, geneSets = pathwayCustom, center=F, p.value=T, method = 'abs_mean')
  
  scores = results[["results"]]; dimnames(scores) = list(geneSetNames,gsub("\\.t","",dimnames(scores)[[2]]))
  pvalues = results[["p.results"]]; dimnames(pvalues) = list(geneSetNames,paste("pval",gsub("\\.t",".pval",dimnames(pvalues)[[2]]),sep="_"))
  valid = apply(scores,1,function(X) { sum(is.na(X)) == 0 })
  insufficientCoverage = paste(geneSetNames[!valid])
  scores = scores[valid,]
  pvalues = pvalues[valid,]
  
  output = data.frame(rownames(scores),scores)
  names(output) = c("Name",colnames(scores))
  output2 = data.frame(rownames(scores),pvalues)
  names(output2) = c("Name",colnames(pvalues))
  write.table(output,sprintf("%sPathway_%s_scores_abstscore.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)
  write.table(output2,sprintf("%sPathway_%s_pval_abstscore.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)
  
  identical(colnames(scores), sampleAnn$Individual_ID)# same position due to zscore computation
  
  sigScores = scores; #sigScores[pvalues>0.005] = NA
  if(any(sampleAnn$Dx == 1)){
    
    cases=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==1,1]) ## mod LT
    pvalCombined = apply(pvalues[,cases],1,function(X) { pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE) })
    # save
    write.table(x = data.frame(chisq = sort(pvalCombined), pathway = names(sort(pvalCombined))),
                file = sprintf('%sPathwaySummaryPvalues_cases_%s_abstscore.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
  }
  
  controls=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==0,1]) ## mod LT
  pvalCombined_controls = apply(pvalues[,controls],1,function(X) {pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE)})
  # save
  write.table(x = data.frame(chisq = sort(pvalCombined_controls), pathway = names(sort(pvalCombined_controls))),
              file = sprintf('%sPathwaySummaryPvalues_controls_%s_abstscore.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
  
  
}



