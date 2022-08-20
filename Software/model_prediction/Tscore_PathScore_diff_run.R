#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(qvalue))

parser <- ArgumentParser(description="Differential pathway analysis")

parser$add_argument("--input_file", type = "character", help = "output path prediction folder (or cov corrected data)")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--originalRNA", type="logical", default = F, help = "if true, original RNA without filtering used [default %(default)s]")
parser$add_argument("--thr_reliableGenes", type = "double", nargs = '*', default = c(0.01, 0), help = "threshold for reliable genes: dev_geno_tot and test_dev_geno [default %(default)s]")
parser$add_argument("--covDat_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--nFolds", type="integer", default = 20, help = "number of fold for control partition [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output folder")

args <- parser$parse_args()
input_file <- args$input_file
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
originalRNA <- args$originalRNA
covDat_file <- args$covDat_file
thr_reliableGenes <- args$thr_reliableGenes
nFolds <- args$nFolds
outFold <- args$outFold

# ####################################################################
# input_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/CG/predictedExpression.txt.gz'
# originalRNA <- F
# GOterms_file <- '/mnt/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/mnt/lucia/refData/ReactomePathways.gmt'
# thr_reliableGenes <- c(0.01, 0)
# nFolds <- 40
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/CG/covariateMatrix.txt'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/CG/devgeno0.01_testdevgeno0/'
# ####################################################################

#########################
### function needed #####
#########################
# modified version of PGSEA that returns the mean t-score of the genes in each gene set
modPGSEA = function (exprs, geneSets, range=c(1,Inf), center=F, p.value=0.005, method="mean") {
  # initialize and prepare results matrix containing mean t-scores for each gene set and sample
  if (is(exprs, "ExpressionSet")) exprs <- exprs(exprs)
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


### load data ###

if(originalRNA){
  
  sampleAnn <- read.table(covDat_file, header = T, stringsAsFactors = F, sep="\t")
  expDat <- read.table(input_file,sep="\t",header=T, check.names = F)
  
  id_samples <- match(sampleAnn$Individual_ID, colnames(expDat))
  eMat <- as.matrix(expDat[,id_samples])
  rownames(eMat) <- expDat$external_gene_name
  geneInfo <- expDat[, -id_samples]
  
  if(!identical(sampleAnn$Individual_ID, colnames(eMat))){print('ERROR: Annotation samples and expression not matching')} # same order
  
}else{
  
  sampleAnn <- read.table(covDat_file, header = T, stringsAsFactors = F,sep="\t")
  sampleAnn$Individual_ID <- as.character(sampleAnn$Individual_ID)
  expDat <- read.table(gzfile(input_file),sep="\t",header=T,  check.names = F)
  
  id_samples <- match(sampleAnn$Individual_ID, colnames(expDat))

  # filter genes
  expDat <- expDat[!(is.na(expDat$dev_geno) | is.na(expDat$test_dev_geno)), ]
  expDat <- expDat[expDat$dev_geno >= thr_reliableGenes[1], ]
  expDat <- expDat[expDat$test_dev_geno > thr_reliableGenes[2], ]
  
  eMat <- as.matrix(expDat[,id_samples])
  rownames(eMat) <- expDat$external_gene_name
  geneInfo <- expDat[, -id_samples]

  if(!identical(sampleAnn$Individual_ID, colnames(eMat))){stop('ERROR: Annotation samples and expression not matching')} 
  
}

# change Individual_ID name if start with a number or '-' or '*' present
sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) >= 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
id_a <-  sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '[*]')[[1]]) >= 2)
id_s <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = ' ')[[1]]) >= 2)
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Temp_ID[id_n], function(x) paste0('X', x))   
}
if(any(id_a)){
  sampleAnn$Temp_ID[id_a] <- sapply(sampleAnn$Temp_ID[id_a], function(x) paste0(strsplit(x, split = '[*]')[[1]], collapse = '_'))   
}
if(any(id_s)){
  sampleAnn$Temp_ID[id_s] <- sapply(sampleAnn$Temp_ID[id_s], function(x) paste0(strsplit(x, split = ' ')[[1]], collapse = '_'))   
}
colnames(eMat) <- sampleAnn$Temp_ID

# add Dx info if not present
if(!'Dx' %in% colnames(sampleAnn)){sampleAnn$Dx <- 0}

inputData <- eMat
sampleGroups <- list()
comparisons <- list()

curType <- "Gene_expression"
sampleIds <- colnames(inputData)
ind <- rownames(inputData)==""
inputData <- inputData[!ind,]
geneNames <- rownames(inputData)
sampleGroups[["all"]] = sampleIds
sampleGroups[["Case"]] = sampleAnn$Temp_ID[sampleAnn$Dx==1] 
sampleGroups[["Ref"]] = sampleAnn$Temp_ID[sampleAnn$Dx==0] 
comparisons[[1]] = c("pairwise","all","Ref")

# Perform Limma analysis to determine sample-specific T-scores relative to a reference population (for each gene)
eMat <- eMat[rowSums(is.na(eMat))==0,]

# produce a genes*samples table of t-scores for each comparison
for (i in 1:length(sampleGroups)){sampleGroups[[i]] = unique(sampleGroups[[i]])}

tscoreTables <- list()
curComparison <- comparisons[[1]]
print(paste("Processing: ",paste(curComparison,collapse="::"),sep=""))
allSamp <- sampleGroups[[curComparison[[2]]]] # all
reference <- sampleGroups[[curComparison[[3]]]] # ref (controls)
overlap <- intersect(allSamp,reference) # ref
nSel <- floor(length(overlap)*0.2) #floor(length(overlap)/nFolds)

values <- inputData
row.names(values) <- geneNames
tTable <- data.frame(geneId=geneNames)

# not balanced 
ovMat=sapply(seq(1,nFolds),function(X){
  vec=rep(F,length(overlap))
  set.seed(42+X)
  vec[sample.int(length(overlap),nSel)]=T
  return(vec)
})


rownames(ovMat) <- overlap # ref division (folds)
tempTable <- tTable

# determine current subset of samples excluding one of the samples that overlap between allSamp and reference
for (i in 1:nFolds){
  
  excluded <- rownames(ovMat)[ovMat[,i]]
  comparisonLabel <- paste("Reference_excluding_fold",i,sep="_")
  print(comparisonLabel)
  curCases <- union(excluded,setdiff(allSamp,overlap)) # union part of controls and cases
  curReference <- setdiff(reference,excluded) # the remaining controls
  curIds <- union(curCases,curReference) # total
  curValues <-values[,curIds]
  # define design matrix for limma analysis
  designFactor <- factor(ifelse(curIds%in%curCases,curIds,"reference")) # remaining control = reference, otherwise ids
  design <- model.matrix(~ -1+designFactor)
  rownames(design) <- curIds
  colnames(design) <- gsub("designFactor","",colnames(design))
  fit <- lmFit(curValues, design)
  contrastMatrixItems <- character(0)
  for (j in 1:length(curCases)) {
    curItem <- paste(curCases[j],"-reference",sep="")
    contrastMatrixItems <- c(contrastMatrixItems,curItem)
  }
  contrastMatrixString <- paste(contrastMatrixItems,collapse=", ")
  contrastMatrixCmd <- paste("makeContrasts(",contrastMatrixString,", levels=design)",sep="")
  contrast.matrix <- eval(parse(text=contrastMatrixCmd))
  contrastNames <- dimnames(contrast.matrix)[["Contrasts"]];# contrastNames
  # perform limma analysis and obtain t-scores
  fit2 <- contrasts.fit(fit, contrast.matrix) # compute estimated coefficients and standard errors for a given set of contrasts
  fit2 <- eBayes(fit2) # compute moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage
  for (k in 1:length(contrastNames)) {
    curContrast <- contrastNames[k];
    curLabel <- paste(comparisonLabel," :: ",curContrast)
    topGenes <- topTable(fit2,coef=curContrast,sort="none",n=Inf)[,"t",drop=F]
    tempTable <- cbind(tempTable,topGenes)
    names(tempTable)[ncol(tempTable)] <- paste(curLabel,"::",names(topGenes)[ncol(topGenes)])
  }
}

# calculate the median t-score from the different leave-one-out reference sets
for (curSample in allSamp) {
  curCols <- names(tempTable)
  curCols <- curCols[grep(paste(curSample,"- reference"),curCols)]
  for (suffix in c("t")) {
    selectedCols <- curCols[grep(paste("::",suffix),curCols)]
    newCol <- paste(curSample,"vs reference","::",suffix)
    print(newCol)
    tTable[,newCol] <- apply(tempTable[,selectedCols,drop=F],1,mean)
  }
}
tscoreTables[[paste(curComparison,collapse="::")]] <- tTable

# save
tscoreTables[[1]] <- tscoreTables[[1]][is.element(tscoreTables[[1]][,1], geneInfo$external_gene_name),] 

tmp <- tscoreTables[[1]]
colnames(tmp)[-1] <- sampleAnn$Individual_ID
write.table(tmp,sprintf("%spredictedTscores.txt", outFold),sep="\t",row.names=F,quote=F)

tmp <- sapply(colnames(tscoreTables[[1]])[-1], function(x) paste(strsplit(x,' ')[[1]][-1], collapse = " "))
colnames(tscoreTables[[1]])[-1] <- mapply(function(x,y) paste(x, y),x= sampleAnn$Individual_ID, y=tmp)

#################################################
#################### Reactome ###################
#################################################

gs <- readGmt(reactome_file)
geneSetName="Reactome"

gs=lapply(gs,function(X){
  X@ids=gsub(",1.0","",X@ids)
  return(X)
})

# lineage scorecard settings
minGenesPerGeneSet = 1
maxGenesPerGeneSet = 50000
# perform parametric gene set analysis on the moderated t statistics that compare each sample with the reference
performClustering = F
pvalScorecard = F
maxGeneSetsPlotted = 100
curGeneSets = gs #scorecardGeneSets
scorecards = list()

i=1
curComparison = names(tscoreTables)[i]
tTable = tscoreTables[[i]]
geneNames=tTable[,1]
print(paste("Processing:",curComparison))
# prepare data for PGSEA analysis
tCols = names(tTable)[grep(":: t",names(tTable))]
exprs = as.matrix(tTable[,tCols])
exprs=exprs[,colSums(is.nan(exprs))==0] # where the NaNs come from?
dimnames(exprs)[[1]] = geneNames
curColNames = sapply(dimnames(exprs)[[2]],function(X) { gsub(" vs.*$","",X)  }); names(curColNames) = NULL
dimnames(exprs)[[2]] = curColNames
geneSetNames = sapply(curGeneSets,function(X) { X@reference })
names(geneSetNames) = NULL    
geneSetLengths = sapply(curGeneSets,function(X) { length(X@ids) })
names(geneSetLengths) = NULL    
# perform PGSEA analysis on the t-scores relative to the reference and subsequently filter out those gene sets that do not fall within the size range
results = modPGSEA(exprs, curGeneSets, range=c(minGenesPerGeneSet, maxGenesPerGeneSet), center=F, p.value=T) 
scores = results[["results"]]; dimnames(scores) = list(geneSetNames,gsub("\\.t","",dimnames(scores)[[2]]))
pvalues = results[["p.results"]]; dimnames(pvalues) = list(geneSetNames,paste("pval",gsub("\\.t",".pval",dimnames(pvalues)[[2]]),sep="_"))
valid = apply(scores,1,function(X) { sum(is.na(X)) == 0 })
insufficientCoverage = paste(geneSetNames[!valid])
scores = scores[valid,]
pvalues = pvalues[valid,]
# for each lineage / group of gene sets, take the mean of the gene set scores of all contributing gene sets
# export lineage scorecard
output = data.frame(rownames(scores),scores)
names(output) = c("pathID",colnames(scores))
output2 = data.frame(rownames(scores),pvalues)
names(output2) = c("pathID",colnames(pvalues))
write.table(output,sprintf("%sPathway_%s_scores.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)     
write.table(output2,sprintf("%sPathway_%s_pval.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)      
scorecards[[curComparison]] = output

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




if(any(sampleAnn$Dx == 1)){
  
  thres=0.05
  rAll=rowSums(pvalues<=thres)
  rCase=rowSums(pvalues[,cases]<=thres)
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%sPathwayEnrichment_%s_FracCases.txt", outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  padj=apply(pvalues,2,function(X){
    return(qvalue(X)$qvalue)
  })
  
  #find gene sets that show high deviation
  thres=0.05
  rAll=rowSums(padj<=thres)
  rCase=rowSums(padj[,cases]<=thres)
  # rAll=rowSums(pvalues<=thres)
  # rCase=rowSums(pvalues[,cases]<=thres)
  
  
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%sPathwayEnrichment_%s_FracCases_pvalAdj.txt",outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  
}


################################################
#################### GO term ###################
################################################

# repeat the same for GO terms
go <- get(load(GOterms_file))
geneSetName <- 'GO'
# lineage scorecard settings
minGenesPerGeneSet = 1
maxGenesPerGeneSet = 50000
# perform parametric gene set analysis on the moderated t statistics that compare each sample with the reference
performClustering = F
pvalScorecard = F
maxGeneSetsPlotted = 100
#curGeneSets = allGeneSets
scorecards = list()

i=1
curComparison = names(tscoreTables)[i]
tTable = tscoreTables[[i]]
geneNames=tTable[,1]
print(paste("Processing:",curComparison))
# prepare data for PGSEA analysis
tCols = names(tTable)[grep(":: t",names(tTable))]
exprs = as.matrix(tTable[,tCols])
exprs=exprs[,colSums(is.nan(exprs))==0] # where the NaNs come from?
dimnames(exprs)[[1]] = geneNames
curColNames = sapply(dimnames(exprs)[[2]],function(X) { gsub(" vs.*$","",X)  }); names(curColNames) = NULL
dimnames(exprs)[[2]] = curColNames

geneSetNames = sapply(go,function(X) { X$GOID })
geneSetLengths = sapply(go,function(X) { length(unique(X$geneIds))})

# perform PGSEA analysis on the t-scores relative to the reference and subsequently filter out those gene sets that do not fall within the size range
results = modPGSEA(exprs, go, range=c(minGenesPerGeneSet, maxGenesPerGeneSet), center=F, p.value=T) 
scores = results[["results"]]; dimnames(scores) = list(geneSetNames,gsub("\\.t","",dimnames(scores)[[2]]))
pvalues = results[["p.results"]]; dimnames(pvalues) = list(geneSetNames,paste("pval",gsub("\\.t",".pval",dimnames(pvalues)[[2]]),sep="_"))
valid = apply(scores,1,function(X) { sum(is.na(X)) == 0 })
insufficientCoverage = paste(geneSetNames[!valid])
scores = scores[valid,]
pvalues = pvalues[valid,]
# for each lineage / group of gene sets, take the mean of the gene set scores of all contributing gene sets
# export lineage scorecard
output = data.frame(rownames(scores),scores)
names(output) = c("pathID",colnames(scores))
output2 = data.frame(rownames(scores),pvalues)
names(output2) = c("pathID",colnames(pvalues))
write.table(output,sprintf("%sPathway_%s_scores.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)     
write.table(output2,sprintf("%sPathway_%s_pval.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)      
scorecards[[curComparison]] = output

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



if(any(sampleAnn$Dx == 1)){
  
  thres=0.05
  rAll=rowSums(pvalues<=thres)
  rCase=rowSums(pvalues[,cases]<=thres)
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%sPathwayEnrichment_%s_FracCases.txt", outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  padj=apply(pvalues,2,function(X){
    return(qvalue(X)$qvalue)
  })
  
  #find gene sets that show high deviation
  thres=0.05
  rAll=rowSums(padj<=thres)
  rCase=rowSums(padj[,cases]<=thres)
  # rAll=rowSums(pvalues<=thres)
  # rCase=rowSums(pvalues[,cases]<=thres)
  
  
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%sPathwayEnrichment_%s_FracCases_pvalAdj.txt",outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  
}


