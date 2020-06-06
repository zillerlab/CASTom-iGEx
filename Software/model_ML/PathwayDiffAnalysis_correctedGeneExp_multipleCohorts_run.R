#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# using corrected gene expression, compute T-scores
# reference dataset externally build

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(matrixStats))

parser <- ArgumentParser(description="compute T-scores and Pathway Scores multiple cohorts combined")

parser$add_argument("--input_file", type = "character", help = "total path corrected predicted gene expression (gz file)")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info to use as reference (only Dx == 0)")
parser$add_argument("--sampleAnn_ref_file", type = "character", nargs = '*', help = "1 file for each cohort, samples to be used as reference, must include Dx == 0")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "vector with name of the cohorts, order must match sampleAnn_ref_file")
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
#parser$add_argument("--nFolds", type="integer", default = 20, help = "number of fold for control partition")
parser$add_argument("--outFold", type="character", help = "Output folder")


args <- parser$parse_args()
input_file <- args$input_file
sampleAnn_file <- args$sampleAnn_file
sampleAnn_ref_file <- args$sampleAnn_ref_file
reactome_file <- args$reactome_file
GOterms_file <- args$GOterms_file
name_cohorts <- args$name_cohorts
#nFolds <- args$nFolds
outFold <- args$outFold

# ####################################################################
# input_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/devgeno0.01_testdevgeno0/corrected_predictedExpression.txt.gz'
# GOterms_file <- '/mnt/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/mnt/lucia/refData/ReactomePathways.gmt'
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/Meta_Analysis_CAD/covariateMatrix.txt'
# name_cohorts <- c('CG', 'German1',  'German2' ,  'German3',  'German4',  'German5', 'LURIC', 'MG', 'WTCCC')
# sampleAnn_ref_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/', cohorts, '/covariateMatrix_base_seed1235.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/devgeno0.01_testdevgeno0/'
# ####################################################################

#################################################################
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

# function to compute t statistic
t_stat <- function(x){
  n <- ncol(x)
  mu <- rowMeans(x)
  sigma <- rowSds(x)
  return(mu/(sigma/sqrt(n)))
}

#################################################################

### load data ###
sampleAnn <- read.table(sampleAnn_file, header = T, stringsAsFactors = F)
expDat <- read.delim(gzfile(input_file),sep="\t",header=T,  check.names = F, stringsAsFactors=F)
id_samples <- match(sampleAnn$Individual_ID, colnames(expDat))

eMat <- as.matrix(expDat[,id_samples])
rownames(eMat) <- expDat$external_gene_name
geneInfo <- expDat[, -id_samples]
identical(sampleAnn$Individual_ID, colnames(eMat)) # same order

# change Individual_ID name if start with a number or '-' present
sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) >= 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Individual_ID[id_n], function(x) paste0('X', x))   
}
colnames(eMat) <- sampleAnn$Temp_ID

# add Dx info if not present
if(!'Dx' %in% colnames(sampleAnn)){sampleAnn$Dx <- 0}

# load ref dataset
sampleAnn$ref <- F
for(i in 1:length(sampleAnn_ref_file)){
  tmp <- read.table(sampleAnn_ref_file[i], header = T, stringsAsFactors = F)
  sampleAnn$ref[sampleAnn$Individual_ID %in% paste(name_cohorts[i], tmp$Individual_ID, sep = '_')] <- T
}

# compute T-scores, extract reference dataset
inputData <- eMat
sampleGroups <- list()
comparisons <- list()

curType <- "corrected_gene_expression"
sampleIds <- colnames(inputData)
ind <- rownames(inputData)==""
inputData <- inputData[!ind,]
geneNames <- rownames(inputData)
sampleGroups[["all"]] = sampleIds
sampleGroups[["Ref"]] = sampleAnn$Temp_ID[sampleAnn$ref & sampleAnn$Dx == 0] # controls that are in ref group
sampleGroups[["Case"]] = sampleAnn$Temp_ID[ (!sampleAnn$ref) | (sampleAnn$ref & sampleAnn$Dx == 1)] # all the samples that are not in ref + cases in ref 
comparisons[[1]] = c("pairwise","all","Ref")

eMat <- eMat[rowSums(is.na(eMat))==0,]

tscoreTables <- list()
curComparison <- comparisons[[1]]
print(paste("Processing: ",paste(curComparison,collapse="::"),sep=""))
allSamp <- sampleGroups[[curComparison[[2]]]] # all
reference <- sampleGroups[[curComparison[[3]]]] # ref (controls)
overlap <- intersect(allSamp,reference) # ref
nSel <- floor(length(overlap)*0.2) #floor(length(overlap)/nFolds)

values <- inputData
row.names(values) <-  geneNames
tTable <-  data.frame(geneId = geneNames)
tempTable <- tTable

# # not balanced 
# ovMat=sapply(seq(1,nFolds),function(X){
#   
#   vec=rep(F,length(overlap))
#   set.seed(42+X)
#   vec[sample.int(length(overlap),nSel)]=T
#   return(vec)
#   
# })
# 
# rownames(ovMat) <- overlap # ref division (folds)
# 
# # determine current subset of samples excluding one of the samples that overlap between allSamp and reference
# for (i in 1:nFolds){
#   
#   excluded <- rownames(ovMat)[ovMat[,i]]
#   comparisonLabel <- paste("Reference_excluding_fold",i,sep="_")#paste("Reference_excluding",excluded,sep="_")
#   print(comparisonLabel)
#   curCases <- union(excluded,setdiff(allSamp,overlap)) # union part of controls and cases
#   curReference <- setdiff(reference,excluded) # the remaining controls
#   curIds <- union(curCases,curReference) # total
#   curValues <- values[,curIds]
#   
#   # use T-statistics (add moderate option)
#   tempTable <- cbind(tempTable, sapply(curCases, function(x) t_stat(curValues[, colnames(curValues) == x] - curValues[,colnames(curValues) %in% curReference])))
#   
#   # # define design matrix for limma analysis
#   # designFactor <- factor(ifelse(curIds%in%curCases,curIds,"reference")) # remaining control = reference, otherwise ids
#   # design <- model.matrix(~ -1+designFactor)
#   # rownames(design) <- curIds
#   # colnames(design) <- gsub("designFactor","",colnames(design))
#   # fit <- lmFit(curValues, design)
#   # 
#   # contrastMatrixItems <- character(0)
#   # for (j in 1:length(curCases)) {
#   #   curItem <- paste(curCases[j],"-reference",sep="")
#   #   contrastMatrixItems <- c(contrastMatrixItems,curItem)
#   # }
#   # contrastMatrixString <- paste(contrastMatrixItems,collapse=", ")
#   # contrastMatrixCmd <- paste("makeContrasts(",contrastMatrixString,", levels=design)",sep="")
#   # contrast.matrix <- eval(parse(text=contrastMatrixCmd))
#   # contrastNames <- dimnames(contrast.matrix)[["Contrasts"]];# contrastNames
#   # # perform limma analysis and obtain t-scores
#   # fit2 <- contrasts.fit(fit, contrast.matrix) # compute estimated coefficients and standard errors for a given set of contrasts
#   # fit2 <- eBayes(fit2) # compute moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage
#   
#   # for (k in 1:length(contrastNames)) {
#   #   curContrast <- contrastNames[k];
#   #   curLabel <- paste(comparisonLabel," :: ",curContrast)
#   #   topGenes <- topTable(fit2,coef=curContrast,sort="none",n=Inf)[,"t",drop=F]
#   #   tempTable <- cbind(tempTable,topGenes)
#   #   names(tempTable)[ncol(tempTable)] <- paste(curLabel,"::",names(topGenes)[ncol(topGenes)])
#   # }
# }
# 
# # calculate the median t-score from the different leave-one-out reference sets
# for (curSample in allSamp) {
#   curCols <- curSample
#   # curCols <- names(tempTable)
#   # curCols <- curCols[grep(paste(curSample,"- reference"),curCols)]
#   for (suffix in c("t")) {
#     # selectedCols <- curCols[grep(paste("::",suffix),curCols)]
#     newCol <- paste(curSample,"vs reference","::",suffix)
#     print(newCol)
#     tTable[,newCol] <- apply(tempTable[, colnames(tempTable) == curCols, drop=F],1,mean)
#     # tTable[,newCol] <- apply(tempTable[,selectedCols,drop=F],1,mean)
#   }
# }
# tscoreTables[[paste(curComparison,collapse="::")]] <- tTable


curCases = sampleGroups[["Case"]] # target samples
curReference = sampleGroups[["Ref"]] # base samples
curIds = union(curCases,curReference) # total
curValues = values[,curIds]
comparisonLabel = "Reference_baseSamples"
print(comparisonLabel)

tempTable <- cbind(tempTable, sapply(curCases, function(x) t_stat(curValues[, colnames(curValues) == x] - curValues[,colnames(curValues) %in% curReference])))

newCol = paste(curCases,"vs reference","::","t")
colnames(tempTable)[-1] <- newCol

tscoreTables[[paste(curComparison,collapse="::")]] = tempTable
curCases_new <- sampleAnn$Individual_ID[match(curCases, sampleAnn$Temp_ID)]

# save
# tscoreTables[[1]] <- tscoreTables[[1]][is.element(tscoreTables[[1]][,1], geneInfo$external_gene_name),] 
# tmp <- sapply(colnames(tscoreTables[[1]])[-1], function(x) paste(strsplit(x,' ')[[1]][-1], collapse = " "))
# colnames(tscoreTables[[1]])[-1] <- mapply(function(x,y) paste(x, y), x = sampleAnn$Individual_ID, y=tmp)

# save
tscoreTables[[1]]=tscoreTables[[1]][is.element(tscoreTables[[1]][,1], geneInfo$external_gene_name),] 
tmp <- sapply(colnames(tscoreTables[[1]])[-1], function(x) paste(strsplit(x,' ')[[1]][-1], collapse = " "))
colnames(tscoreTables[[1]])[-1] <- mapply(function(x,y) paste(x, y), x = curCases_new, y=tmp)
write.table(tscoreTables[[1]],sprintf("%scorrected_predictedTscores.txt", outFold),sep="\t",row.names=F,quote=F)

##########################
### pathScore reactome ###
##########################

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
names(output) = c("sampleName",colnames(scores))
output2 = data.frame(rownames(scores),pvalues)
names(output2) = c("sampleName",colnames(pvalues))
write.table(output,sprintf("%scorrected_Pathway_%s_scores.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)     
write.table(output2,sprintf("%scorrected_Pathway_%s_pval.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)      
scorecards[[curComparison]] = output

sigScores = scores; #sigScores[pvalues>0.005] = NA
if(any(sampleAnn$Dx == 1)){
  
  cases=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==1,1]) ## mod LT
  pvalCombined = apply(pvalues[,cases],1,function(X) { pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE) })
  # save
  write.table(x = data.frame(chisq = sort(pvalCombined), pathway = names(sort(pvalCombined))), 
              file = sprintf('%scorrected_PathwaySummaryPvalues_cases_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
}

controls=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==0,1]) ## mod LT
pvalCombined_controls = apply(pvalues[,controls],1,function(X) {pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE)})
# save
write.table(x = data.frame(chisq = sort(pvalCombined_controls), pathway = names(sort(pvalCombined_controls))), 
            file = sprintf('%scorrected_PathwaySummaryPvalues_controls_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)

if(any(sampleAnn$Dx == 1)){
  
  thres=0.05
  rAll=rowSums(pvalues<=thres)
  rCase=rowSums(pvalues[,cases]<=thres)
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%scorrected_PathwayEnrichment_%s_FracCases.txt", outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
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
              sprintf("%scorrected_PathwayEnrichment_%s_FracCases_pvalAdj.txt",outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  
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
names(output) = c("sampleName",colnames(scores))
output2 = data.frame(rownames(scores),pvalues)
names(output2) = c("sampleName",colnames(pvalues))
write.table(output,sprintf("%scorrected_Pathway_%s_scores.txt", outFold, geneSetName),sep="\t",quote=F,row.names=F)     
write.table(output2,sprintf("%scorrected_Pathway_%s_pval.txt", outFold,geneSetName),sep="\t",quote=F,row.names=F)      
scorecards[[curComparison]] = output

sigScores = scores; #sigScores[pvalues>0.005] = NA
if(any(sampleAnn$Dx == 1)){
  
  cases=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==1,1]) ## mod LT
  pvalCombined = apply(pvalues[,cases],1,function(X) { pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE) })
  # save
  write.table(x = data.frame(chisq = sort(pvalCombined), pathway = names(sort(pvalCombined))), 
              file = sprintf('%scorrected_PathwaySummaryPvalues_cases_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)
}

controls=is.element(colnames(scores),sampleAnn[sampleAnn$Dx==0,1]) ## mod LT
pvalCombined_controls = apply(pvalues[,controls],1,function(X) {pchisq(-2*sum(log(X)),2*length(X),lower.tail=FALSE)})
# save
write.table(x = data.frame(chisq = sort(pvalCombined_controls), pathway = names(sort(pvalCombined_controls))), 
            file = sprintf('%scorrected_PathwaySummaryPvalues_controls_%s.txt', outFold, geneSetName), sep = '\t', quote = F, col.names = T, row.names = F)

if(any(sampleAnn$Dx == 1)){
  
  thres=0.05
  rAll=rowSums(pvalues<=thres)
  rCase=rowSums(pvalues[,cases]<=thres)
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%scorrected_PathwayEnrichment_%s_FracCases.txt", outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  padj=apply(pvalues,2,function(X){
    return(qvalue(X)$qvalue)
  })
  
  #find gene sets that show high deviation
  thres=0.05
  rAll=rowSums(padj<=thres)
  rCase=rowSums(padj[,cases]<=thres)
  
  rFrac=rCase/rAll
  ind=!is.nan(rFrac)
  pEn=cbind(rAll,rCase,rFrac)
  
  write.table(data.frame(cbind(pathway=rownames(pEn[ind,]), pEn[ind,])),
              sprintf("%scorrected_PathwayEnrichment_%s_FracCases_pvalAdj.txt",outFold, geneSetName),sep="\t",quote=F, row.names = F)
  
  
}


