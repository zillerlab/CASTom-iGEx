# from gwas results, create list of genes in significant loci
# match with PriLer genes list

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library('argparse'))
suppressPackageStartupMessages(library('LDlinkR'))
suppressPackageStartupMessages(library('biomaRt'))
suppressPackageStartupMessages(library('proxysnps'))

parser <- ArgumentParser(description="extract loci GWAS and match with PriLer results")
parser$add_argument("--token_number", type = "integer", help = "token for LDlinkR")
parser$add_argument("--GWAS_res", type = "character", help = "GWAS results file, must contain all chromosomes and be zipped")
parser$add_argument("--windows_bp", type = "integer", default = 1000000, help = "windows to find independent association")
parser$add_argument("--pval_corr_thr", type = "double", default = 0.05, help = "correct for pvalues thr")
parser$add_argument("--outFold", type = "character", help = "")
parser$add_argument("--refFold", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--tscore_discovery_file", type = "character", help = "file with PriLer result for discovery dataset")
parser$add_argument("--tscore_replication_file", type = "character",default = 'NA', help = "file with PriLer result for replication dataset")

args <- parser$parse_args()
token_number <- args$token_number
GWAS_res <- args$GWAS_res
windows_bp <- args$windows_bp
pval_corr_thr <- args$pval_corr_thr
outFold <- args$outFold
refFold <- args$refFold
tscore_discovery_file <- args$tscore_discovery_file
tscore_replication_file <- args$tscore_replication_file
pheno_name <- args$pheno_name

#####################################################################################################
# token_number <- '943676852910'
# GWAS_res <- '/psycl/g/mpsziller/lucia/refData/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
# windows_bp <- 1000000
# pval_corr_thr <- 0.05
# pheno_name <- 'CAD'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/'
# refFold <- '/psycl/g/mpsziller/lucia/priler_project/refData/'
# tscore_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/tscore_pval_CAD_HARD_covCorr.txt'
#####################################################################################################

HLA_reg <- c(26000000, 34000000)

gwas_res <- read.table(gzfile(GWAS_res), header = T, stringsAsFactors = F)
gwas_res <- gwas_res[!is.na(gwas_res$chr), ]
# correct pvalues
gwas_res$pval_corr <- p.adjust(gwas_res$p.value_gc, method = 'BH')
# order based on position and chr
gwas_res <- gwas_res[order(gwas_res$bp_hg19), ]
gwas_res <- gwas_res[order(gwas_res$chr), ]

# filter out
gwas_res_sign <- gwas_res[gwas_res$pval_corr<=pval_corr_thr, ]

df_loci <- data.frame(chr = NULL, start = NULL, end = NULL, n_sign_SNPs = NULL, best_SNP = NULL, best_SNP_pval = NULL, best_SNP_BHpval = NULL, stringsAsFactors = F)

# filter out for MHC locus, include only 1 snp (if any)
sign_MHC <- gwas_res_sign[gwas_res_sign$chr == 6, ]
sign_MHC <- sign_MHC[sign_MHC$bp_hg19 <= HLA_reg[2] & sign_MHC$bp_hg19 >= HLA_reg[1],] 

df_ind_loci <- rbind(df_loci, data.frame(chr = 6, start = HLA_reg[1], end = HLA_reg[2], n_sign_SNPs = nrow(sign_MHC), 
                                     best_SNP = sign_MHC$Markername[which.min(sign_MHC$p.value_gc)], best_SNP_pval = sign_MHC$p.value_gc[which.min(sign_MHC$p.value_gc)], 
                                     best_SNP_BHpval = sign_MHC$pval_corr[which.min(sign_MHC$p.value_gc)], stringsAsFactors = F))

tmp <- list()
# exclude MHC region
gwas_res_sign <- gwas_res_sign[!gwas_res_sign$Markername %in% sign_MHC$Markername, ]
snps_pred <- NULL

for(i in 1:nrow(gwas_res_sign)){
  
  print(i)
  
  tmp_snp <- gwas_res_sign[i, ]
  adj_snps <- gwas_res_sign[gwas_res_sign$chr == tmp_snp$chr &  abs(tmp_snp$bp_hg19 - gwas_res_sign$bp_hg19) <= windows_bp, ]
  
  if(nrow(adj_snps)==1){
    tmp[[i]] <- adj_snps
    
  }else{
  
    snp_name <- adj_snps$snptestid
    not_rs_ss <- sapply(snp_name, function(x) substr(x, start = 1, stop = 2) != 'rs' & substr(x, start = 1, stop = 2) != 'ss')
    ss_name <- sapply(snp_name, function(x) substr(x, start = 1, stop = 2) == 'ss')
    snp_name[ss_name] <- paste0('chr', adj_snps$chr, ':', adj_snps$bp_hg19)
    
    # remove not rs or ss snps
    snp_name <- snp_name[!not_rs_ss]
    adj_snps <- adj_snps[!not_rs_ss,]
    print(length(snp_name))
    
    if(!identical(snp_name, snps_pred)){
      
      ld_mat <- LDmatrix(snp_name, pop = "EUR", r2d = "r2", token = token_number)
      
      if(!identical(as.character(ld_mat[,1]), snp_name)){
        if(!any(ss_name)){
          print('not all snps found, remove not matched')
          adj_snps <- adj_snps[match(ld_mat$RS_number,adj_snps$snptestid), ]
        }else{
          print(paste0('name changed (ss) but same length: ', length(snp_name) == nrow(ld_mat)))
        }
      }
      
      ld_mat <- ld_mat[,-1]
      ld_mat[is.na(ld_mat)] <- 0
      
      if(any(ld_mat[upper.tri(ld_mat)] > 0.1)){
        locus_list <- apply(ld_mat, 1, function(x) abs(x)>0.1)
        len_snp <- c()
        keep_snp <- c()
        for(j in 1:nrow(locus_list)){
          tmp_sel <-  adj_snps[locus_list[j,],]
          tmp_sel <- tmp_sel[!tmp_sel$Markername %in% len_snp, ]
          len_snp <- unique(c(len_snp, tmp_sel$Markername))
          keep_snp <- unique(c(keep_snp, tmp_sel$Markername[which.min(tmp_sel$p.value_gc)]))
        }
        #print(keep_gene)
        tmp[[i]] <- adj_snps[adj_snps$Markername %in% keep_snp, ]
      }
    }
    
    snps_pred <- snp_name
    
  }
}  

LDind_snps <- do.call(rbind, tmp)
LDind_snps <- LDind_snps[!duplicated(LDind_snps$Markername),]
  
# investigate single indipendent snps to define physical regions
for(i in 1:nrow(LDind_snps)){
  
  print(i)
  if(substr(LDind_snps$snptestid[i], start = 1, stop = 2) == 'ss'){
    LDind_snps$snptestid[i] <-  paste0('chr', LDind_snps$chr[i], ':', LDind_snps$bp_hg19[i])
  }
  
  try(proxy_tmp <- LDproxy(LDind_snps$snptestid[i], pop = "EUR", r2d = "r2", token = token_number))
  
  if(!grepl('error',as.character(proxy_tmp[1,1]))){
  
    tmp <- sapply(as.character(proxy_tmp$Coord), function(x) strsplit(x, split = 'chr')[[1]][2])
    proxy_tmp$chr <- sapply(tmp, function(x) strsplit(x, split = '[:]')[[1]][1])
    proxy_tmp$pos <- sapply(tmp, function(x) strsplit(x, split = '[:]')[[1]][2])
    proxy_tmp_red <- proxy_tmp[proxy_tmp$R2 >= 0.6, ]
  
    print(paste0('n. snps not in LD: ', nrow(proxy_tmp)-nrow(proxy_tmp_red), '/', nrow(proxy_tmp)))
  
    df_ind_loci <- rbind(df_ind_loci,data.frame(chr = LDind_snps$chr[i], start = min(proxy_tmp_red$pos), end = max(proxy_tmp_red$pos), 
                                              n_sign_SNPs = sum(gwas_res_sign$chr == LDind_snps$chr[i] & gwas_res_sign$bp_hg19 <= max(proxy_tmp_red$pos) & gwas_res_sign$bp_hg19 >= min(proxy_tmp_red$pos)), 
                                              best_SNP = LDind_snps$Markername[i], best_SNP_pval = LDind_snps$p.value_gc[i], 
                                              best_SNP_BHpval = LDind_snps$pval_corr[i], stringsAsFactors = F))

  }
}

df_ind_loci$start <- as.numeric(df_ind_loci$start)
df_ind_loci$end <- as.numeric(df_ind_loci$end)
df_ind_loci$size <- df_ind_loci$end - df_ind_loci$start

# put together loci close by 250kb
comb_res <- vector(mode = 'list', length = 22)
new_loci <- df_ind_loci

for(i in 1:22){
  
  print(i)
  tmp_res <- new_loci[new_loci$chr == i, ]
  comb_res[[i]] <- tmp_res
  
  # repeat iteratively
  while(any(tmp_res$start[-1] - tmp_res$end[-nrow(tmp_res)] < 0 | 
            abs(tmp_res$start[-1] - tmp_res$end[-nrow(tmp_res)]) <=250000)){
    
    tmp <- tmp_res[tmp_res$chr == i, ]
    if(nrow(tmp) == 1){
      comb_res[[i]] <- data.frame(chr = tmp$chr, start = tmp$start, end = tmp$end, size = tmp$size, n_sign_SNPs = tmp$n_sign_SNPs, 
                                best_SNP = tmp$best_SNP, best_SNP_pval = tmp$best_SNP_pval, best_SNP_BHpval = tmp$best_SNP_BHpval, stringsAsFactors = F)
    
    }else{
    
      comb_res[[i]] <- data.frame(chr = NULL, start = NULL, end = NULL, size = NULL, n_sign_SNPs = NULL, best_SNP = NULL, best_SNP_pval = NULL, best_SNP_BHpval = NULL, stringsAsFactors = F)
      while(nrow(tmp)>0){
        new <- tmp[c(T, tmp$start[-1] - tmp$end[1]<=250000), ]
        comb_res[[i]] <- rbind(comb_res[[i]], data.frame(chr = new$chr[1], start = min(new$start), end = max(new$end), size = max(new$end) - min(new$start),
                                                           n_sign_SNPs = sum(new$n_sign_SNPs), 
                                                           best_SNP = new$best_SNP[which.min(new$best_SNP_pval)], best_SNP_pval = new$best_SNP_pval[which.min(new$best_SNP_pval)],
                                                           best_SNP_BHpval = new$best_SNP_BHpval[which.min(new$best_SNP_pval)], stringsAsFactors = F))
        tmp <- tmp[!tmp$best_SNP %in% new$best_SNP, ]
      }
    }
    tmp_res <- comb_res[[i]]
  }
  
}

df_loci <- do.call(rbind, comb_res) 
df_loci <- df_loci[order(df_loci$start),]
df_loci <- df_loci[order(df_loci$chr),]

# save:
write.table(df_loci, file = sprintf('%sGWAS_%s_loci.txt', outFold, pheno_name), col.names = T, row.names = F, quote = F, sep = '\t')

# get genes in the loci, find genes st +/-200kb TSS intersect with the region
# already annotated to build for PriLer (load)
biomart_annTSS <- read.table(sprintf("%s/hg19.ENSEMBL_geneTSS_biomart_correct.txt", refFold), h=T, stringsAsFactors = F, sep = '\t')
biomart_gene <- read.table(sprintf("%s/hg19.ENSEMBL_genes_biomart.txt", refFold), h=T, stringsAsFactors = F, sep = '\t')

# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# filterlist <- paste(df_loci$chr, df_loci$start, df_loci$end, sep = ':')
# gene_int <- getBM(attributes = c( "ensembl_gene_id","external_gene_name","chromosome_name", "start_position", "end_position","gene_biotype"),
#                  filters = c("chromosomal_region"),values = list(chromosomal_region=filterlist), 
#                  mart = ensembl)

biomart_annTSS <- biomart_annTSS[order(biomart_annTSS$chromstart),]
biomart_annTSS$chr_id <- sapply(biomart_annTSS$chrom, function(x) strsplit(x, split = 'chr')[[1]][2])
biomart_annTSS$chr_id <- as.numeric(biomart_annTSS$chr_id)
biomart_annTSS <- biomart_annTSS[!is.na(biomart_annTSS$chr_id),]
biomart_annTSS <- biomart_annTSS[order(biomart_annTSS$chr_id),]
biomart_annTSS$start_reg <- biomart_annTSS$chromstart - 200000
biomart_annTSS$end_reg <- biomart_annTSS$chromstart + 200000
biomart_annTSS$start_reg[biomart_annTSS$start_reg<0] <- 0

# annotate loci
df_loci$gene_int <- NA 
df_loci$ngenes <- NA 

for(i in 1:nrow(df_loci)){
  
  start_loci <- df_loci$start[i]
  end_loci <- df_loci$end[i]
  id1 <- biomart_annTSS$chr_id == df_loci$chr[i] & start_loci <= biomart_annTSS$end_reg & start_loci >= biomart_annTSS$start_reg
  id2 <- biomart_annTSS$chr_id == df_loci$chr[i] & end_loci <= biomart_annTSS$end_reg & end_loci >= biomart_annTSS$start_reg
  id3 <- biomart_annTSS$chr_id == df_loci$chr[i] & end_loci >= biomart_annTSS$end_reg & start_loci <= biomart_annTSS$start_reg
  gene_tmp <- biomart_annTSS[id1 | id2 | id3, ]
  df_loci$gene_int[i] <- paste0(gene_tmp$external_gene_name, collapse = ',')
  df_loci$ngenes[i] <- nrow(gene_tmp)
}

# save
write.table(df_loci, file = sprintf('%sGWAS_%s_loci_annotatedGenes.txt', outFold, pheno_name), col.names = T, row.names = F, 
            quote = F, sep = '\t')

tot_genes <- unlist(lapply(df_loci$gene_int, function(x) strsplit(x, split = '[,]')[[1]]))
tot_genes <- unique(tot_genes)

# load PriLer result and intersect genes: new association?
tscore_tot <- read.table(tscore_discovery_file, h=T, stringsAsFactors=F)
tscore_sign <- tscore_tot[tscore_tot[, 10] <= pval_corr_thr,]
genes <- unique(tscore_sign$external_gene_name)
new_genes <- genes[!genes %in% tot_genes]
tscore_new <- tscore_sign[tscore_sign$external_gene_name %in% new_genes, ]

new_genes_annotation <- biomart_annTSS[biomart_annTSS$external_gene_name %in% new_genes, c('chrom', 'chromend', 'ensembl_gene_id', 'external_gene_name') ]
new_genes_annotation <- cbind(new_genes_annotation, 
                        biomart_gene[match(new_genes_annotation$ensembl_gene_id, biomart_gene$ensembl_gene_id), c('start_position', 'end_position')])
colnames(new_genes_annotation)[colnames(new_genes_annotation) == 'chromend'] <- 'TSS_start'

# merge loci
new_list <- list()
for(i in 1:22){
  print(i)
  
  tmp <- new_genes_annotation[new_genes_annotation$chrom == paste0('chr', i), ]
  new_list[[i]] <- data.frame(chrom = NULL, 
  ensembl_gene_id = NULL, external_gene_name = NULL, 
  start_position= NULL, end_position = NULL, ngenes = NULL)
  
  if(nrow(tmp)>1){
    dist_mat <- sapply(tmp$start_position, function(x) abs(x - tmp$start_position) < 250000)
    colnames(dist_mat) <- tmp$external_gene_name
    tmp_name <- unique(apply(dist_mat, 1, function(x) paste0(colnames(dist_mat)[x], collapse = ',')))
    for(j in 1:length(tmp_name)){
      gene_id <- strsplit(tmp_name[j], split = '[,]')[[1]]
      new_list[[i]] <- rbind(new_list[[i]], 
                             data.frame(chrom = paste0('chr', i), 
                             ensembl_gene_id = paste0(tmp$ensembl_gene_id[match(gene_id, tmp$external_gene_name)], collapse = ','),
                             external_gene_name = tmp_name[j], 
                             start_position= min(tmp$start_position[tmp$external_gene_name %in% gene_id]), 
                             end_position= max(tmp$end_position[tmp$external_gene_name %in% gene_id]), 
                             ngenes = length(gene_id), stringsAsFactors = F))  
      
      
    }
  }else{
    if(nrow(tmp)==1){new_list[[i]] <- cbind(tmp[, !colnames(tmp) %in% 'TSS_start'], data.frame(ngenes = 1))}
  }
}

new_list <- do.call(rbind,new_list)

# merge position that share genes
new_list_genes <- lapply(new_list$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
update_list <- new_list
tmp <- new_list

while(nrow(tmp)>0){
  
  print(nrow(tmp))
  tmp_row <- tmp[1, ]
  gene_tmp <- strsplit(tmp_row$external_gene_name, split = '[,]')[[1]]
  keep <- sapply(new_list_genes, function(x) any(sapply(gene_tmp, function(y) y %in% x)))
  if(any(keep[-1])){
    tmp_list <- tmp[keep, ]
    new <- data.frame(chrom = tmp_list$chrom[1], 
                      ensembl_gene_id = paste0(unique(unlist(lapply(tmp_list$ensembl_gene_id, function(x) strsplit(x , split = '[,]')[[1]]))), collapse = ','), 
                      external_gene_name = paste0(unique(unlist(lapply(tmp_list$external_gene_name, function(x) strsplit(x , split = '[,]')[[1]]))), collapse = ','), 
                      start_position = min(tmp_list$start_position), end_position = max(tmp_list$end_position), 
                      ngenes = length(unique(unlist(lapply(tmp_list$ensembl_gene_id, function(x) strsplit(x , split = '[,]')[[1]])))))
    
    update_list <- update_list[!update_list$external_gene_name %in% tmp_list$external_gene_name,]
    update_list <- rbind(update_list, new)
  }
  tmp <- tmp[-c(1:sum(keep)),]
  new_list_genes <- lapply(tmp$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
  
}

update_list <- update_list[order(update_list$start_position), ]
update_list <- update_list[order(as.numeric(sapply(update_list$chrom, 
                                        function(x) strsplit(x, split = 'chr')[[1]][2]))), ]
new_list <- update_list

# repeat  
new_list_genes <- lapply(new_list$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
update_list <- new_list
tmp <- new_list

while(nrow(tmp)>0){
  
  print(nrow(tmp))
  tmp_row <- tmp[1, ]
  gene_tmp <- strsplit(tmp_row$external_gene_name, split = '[,]')[[1]]
  keep <- sapply(new_list_genes, function(x) any(sapply(gene_tmp, function(y) y %in% x)))
  if(any(keep[-1])){
    tmp_list <- tmp[keep, ]
    new <- data.frame(chrom = tmp_list$chrom[1], 
                      ensembl_gene_id = paste0(unique(unlist(lapply(tmp_list$ensembl_gene_id, function(x) strsplit(x , split = '[,]')[[1]]))), collapse = ','), 
                      external_gene_name = paste0(unique(unlist(lapply(tmp_list$external_gene_name, function(x) strsplit(x , split = '[,]')[[1]]))), collapse = ','), 
                      start_position = min(tmp_list$start_position), end_position = max(tmp_list$end_position), 
                      ngenes = length(unique(unlist(lapply(tmp_list$ensembl_gene_id, function(x) strsplit(x , split = '[,]')[[1]])))))
    
    update_list <- update_list[!update_list$external_gene_name %in% tmp_list$external_gene_name,]
    update_list <- rbind(update_list, new)
  }
  tmp <- tmp[-c(1:sum(keep)),]
  new_list_genes <- lapply(tmp$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
  
}

update_list <- update_list[order(update_list$start_position), ]
update_list <- update_list[order(as.numeric(sapply(update_list$chrom, 
                                                   function(x) strsplit(x, split = 'chr')[[1]][2]))), ]

# check GWAS for new loci:
update_list$best_GWAS_pvalue <- NA
update_list$best_GWAS_sign <- NA
update_list$best_GWAS_Markername <- NA

for(i in 1:nrow(update_list)){
  print(i)
  tmp <- gwas_res[gwas_res$chr == as.numeric(strsplit(update_list$chrom[i], split = 'chr')[[1]][2])
                  & gwas_res$bp_hg19 <= update_list$end_position[i] + 500000 & gwas_res$bp_hg19 >= update_list$start_position[i] - 500000, ]
  
  update_list$best_GWAS_pvalue[i] <- min(tmp$p.value_gc)
  update_list$best_GWAS_Markername[i] <- tmp$Markername[which.min(tmp$p.value_gc)]
  update_list$best_GWAS_sign[i] <- tmp$pval_corr[which.min(tmp$p.value_gc)] <= pval_corr_thr
}

# save:
write.table(update_list, file = sprintf('%snewloci_Priler_annotatedGenes_%s.txt', outFold, pheno_name), 
            col.names = T, row.names = F, sep = '\t', quote = F)


# are the new loci reproduced on discovery?
# combine dataset correctly
list_genes <- unique(unlist(lapply(update_list$external_gene_name[!update_list$best_GWAS_sign], function(x) strsplit(x, split = '[,]')[[1]])))
if(tscore_replication_file != 'NA'){
  
  tscore_rep <- read.table(tscore_replication_file, h=T, stringsAsFactors=F)
  tscore_rep$id <- paste0(tscore_rep[,1], '_tissue_', tscore_rep$tissue)
  
  tscore_sign_red <- tscore_sign[tscore_sign$external_gene_name %in% list_genes, ]
  tscore_sign_red$id <- paste0(tscore_sign_red[,1], '_tissue_', tscore_sign_red$tissue)
  
  tscore_rep <- tscore_rep[match(tscore_sign_red$id,tscore_rep$id),]
  tscore_sign_red$rep_zstat <- tscore_rep[,7]
  tscore_sign_red$rep <- F
  tscore_sign_red$rep[tscore_sign_red$rep_zstat*tscore_sign_red[,7]>0 & tscore_rep[,8]<=0.05] <- T
  
}else{
  
  tscore_sign_red <- tscore_sign[tscore_sign$external_gene_name %in% list_genes, ]
  tscore_sign_red$id <- paste0(tscore_sign_red[,1], '_tissue_', tscore_sign_red$tissue)
  
}

# save 
write.table(tscore_sign_red, file = sprintf('%snewloci_Priler_tscore_%s.txt', outFold, pheno_name), 
            col.names = T, row.names = F, sep = '\t', quote = F)



