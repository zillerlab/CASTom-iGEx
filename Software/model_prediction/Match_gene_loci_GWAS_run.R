# from gwas loci extract genes
# match with PriLer genes list

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library('argparse'))
# suppressPackageStartupMessages(library('LDlinkR'))
suppressPackageStartupMessages(library('biomaRt'))
suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('R.utils'))

parser <- ArgumentParser(description="match GWAS loci and match with PriLer results")
parser$add_argument("--loci_GWAS_res", type = "character", help = "processed GWAS loci results")
parser$add_argument("--GWAS_res", type = "character", help = "GWAS original results")
parser$add_argument("--pval_corr_thr", type = "double", default = 0.05, help = "correct for pvalues thr")
parser$add_argument("--outFold", type = "character", help = "")
parser$add_argument("--refFold", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--tscore_discovery_file", type = "character", help = "file with PriLer result for discovery dataset")
parser$add_argument("--tscore_replication_file", type = "character", default = NULL, help = "file with PriLer result for replication dataset")

args <- parser$parse_args()
loci_GWAS_res <- args$loci_GWAS_res
GWAS_res <- args$GWAS_res
pval_corr_thr <- args$pval_corr_thr
outFold <- args$outFold
refFold <- args$refFold
tscore_discovery_file <- args$tscore_discovery_file
tscore_replication_file <- args$tscore_replication_file
pheno_name <- args$pheno_name

#####################################################################################################
# loci_GWAS_res <- '/psycl/g/mpsziller/lucia/refData/108lociSCZ_SuppTable2_mergedLoci.txt'
# GWAS_res <- '/psycl/g/mpsziller/lucia/refData/Original_SCZ_variants_PGC.txt.gz'
# pval_corr_thr <- 0.05
# pheno_name <- 'SCZ'
# outFold <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/'
# refFold <- '/psycl/g/mpsziller/lucia/priler_project/refData/'
# tscore_discovery_file <- '/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/tscore_pval_Dx_pheno_covCorr_seltissues.txt'
# tscore_replication_file <- NULL
# ####################################################################################################

HLA_reg <- c(26000000, 34000000)
gwas_res <- fread(GWAS_res, header = T, stringsAsFactors = F, data.table = F)

if(pheno_name == 'SCZ'){
  colnames(gwas_res)[colnames(gwas_res) == 'CHR'] <- 'chr'
  colnames(gwas_res)[colnames(gwas_res) == 'P'] <- 'p.value_gc'
  colnames(gwas_res)[colnames(gwas_res) == 'BP'] <- 'bp_hg19'
  # colnames(gwas_res)[colnames(gwas_res) == 'SNP'] <- 'Markername'
  gwas_res$chr <- substr(gwas_res$chr, 4, 6)
  gwas_res$Markername <- paste0(gwas_res$chr, ':', gwas_res$bp_hg19, '_', gwas_res$A1, '_', gwas_res$A2)
}
if(grepl('CAD',pheno_name)){
  colnames(gwas_res)[colnames(gwas_res) == 'p-value_gc'] <- 'p.value_gc'
}

gwas_res <- gwas_res[!is.na(gwas_res$chr), ]
# correct pvalues
gwas_res$pval_corr <- p.adjust(gwas_res$p.value_gc, method = 'BH')

# load GWAS loci res and annotate with genes:
loci_GWAS <- read.table(loci_GWAS_res, h=T, stringsAsFactors = F, sep = '\t')

# get genes in the loci, find genes st +/-200kb TSS intersect with the region
# already annotated to build for PriLer (load)
biomart_annTSS <- read.table(sprintf("%s/hg19.ENSEMBL_geneTSS_biomart_correct.txt", refFold), h=T, stringsAsFactors = F, sep = '\t')
biomart_gene <- read.table(sprintf("%s/hg19.ENSEMBL_genes_biomart.txt", refFold), h=T, stringsAsFactors = F, sep = '\t')

biomart_annTSS <- biomart_annTSS[order(biomart_annTSS$chromstart),]
biomart_annTSS$chr_id <- sapply(biomart_annTSS$chrom, function(x) strsplit(x, split = 'chr')[[1]][2])
biomart_annTSS$chr_id <- as.numeric(biomart_annTSS$chr_id)
biomart_annTSS <- biomart_annTSS[!is.na(biomart_annTSS$chr_id),]
biomart_annTSS <- biomart_annTSS[order(biomart_annTSS$chr_id),]
biomart_annTSS$start_reg <- biomart_annTSS$chromstart - 200000
biomart_annTSS$end_reg <- biomart_annTSS$chromstart + 200000
biomart_annTSS$start_reg[biomart_annTSS$start_reg<0] <- 0

# annotate loci
loci_GWAS$gene_int <- NA 
loci_GWAS$ngenes <- NA 

for(i in 1:nrow(loci_GWAS)){
  
  start_loci <- loci_GWAS$start_reg[i]
  end_loci <- loci_GWAS$end_reg[i]
  id1 <- biomart_annTSS$chr_id == loci_GWAS$CHR[i] & start_loci <= biomart_annTSS$end_reg & start_loci >= biomart_annTSS$start_reg
  id2 <- biomart_annTSS$chr_id == loci_GWAS$CHR[i] & end_loci <= biomart_annTSS$end_reg & end_loci >= biomart_annTSS$start_reg
  id3 <- biomart_annTSS$chr_id == loci_GWAS$CHR[i] & end_loci >= biomart_annTSS$end_reg & start_loci <= biomart_annTSS$start_reg
  gene_tmp <- biomart_annTSS[id1 | id2 | id3, ]
  loci_GWAS$gene_int[i] <- paste0(gene_tmp$external_gene_name, collapse = ',')
  loci_GWAS$ngenes[i] <- nrow(gene_tmp)
}


# save
write.table(loci_GWAS, file = sprintf('%sGWAS_%s_loci_annotatedGenes.txt', outFold, pheno_name), col.names = T, row.names = F, 
            quote = F, sep = '\t')

tot_genes <- unlist(lapply(loci_GWAS$gene_int, function(x) strsplit(x, split = '[,]')[[1]]))
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
notmerged_list <- new_list

# merge position that share genes
stop_cond <- F
new_list_genes <- lapply(new_list$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
n_cycle <- 0

while(!stop_cond){
  
  n_cycle <- n_cycle +1
  print(paste0("##########", n_cycle, "##########"))
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
  new_list_genes <- lapply(new_list$external_gene_name, function(x) strsplit(x, split = '[,]')[[1]])
  
  matrix_int <- sapply(1:length(new_list_genes), function(x) any(new_list_genes[[x]] %in% 
                                                                   unlist(new_list_genes[-x])))
  stop_cond <- all(!matrix_int)
  
}

update_list <- update_list[order(update_list$start_position), ]
update_list <- update_list[order(as.numeric(sapply(update_list$chrom, 
                                                   function(x) strsplit(x, split = 'chr')[[1]][2]))), ]


# check GWAS for new loci:
update_list$best_GWAS_Markername <- NA
update_list$best_GWAS_pvalue <- NA
update_list$best_GWAS_sign <- NA
update_list$best_GWAS_signGW <- NA

for(i in 1:nrow(update_list)){
  print(i)
  tmp <- gwas_res[gwas_res$chr == as.numeric(strsplit(update_list$chrom[i], split = 'chr')[[1]][2])
                  & gwas_res$bp_hg19 <= update_list$end_position[i] + 500000 & gwas_res$bp_hg19 >= update_list$start_position[i] - 500000, ]
  
  update_list$best_GWAS_pvalue[i] <- min(tmp$p.value_gc)
  update_list$best_GWAS_Markername[i] <- tmp$Markername[which.min(tmp$p.value_gc)]
  update_list$best_GWAS_sign[i] <- tmp$pval_corr[which.min(tmp$p.value_gc)] <= pval_corr_thr
  update_list$best_GWAS_signGW[i] <- min(tmp$p.value_gc) <= 5*10^-8
  
}

# save:
write.table(update_list, file = sprintf('%snewloci_Priler_annotatedGenes_%s.txt', outFold, pheno_name), 
            col.names = T, row.names = F, sep = '\t', quote = F)


# are the new loci reproduced on discovery?
# combine dataset correctly
list_genes <- unique(unlist(lapply(update_list$external_gene_name[!update_list$best_GWAS_sign], function(x) strsplit(x, split = '[,]')[[1]])))
if(!is.null(tscore_replication_file)){
  
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



