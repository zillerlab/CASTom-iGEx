# combine all tissues results (CAD)
library(qvalue)
library(data.table)
library(rlist)

setwd('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/')
tissues_name <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage','Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
bp_loci <- 1000000
cis_size <- 200000

df_tscore_HARD <- df_pathR_HARD <- df_pathGO_HARD <- list()
df_tscore_SOFT <- df_pathR_SOFT <- df_pathGO_SOFT <- list()
df_pathWiki <- list()

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  
  tmp <- get(load(sprintf('predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData', t)))
  df_tscore_HARD[[i]] <- tmp$tscore[[1]]
  df_tscore_SOFT[[i]] <- tmp$tscore[[2]]
  df_tscore_HARD[[i]]$tissue <-  df_tscore_SOFT[[i]]$tissue <- t
  
  # add genes in the pathway and if there is an improvment in significance
  df_pathR_HARD[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR_SOFT[[i]] <- tmp$pathScore_reactome[[2]]
  df_pathR_HARD[[i]]$genes_path <- df_pathR_SOFT[[i]]$genes_path <- NA
  df_pathR_HARD[[i]]$improvement_sign <- df_pathR_SOFT[[i]]$improvement_sign <- NA
  for(j in 1:nrow(df_pathR_HARD[[i]])){
    df_pathR_HARD[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_reactome[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathR_SOFT[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_reactome[[2]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathR_HARD[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_reactome[[1]][[j]]$tscore[,8] > df_pathR_HARD[[i]][j,13])
    df_pathR_SOFT[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_reactome[[2]][[j]]$tscore[,8] > df_pathR_SOFT[[i]][j,13])
  }
  df_pathR_HARD[[i]]$tissue <-  df_pathR_SOFT[[i]]$tissue <- t
  
  df_pathGO_HARD[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO_SOFT[[i]] <- tmp$pathScore_GO[[2]]
  df_pathGO_HARD[[i]]$genes_path <- df_pathGO_SOFT[[i]]$genes_path <- NA
  df_pathGO_HARD[[i]]$improvement_sign <- df_pathGO_SOFT[[i]]$improvement_sign <- NA
  for(j in 1:nrow(df_pathGO_HARD[[i]])){
    df_pathGO_HARD[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_GO[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathGO_SOFT[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_GO[[2]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathGO_HARD[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_GO[[1]][[j]]$tscore[,8] > df_pathGO_HARD[[i]][j,15])
    df_pathGO_SOFT[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_GO[[2]][[j]]$tscore[,8] > df_pathGO_SOFT[[i]][j,15])
  }
  df_pathGO_HARD[[i]]$tissue <-  df_pathGO_SOFT[[i]]$tissue <- t
 
  tmp <- get(load(sprintf('predict_CAD/%s/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr_customPath_WikiPath2019Human.RData', t)))
  df_pathWiki[[i]] <- tmp$pathScore[[1]]
  df_pathWiki[[i]]$genes_path <- NA
  df_pathWiki[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathWiki[[i]])){
    df_pathWiki[[i]]$genes_path[j] <- paste0(tmp$info_pathScore[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathWiki[[i]]$improvement_sign[j] <- all(tmp$info_pathScore[[1]][[j]]$tscore[,8] > df_pathWiki[[i]][j,13])
  }
  df_pathWiki[[i]]$tissue <- t
   
}

df_tscore_HARD <- do.call(rbind, df_tscore_HARD)
df_tscore_SOFT <- do.call(rbind, df_tscore_SOFT)
df_pathR_HARD <- do.call(rbind, df_pathR_HARD)
df_pathR_SOFT <- do.call(rbind, df_pathR_SOFT)
df_pathGO_HARD <- do.call(rbind, df_pathGO_HARD)
df_pathGO_SOFT <- do.call(rbind, df_pathGO_SOFT)
df_pathWiki <- do.call(rbind, df_pathWiki)

# add overall corrected pval
df_tscore_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_tscore_HARD[, 8], method = 'BH')
df_tscore_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_tscore_SOFT[, 8], method = 'BH')
df_pathR_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_pathR_HARD[, 13], method = 'BH')
df_pathR_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_pathR_SOFT[, 13], method = 'BH')
df_pathGO_HARD$CAD_HARD_BHcorr_overall <- p.adjust(df_pathGO_HARD[, 15], method = 'BH')
df_pathGO_SOFT$CAD_SOFT_BHcorr_overall <- p.adjust(df_pathGO_SOFT[, 15], method = 'BH')
df_pathWiki$CAD_HARD_BHcorr_overall <- p.adjust(df_pathWiki[, 13], method = 'BH')


# create a function to remove pathway with 1 gene and recompute pvalues
recompte_path <- function(tissues_name, res, id_pval){
  tmp <- lapply(tissues_name, function(x) res[res$tissue == x & res$ngenes_tscore>1,])
  for(i in 1:length(tmp)){
    tmp[[i]][, id_pval+1] <- qvalue(tmp[[i]][, id_pval])$qvalue
    tmp[[i]][, id_pval+2] <- p.adjust(tmp[[i]][, id_pval], method = 'BH')
  }
  tmp <- do.call(rbind, tmp)
  tmp[, id_pval+6] <- p.adjust(tmp[, id_pval], method = 'BH')
  return(tmp)
}
df_pathR_HARD_red <- recompte_path(res = df_pathR_HARD, tissues_name = tissues_name, id_pval = 13)
df_pathR_SOFT_red <- recompte_path(res = df_pathR_SOFT, tissues_name = tissues_name, id_pval = 13)
df_pathGO_HARD_red <- recompte_path(res = df_pathGO_HARD, tissues_name = tissues_name, id_pval = 15)
df_pathGO_SOFT_red <- recompte_path(res = df_pathGO_SOFT, tissues_name = tissues_name, id_pval = 15)
df_pathWiki_red <- recompte_path(res = df_pathWiki, tissues_name = tissues_name, id_pval = 13)


### save results
write.table(x = df_tscore_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/tscore_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_tscore_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/tscore_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_HARD, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_SOFT, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_SOFT_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_WikiPath2019Human_pval_CAD_HARD_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_HARD_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_HARD_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR_SOFT_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_Reactome_pval_CAD_SOFT_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_HARD_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_HARD_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_SOFT_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_GO_pval_CAD_SOFT_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki_red, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/path_WikiPath2019Human_pval_CAD_HARD_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)

#####################################
####### divide genes per loci #######
#### tissue specific and overall ####
#####################################

# load gene location
gene_loc <- list()
for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  tmp <- fread(sprintf('train_GTEx/%s/200kb/CAD_GWAS_bin5e-2/resPrior_regEval_allchr.txt', t), data.table = F, stringsAsFactors = F, h=T)
  gene_loc[[i]] <- tmp[, 1:9]
  gene_loc[[i]]$tissue <- t
} 
gene_loc <- do.call(rbind,gene_loc)
gene_loc$new_id <- paste(gene_loc$ensembl_gene_id, gene_loc$tissue, sep= '_')

# consider only significant results
gene_table <- df_tscore_HARD[df_tscore_HARD[,10] <= 0.05,]
gene_table$new_id <- paste(gene_table$ensembl_gene_id, gene_table$tissue, sep= '_')
# filter
if(all(gene_table$new_id %in% gene_loc$new_id)){
  gene_loc_filt <- gene_loc[match(gene_table$new_id, gene_loc$new_id), ]   
}else{
  stop('problem in matching genes!')
}
gene_table <- cbind(gene_loc_filt[,c('type', 'chrom', 'TSS_start', 'start_position', 'end_position')],gene_table)

# function merge loci
merge_loci_genes <- function(gene_table, cis_size, bp_loci, tissue = 'combined'){
  
  tmp <- gene_table
  tmp_loci <- data.frame(chrom = c(), start = c(), end = c(), ngenes = c(), gene = c(), tissue = c())
  
  # divide per chr
  chr_id <- unique(tmp$chrom)
  tmp_chr <- lapply(chr_id, function(x) tmp[tmp$chrom == x,])
  
  for(j in 1:length(chr_id)){
    print(j)
    if(nrow(tmp_chr[[j]]) == 1){
      
      tmp_loci <- rbind(tmp_loci, data.frame(chrom = chr_id[j], start = tmp_chr[[j]]$TSS_start - cis_size, end = tmp_chr[[j]]$TSS_start + cis_size, 
                                             ngenes = 1, gene = tmp_chr[[j]]$external_gene_name, ensembl_gene = tmp_chr[[j]]$ensembl_gene_id,
                                             tissue = tissue))  
    }else{
      
      tmp_chr[[j]] <- tmp_chr[[j]][order(tmp_chr[[j]]$TSS_start), ]
      reg_gene <- data.frame(start = tmp_chr[[j]]$TSS_start - cis_size,  end = tmp_chr[[j]]$TSS_start + cis_size)
      merg_cond <- sapply(reg_gene$end, function(x) abs(x-reg_gene$start) < bp_loci) # the end of the second genes is close to the start of the first gene 1Mb
      
      merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
      merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
      merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
      
      merge_pos <- lapply(merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
      new_merge_pos <- list()
      all_merg <- F
      it <- 0
      
      if(length(merge_pos)>1){
        while(!all_merg){
          
          it <- it+1
          # print(it)
          
          for(l in 1:(length(merge_pos)-1)){
            
            if(!all(is.na(merge_pos[[l]]))){
              
              if(all(!merge_pos[[l]] %in% merge_pos[[l+1]])){
                new_merge_pos <- list.append(new_merge_pos, merge_pos[[l]])
              }else{
                if(!(all(merge_pos[[l]] %in% merge_pos[[l+1]]) | all(merge_pos[[l+1]] %in% merge_pos[[l]]))){
                  new_merge_pos <- list.append(new_merge_pos, unique(c(merge_pos[[l]], merge_pos[[l+1]])))
                }else{
                  if(all(merge_pos[[l+1]] %in% merge_pos[[l]])){
                    merge_pos[[l+1]] <- NA
                    new_merge_pos <- list.append(new_merge_pos, merge_pos[[l]])
                  }
                }
              }
              
            }
            
          }
          
          new_merge_pos <- list.append(new_merge_pos, merge_pos[[length(merge_pos)]])
          
          all_merg <- all(!duplicated(unlist(new_merge_pos)))
          merge_pos <- new_merge_pos
          new_merge_pos <- list() 
          
        }
        
        # remove NA
        merge_pos <- merge_pos[!sapply(merge_pos, function(x) all(is.na(x)))]
      }
      tmp_res <-  lapply(merge_pos, function(x) data.frame(chrom = chr_id[j], start = min(tmp_chr[[j]]$TSS_start[x] - cis_size), 
                                                           end = max(tmp_chr[[j]]$TSS_start[x] + cis_size), 
                                                           ngenes = length(x),
                                                           gene = paste0(unique(tmp_chr[[j]]$external_gene_name[x]), collapse = ','), 
                                                           ensembl_gene = paste0(unique(tmp_chr[[j]]$ensembl_gene_id[x]), collapse = ','), 
                                                           tissue = tissue))
      tmp_loci <-  rbind(tmp_loci,do.call(rbind, tmp_res))
      
    }
    
  }
  tmp_loci$start[tmp_loci$start < 0] <- 0
  tmp_loci$loci_id <- paste0(tmp_loci$chrom,':',round(tmp_loci$start/1000000, digits = 1), '-', round(tmp_loci$end/1000000, digits = 1), 'Mb')
  tmp_loci$loci_complete <- paste0(tmp_loci$chrom,':',tmp_loci$start,'-',tmp_loci$end)
  
  return(tmp_loci)
}
# tissue specific
tissue_spec_loci <- list()
for(i in 1:length(tissues_name)){
  
  print(tissues_name[i])
  tissue_spec_loci[[i]] <- merge_loci_genes(gene_table = gene_table[gene_table$tissue %in% tissues_name[i], ],
                                            cis_size = cis_size, bp_loci = bp_loci, tissue = tissues_name[i])
}

tissue_spec_loci <- do.call(rbind, tissue_spec_loci)
# annotate significant table for loci id
gene_table$loci_tissue_specific <- NA
for(t in tissues_name){
  id <- which(gene_table$tissue %in% t)
  for(i in id){
    gene_table$loci_tissue_specific[i] <- tissue_spec_loci$loci_id[tissue_spec_loci$tissue == t & 
                                                                     grepl(gene_table$ensembl_gene_id[i], tissue_spec_loci$ensembl_gene)]  
  }
}

# all tissues
all_loci <- merge_loci_genes(gene_table = gene_table, cis_size = cis_size, bp_loci = bp_loci)
# annotate significant table for loci id
gene_table$loci <- NA
gene_table$loci_complete <- NA
for(i in 1:nrow(gene_table)){
  gene_table$loci[i] <- all_loci$loci_id[grepl(gene_table$ensembl_gene_id[i], all_loci$ensembl_gene)]  
  gene_table$loci_complete[i] <- all_loci$loci_complete[grepl(gene_table$ensembl_gene_id[i], all_loci$ensembl_gene)] 
}
gene_table <- gene_table[, !colnames(gene_table) %in% 'new_id']

# save
write.table(x = gene_table, file = 'predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/GeneTscores_CADHARD_annotated.txt', col.names=T, row.names=F, sep = '\t', quote = F)


