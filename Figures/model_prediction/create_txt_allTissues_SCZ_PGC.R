# combine all tissues results (SCZ)
library(qvalue)
library(data.table)
library(rlist)

setwd('/psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ')
tissues_name <- read.table('Tissues_PGC_red2', h=F, stringsAsFactors = F)$V1
cis_size <- 200000
bp_loci <- 1000000

df_tscore <- df_pathR <- df_pathGO <- df_pathWiki <- list()

for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr.RData', t)))
  df_tscore[[i]] <- tmp$tscore[[1]]
  df_tscore[[i]]$tissue <- t
  
  # add genes in the pathway and if there is an improvment in significance
  df_pathR[[i]] <- tmp$pathScore_reactome[[1]]
  df_pathR[[i]]$genes_path <- NA
  df_pathR[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathR[[i]])){
    df_pathR[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_reactome[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathR[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_reactome[[1]][[j]]$tscore[,8] > df_pathR[[i]][j,13])
  }
  df_pathR[[i]]$tissue <- t
  
  df_pathGO[[i]] <- tmp$pathScore_GO[[1]]
  df_pathGO[[i]]$genes_path <- NA
  df_pathGO[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathGO[[i]])){
    df_pathGO[[i]]$genes_path[j] <- paste0(tmp$info_pathScore_GO[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathGO[[i]]$improvement_sign[j] <- all(tmp$info_pathScore_GO[[1]][[j]]$tscore[,8] > df_pathGO[[i]][j,15])
  }
  df_pathGO[[i]]$tissue <- t
  
  
  tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr_customPath_WikiPath2019Human.RData', t)))
  df_pathWiki[[i]] <- tmp$pathScore[[1]]
  df_pathWiki[[i]]$genes_path <- NA
  df_pathWiki[[i]]$improvement_sign  <- NA
  for(j in 1:nrow(df_pathWiki[[i]])){
    df_pathWiki[[i]]$genes_path[j] <- paste0(tmp$info_pathScore[[1]][[j]]$tscore$external_gene_name, collapse = ',')
    df_pathWiki[[i]]$improvement_sign[j] <- all(tmp$info_pathScore[[1]][[j]]$tscore[,8] > df_pathWiki[[i]][j,13])
  }
  df_pathWiki[[i]]$tissue <- t
}

df_tscore <- do.call(rbind, df_tscore)
df_pathR <- do.call(rbind, df_pathR)
df_pathGO <- do.call(rbind, df_pathGO)
df_pathWiki <- do.call(rbind, df_pathWiki)


# add overall corrected pval
df_tscore$Dx_pval_BHcorr_overall <- p.adjust(df_tscore[, 8], method = 'BH')
df_pathR$Dx_pval_BHcorr_overall <- p.adjust(df_pathR[, 13], method = 'BH')
df_pathGO$Dx_pval_BHcorr_overall <- p.adjust(df_pathGO[, 15], method = 'BH')
df_pathWiki$Dx_pval_BHcorr_overall <- p.adjust(df_pathWiki[, 13], method = 'BH')


# create a function to remove pathway with 1 gene and recompute pvalues
recompte_path <- function(tissues_name, res, id_pval){
  tmp <- lapply(tissues_name, function(x) res[res$tissue == x & res$ngenes_tscore>1,])
  for(i in 1:length(tmp)){
    tmp[[i]][, id_pval+1] <- qvalue(tmp[[i]][, id_pval])$qvalue
    tmp[[i]][, id_pval+2] <- p.adjust(tmp[[i]][, id_pval], method = 'BH')
  }
  tmp <- do.call(rbind, tmp)
  tmp[, id_pval+9] <- p.adjust(tmp[, id_pval], method = 'BH')
  colnames(tmp)[ncol(tmp)] <- 'Dx_pval_BHcorr_overall'
  return(tmp)
}

df_pathR_red <- recompte_path(res = df_pathR, tissues_name = tissues_name, id_pval = 13)
df_pathGO_red <- recompte_path(res = df_pathGO, tissues_name = tissues_name, id_pval = 15)
df_pathWiki_red <- recompte_path(res = df_pathWiki, tissues_name = tissues_name, id_pval = 13)


### save results
write.table(x = df_tscore, file = 'OUTPUT_all/tscore_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathR, file = 'OUTPUT_all/path_Reactome_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO, file = 'OUTPUT_all/path_GO_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki, file = 'OUTPUT_all/customPath_WikiPath2019Human_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)

write.table(x = df_pathR_red, file = 'OUTPUT_all/path_Reactome_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathGO_red, file = 'OUTPUT_all/path_GO_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathWiki_red, file = 'OUTPUT_all/customPath_WikiPath2019Human_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)

# CMC pathways
t <- tissues_name[1]
tmp <- get(load(sprintf('%s/pval_Dx_pheno_covCorr_customPath_CMC_GeneSets.RData', t)))
df_pathCMC <- tmp$pathScore[[1]]
df_pathCMC$genes_path <- NA
df_pathCMC$improvement_sign  <- NA
for(j in 1:nrow(df_pathCMC)){
  df_pathCMC$genes_path[j] <- paste0(tmp$info_pathScore[[1]][[j]]$tscore$external_gene_name, collapse = ',')
  df_pathCMC$improvement_sign[j] <- all(tmp$info_pathScore[[1]][[j]]$tscore[,8] > df_pathCMC[j,13])
}
df_pathCMC$tissue <- t

df_pathCMC_red <- recompte_path(res = df_pathCMC, tissues_name = t, id_pval = 13)
write.table(x = df_pathCMC, file = 'DLPC_CMC/customPath_CMC_GeneSets_pval_SCZ_covCorr.txt', col.names=T, row.names=F, sep = '\t', quote = F)
write.table(x = df_pathCMC_red, file = 'DLPC_CMC/customPath_CMC_GeneSets_pval_SCZ_covCorr_filt.txt', col.names=T, row.names=F, sep = '\t', quote = F)

#####################################
####### divide genes per loci #######
#### tissue specific and overall ####
#####################################

# load gene location
gene_loc <- list()
for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  tmp <- fread(sprintf('%s/resPrior_regEval_allchr.txt', t), data.table = F, stringsAsFactors = F, h=T)
  gene_loc[[i]] <- tmp[, 1:9]
  gene_loc[[i]]$tissue <- t
} 
gene_loc <- do.call(rbind,gene_loc)
gene_loc$new_id <- paste(gene_loc$ensembl_gene_id, gene_loc$tissue, sep= '_')

# consider only significant results
gene_table <- df_tscore[df_tscore[,10] <= 0.05,]
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
  tmp_loci$loci_complete <- paste0(tmp_loci$chrom,':',tmp_loci$start, '-', tmp_loci$end)
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
write.table(x = gene_table, file = 'OUTPUT_all/GeneTscores_SCZ_annotated.txt', col.names=T, row.names=F, sep = '\t', quote = F)


