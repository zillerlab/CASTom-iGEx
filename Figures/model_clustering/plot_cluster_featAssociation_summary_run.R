options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rlist))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="overal results phenotype association")
parser$add_argument("--gene_feat_file", type = "character", help = "")
parser$add_argument("--gene_info_file", type = "character", help = "")
parser$add_argument("--pathR_feat_file", type = "character", help = "")
parser$add_argument("--pathR_info_file", type = "character", help = "")
parser$add_argument("--pathGO_feat_file", type = "character", help = "")
parser$add_argument("--pathGO_info_file", type = "character", help = "")
parser$add_argument("--pheno_name", type = "character", help = "")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--type_input", type = "character", help = "")
parser$add_argument("--type_cluster", type = "character", help = "")
parser$add_argument("--cis_size", default = 200000, type = "character", help = "")
parser$add_argument("--color_tissues_file", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
gene_feat_file <- args$gene_feat_file
gene_info_file <- args$gene_info_file
pathR_feat_file <- args$pathR_feat_file
pathR_info_file <- args$pathR_info_file
pathGO_feat_file <- args$pathGO_feat_file
pathGO_info_file <- args$pathGO_info_file
pheno_name <- args$pheno_name
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
color_tissues_file <- args$color_tissues_file
cis_size <- args$cis_size
outFold <- args$outFold

########################################################################################################################
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# pheno_name <- 'CAD_HARD'
# color_tissues_file <- '/psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt'
# gene_info_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_tscoreOriginal_tscoreClusterCases_infoGenes.txt'
# gene_feat_file  <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_tscoreOriginal_tscoreClusterCases_featAssociation.txt'
# pathR_info_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_ReactomeOriginal_tscoreClusterCases_infoGenes.txt'
# pathR_feat_file  <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_ReactomeOriginal_tscoreClusterCases_featAssociation.txt'
# pathGO_info_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_GOOriginal_tscoreClusterCases_infoGenes.txt'
# pathGO_feat_file  <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_GOOriginal_tscoreClusterCases_featAssociation.txt'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_'
# cis_size <- 200000
########################################################################################################################

# load results
gene_loc = read.delim(gene_info_file, h=T, stringsAsFactors = F, sep = '\t')
gene_feat = read.delim(gene_feat_file, h=T, stringsAsFactors = F, sep = '\t')
gene_loc$new_id <- paste0(gene_loc$external_gene_name, '_tissue_', gene_loc$tissue)
gene_feat$new_id <- paste0(gene_feat$feat, '_tissue_', gene_feat$tissue)

pathR_info = read.delim(pathR_info_file, h=T, stringsAsFactors = F, sep = '\t')
pathR_feat = read.delim(pathR_feat_file, h=T, stringsAsFactors = F, sep = '\t')
pathR_info$new_id <- paste0(pathR_info$path, '_tissue_', pathR_info$tissue)
pathR_feat$new_id <- paste0(pathR_feat$feat, '_tissue_', pathR_feat$tissue)

pathGO_info = read.delim(pathGO_info_file, h=T, stringsAsFactors = F, sep = '\t')
pathGO_feat = read.delim(pathGO_feat_file, h=T, stringsAsFactors = F, sep = '\t')
pathGO_info$new_id <- paste0(pathGO_info$path, '_tissue_', pathGO_info$tissue)
pathGO_feat$new_id <- paste0(pathGO_feat$feat, '_tissue_', pathGO_feat$tissue)

tissues <- unique(gene_feat$tissue)
color_tissues <- read.table(color_tissues_file, h=T, stringsAsFactors = F)

# for gene and each tissue, put count the number of associations and divide per loci
### gene ###
tissue_spec_loci <- list()
for(i in 1:length(tissues)){
  
  print(tissues[i])
  tmp_info <- gene_loc[gene_loc$tissue == tissues[i], ]
  tmp_feat <- gene_feat[gene_feat$tissue == tissues[i], ]
  # consider only significant differences
  tmp_feat <- tmp_feat[tmp_feat$pval_corr <= 0.05,]
  tmp_info <- tmp_info[tmp_info$new_id %in% tmp_feat$new_id,]
  
  # combine loci
  # MHC_loci <- tmp_info[tmp_info$chrom == 'chr6' & (tmp_info$TSS_start - cis_size > 26000000 & tmp_info$TSS_start + cis_size < 34000000), ]
  MHC_loci <- data.frame() # MHC will be automatically detected
  if(nrow(MHC_loci)>0){
    tmp_info_noMHC <- tmp_info[!tmp_info$new_id %in% MHC_loci$new_id, ] 
    tmp_loci <- data.frame(chrom = 'chr6', start = min(MHC_loci$TSS_start - cis_size), end = max(MHC_loci$TSS_start + cis_size), 
                           ngenes = nrow(MHC_loci), gene = paste0(MHC_loci$external_gene_name, collapse = ','), mean_Zstat = mean(MHC_loci$Zstat), 
                           sd_Zstat = sd(MHC_loci$Zstat), highest_Zstat = MHC_loci$Zstat[which.max(abs(MHC_loci$Zstat))],tissue = tissues[i]) 
  }else{
    tmp_info_noMHC <- tmp_info  
    tmp_loci <- data.frame(chrom = c(), start = c(), end = c(), ngenes = c(), gene = c(),  mean_Zstat = c(), sd_Zstat = c(), highest_Zstat = c(), tissue = c())
  }
  
  # divide per chr
  chr_id <- unique(tmp_info_noMHC$chrom)
  tmp_info_noMHC <- lapply(chr_id, function(x) tmp_info_noMHC[tmp_info_noMHC$chrom == x,])
  
  for(j in 1:length(chr_id)){
    # print(j)
    if(nrow(tmp_info_noMHC[[j]]) == 1){
      
      tmp_loci <- rbind(tmp_loci, data.frame(chrom = chr_id[j], start = tmp_info_noMHC[[j]]$TSS_start - cis_size, end = tmp_info_noMHC[[j]]$TSS_start + cis_size, 
                                             ngenes = 1, gene = tmp_info_noMHC[[j]]$external_gene_name, mean_Zstat = tmp_info_noMHC[[j]]$Zstat, 
                                             sd_Zstat = NA, highest_Zstat = tmp_info_noMHC[[j]]$Zstat, tissue = tissues[i]))  
    }else{
      
      tmp_info_noMHC[[j]] <- tmp_info_noMHC[[j]][order(tmp_info_noMHC[[j]]$TSS_start), ]
      reg_gene <- data.frame( start = tmp_info_noMHC[[j]]$TSS_start - cis_size,  end = tmp_info_noMHC[[j]]$TSS_start + cis_size)
      merg_cond <- sapply(reg_gene$end, function(x) abs(x-reg_gene$start) < 1000000) # the end of the second genes is close to the start of the first gene 1Mb
      
      merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
      merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
      merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
      
      merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
      new_merge_pos <- list()
      all_merg <- F
      it <- 0
      
      if(length(merge_pos)>1){
        while(!all_merg){
          
          it <- it+1
          print(it)
          
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
      tmp <-  lapply(merge_pos, function(x) data.frame(chrom = chr_id[j], start = min(tmp_info_noMHC[[j]]$TSS_start[x] - cis_size), end = max(tmp_info_noMHC[[j]]$TSS_start[x] + cis_size), 
                                                       ngenes = length(x), gene = paste0(tmp_info_noMHC[[j]]$external_gene_name[x], collapse = ','), 
                                                       mean_Zstat = mean(tmp_info_noMHC[[j]]$Zstat[x]), sd_Zstat = sd(tmp_info_noMHC[[j]]$Zstat[x]), 
                                                       highest_Zstat = tmp_info_noMHC[[j]]$Zstat[x][which.max(abs(tmp_info_noMHC[[j]]$Zstat[x]))], tissue = tissues[i]))
      tmp_loci <-  rbind(tmp_loci,do.call(rbind, tmp))
      
    }
    
  }
  
  tmp_loci$loci_id <- paste0(tmp_loci$chrom,':',round(tmp_loci$start/1000000, digits = 1), '-', round(tmp_loci$end/1000000, digits = 1), 'Mb')
  tissue_spec_loci[[i]] <- tmp_loci
  
  tissue_spec_loci[[i]]$comp_sign <- NA
  # for each loci, find groups that are significantly different
  for(l in 1:nrow(tissue_spec_loci[[i]])){
    genes <- strsplit(tissue_spec_loci[[i]]$gene[l], split = ',')[[1]]
    tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% genes & tmp_feat$pval_corr <= 0.05]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
    tissue_spec_loci[[i]]$comp_sign[l] <-  paste0(tmp_gr, collapse = ',') 
  }
  
  
}

tissue_spec_loci <- do.call(rbind, tissue_spec_loci)
# save
write.table(x = tissue_spec_loci, file = sprintf('%s%s_%s_cluster%s_summary_geneLoci_tissueSpec.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

#### all tisues ####
tmp_info <- gene_loc
tmp_feat <- gene_feat
# consider only significant differences
tmp_feat <- tmp_feat[tmp_feat$pval_corr <= 0.05,]
tmp_info <- tmp_info[tmp_info$new_id %in% tmp_feat$new_id,]

# combine loci
# MHC_loci <- tmp_info[tmp_info$chrom == 'chr6' & (tmp_info$TSS_start - cis_size > 26000000 & tmp_info$TSS_start + cis_size < 34000000), ]
MHC_loci <- data.frame()
if(nrow(MHC_loci)>0){
  tmp_info_noMHC <- tmp_info[!tmp_info$new_id %in% MHC_loci$new_id, ] 
  tmp_loci <- data.frame(chrom = 'chr6', start = min(MHC_loci$TSS_start - cis_size), end = max(MHC_loci$TSS_start + cis_size), 
                         ngenes_withrep = nrow(MHC_loci), ngenes_unique = length(unique(MHC_loci$external_gene_name)), gene = paste0(unique(MHC_loci$external_gene_name), collapse = ','), mean_Zstat = mean(MHC_loci$Zstat), 
                         sd_Zstat = sd(MHC_loci$Zstat), highest_Zstat = MHC_loci$Zstat[which.max(abs(MHC_loci$Zstat))]) 
}else{
  tmp_info_noMHC <- tmp_info  
  tmp_loci <- data.frame(chrom = c(), start = c(), end = c(),  ngenes_withrep = c(), ngenes_unique = c(), gene = c(),  mean_Zstat = c(), sd_Zstat = c(), highest_Zstat = c())
}

# divide per chr
chr_id <- unique(tmp_info_noMHC$chrom)
tmp_info_noMHC <- lapply(chr_id, function(x) tmp_info_noMHC[tmp_info_noMHC$chrom == x,])

for(j in 1:length(chr_id)){
  print(chr_id[j])
  if(nrow(tmp_info_noMHC[[j]]) == 1){
    
    tmp_loci <- rbind(tmp_loci, data.frame(chrom = chr_id[j], start = tmp_info_noMHC[[j]]$TSS_start - cis_size, end = tmp_info_noMHC[[j]]$TSS_start + cis_size, 
                                           ngenes_withrep = 1, ngenes_unique = 1, gene = tmp_info_noMHC[[j]]$external_gene_name, mean_Zstat = tmp_info_noMHC[[j]]$Zstat, 
                                           sd_Zstat = NA, highest_Zstat = tmp_info_noMHC[[j]]$Zstat))  
  }else{
    
    tmp_info_noMHC[[j]] <- tmp_info_noMHC[[j]][order(tmp_info_noMHC[[j]]$TSS_start), ]
    reg_gene <- data.frame( start = tmp_info_noMHC[[j]]$TSS_start - cis_size,  end = tmp_info_noMHC[[j]]$TSS_start + cis_size)
    merg_cond <- sapply(reg_gene$end, function(x) abs(x-reg_gene$start) < 1000000) # the end of the second genes is close to the start of the first gene 1Mb
    
    merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
    merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
    merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
    
    merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
    new_merge_pos <- list()
    all_merg <- F
    it <- 0
    
    while(!all_merg){
      
      it <- it+1
      print(it)
      
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
    
    tmp <-  lapply(merge_pos, function(x) data.frame(chrom = chr_id[j], start = min(tmp_info_noMHC[[j]]$TSS_start[x] - cis_size), end = max(tmp_info_noMHC[[j]]$TSS_start[x] + cis_size), 
                                                     ngenes_withrep = length(x), ngenes_unique = length(unique(tmp_info_noMHC[[j]]$external_gene_name[x])), gene = paste0(unique(tmp_info_noMHC[[j]]$external_gene_name[x]), collapse = ','), 
                                                     mean_Zstat = mean(tmp_info_noMHC[[j]]$Zstat[x]), sd_Zstat = sd(tmp_info_noMHC[[j]]$Zstat[x]), 
                                                     highest_Zstat = tmp_info_noMHC[[j]]$Zstat[x][which.max(abs(tmp_info_noMHC[[j]]$Zstat[x]))]))
    tmp_loci <-  rbind(tmp_loci,do.call(rbind, tmp))
    
  }
  
}

tmp_loci$loci_id <- paste0(tmp_loci$chrom,':',round(tmp_loci$start/1000000, digits = 1), '-', round(tmp_loci$end/1000000, digits = 1), 'Mb')
alltissues_loci <- tmp_loci 
alltissues_loci$comp_sign <- NA
# for each loci, find groups that are significantly different
for(i in 1:nrow(alltissues_loci)){
  genes <- strsplit(alltissues_loci$gene[i], split = ',')[[1]]
  tmp_gr <- sapply(unique(gene_feat$comp[gene_feat$feat %in% genes & gene_feat$pval_corr <= 0.05]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
  alltissues_loci$comp_sign[i] <-  paste0(tmp_gr, collapse = ',') 
}
# save
write.table(x = alltissues_loci, file = sprintf('%s%s_%s_cluster%s_summary_geneLoci_allTissues.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)


#############################################################################################
### plot distributions ###
gr_tot <- unique(gene_feat$comp)
gr_tot <- sapply(gr_tot, function(x) strsplit(x, split = '_vs_all')[[1]][1])
df <- data.frame(ngenes = c(), tissue = c(), comp = c(), loci = c())

for(i in 1:length(gr_tot)){
  tmp <- alltissues_loci[grepl(gr_tot[i], alltissues_loci$comp),]
  df <- rbind(df, data.frame(ngenes = tmp$ngenes_unique, tissue = rep('All', nrow(tmp)), comp = rep(gr_tot[i], nrow(tmp)), loci = tmp$loci_id))
}

for(j in 1:length(tissues)){
  for(i in 1:length(gr_tot)){
    tmp <- tissue_spec_loci[grepl(gr_tot[i], tissue_spec_loci$comp) & tissue_spec_loci$tissue %in% tissues[j],]
    df <- rbind(df, data.frame(ngenes = tmp$ngenes, tissue = rep(tissues[j], nrow(tmp)), comp = rep(gr_tot[i], nrow(tmp)), loci = tmp$loci_id))
  }
}

df$comp <- factor(df$comp, levels = unname(gr_tot))
df$tissue <- factor(df$tissue, levels = c('All', tissues))
newcolours <- c('grey', color_tissues$color[match(tissues, color_tissues$tissues)])
# put loci in id based on overall comparison
df$loci_name <- df$loci
id_notall <- which(df$tissue !='All')
for(i in id_notall){
  tmp <- tissue_spec_loci[tissue_spec_loci$loci_id == df$loci[i] & grepl(df$comp[i], tissue_spec_loci$comp_sign) & tissue_spec_loci$tissue == df$tissue[i] ,]
  tmp_all <- alltissues_loci[alltissues_loci$chrom == tmp$chrom,]
  tmp_all <- tmp_all[which.min(abs(tmp_all$start - tmp$start)),]
  df$loci_name[i] <- tmp_all$loci_id
}
df$loci_name <- factor(df$loci_name, levels =unique( alltissues_loci$loci_id))

# same plot but n. of loci and not genes
df_loci <- data.frame(tissue = unlist(lapply(c('All', tissues), function(x) rep(x, length(gr_tot)))), comp = rep(gr_tot, length(tissues)+1), 
                      nloci = as.vector(table(df$comp, df$tissue)))
df_loci$comp <- factor(df_loci$comp, levels = unname(gr_tot))
df_loci$tissue <- factor(df_loci$tissue, levels = c('All', tissues))
coul <- colorRampPalette(brewer.pal(11, "Spectral"))(length(levels(df$loci_name)))
color_gr <- pal_d3("category10")(length(gr_tot))

pl1 <- ggplot(data = subset(df, tissue != 'All'), aes(x = tissue, y = ngenes, fill = loci_name))+
  geom_bar(alpha = 0.9, width = 0.5, stat = 'identity')+
  facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. genes cluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'right', legend.key.size = unit(0.2, "cm"), 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours[-1]))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values = coul)+
  # scale_fill_d3()+
  coord_flip() 

# pl1 <- ggplot_gtable(ggplot_build(pl1))
# stripr <- which(grepl('strip-t', pl1$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl1$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

pl2 <- ggplot(data = subset(df_loci, tissue != 'All'), aes(x = tissue, y = nloci))+
  geom_bar(alpha = 0.9, width = 0.5, stat = 'identity', fill = 'darkgrey')+
  facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. loci cluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'right', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours[-1]))+
  # guides(fill=guide_legend(ncol=2))+
  # scale_fill_d3()+
  coord_flip() 

# pl2 <- ggplot_gtable(ggplot_build(pl2))
# stripr <- which(grepl('strip-t', pl2$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=1, nrow=2, align = 'v', common.legend = T, legend = 'right')
# tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=1, nrow=2, align = 'v')
ggsave(filename =  sprintf('%scluster_ngenes_nloci_pertissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 10, height = 7, dpi = 500)
ggsave(filename = sprintf('%scluster_ngenes_nloci_pertissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 10, height = 7, dpi = 500, compress = F)


pl1 <- ggplot(data = subset(df, tissue == 'All'), aes(x = comp, y = ngenes, fill = loci_name))+
  geom_bar(alpha = 0.9, width = 0.5, stat = 'identity')+
  ylab('n. genes cluster relevant')+ 
  theme_bw()+ 
  # ggtitle('All tissues combined')+
  theme(legend.position = 'none', legend.key.size = unit(0.2, "cm"), plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(colour = color_gr))+
  scale_fill_manual(values = coul)+
  guides(fill=guide_legend(ncol=2))
  # scale_fill_d3()+

pl2 <- ggplot(data = subset(df_loci, tissue == 'All'), aes(x = comp, y = nloci))+
  geom_bar(alpha = 0.9, width = 0.5, stat = 'identity', fill = 'darkgrey')+
  # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. loci cluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(colour = color_gr))+
  # guides(fill=guide_legend(ncol=2))+
  scale_fill_d3()

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h')
ggsave(filename =  sprintf('%scluster_ngenes_nloci_alltissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 4, height = 5, dpi = 500)
ggsave(filename = sprintf('%scluster_ngenes_nloci_alltissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 4, height = 5, dpi = 500, compress = F)



# # tot_pl <- ggarrange(plotlist = list(pl, pl_side), ncol=2, nrow=1, widths=c(1, 0.3), align = 'h', common.legend = TRUE)
# tot_pl <- pl
# ggsave(filename =  sprintf('%scluster_nloci_pertissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 10, height = 4, dpi = 500)
# ggsave(filename = sprintf('%scluster_nloci_pertissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 10, height = 4, dpi = 500, compress = F)

# # heatmap version
# save_pheatmap_png <- function(x, filename, width=10, height=10, res = 200) {
#   png(filename, width = width, height = height, res = res, units = 'in')
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
#   pdf(filename, width = width, height = height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# 
# mat_loci <- t(table(df$comp, df$tissue))
# coul <- colorRampPalette(brewer.pal(9, "Blues"))(100)
# hm_pl <- pheatmap(mat=mat_loci, show_colnames = T, color=coul, show_rownames = T, 
#                   cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', display_numbers = T, 
#                   drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
#                   main = sprintf(''), treeheight_row = 0, treeheight_col = 0, cellwidth = 25)
# save_pheatmap_png(hm_pl,sprintf('%scluster_nloci_pertissues_heatmap_%s_%s.png', outFold, type_data, type_input), height = 5, width =7, res = 200)
# save_pheatmap_pdf(hm_pl, sprintf('%scluster_nloci_pertissues_heatmap_%s_%s.pdf', outFold, type_data, type_input), height = 5, width =7)

#################################################################################################################################################
## find pathways differences
tissue_spec_pathR <- list()

for(i in 1:length(tissues)){
  
  print(tissues[i])
  tmp_info <- pathR_info[pathR_info$tissue == tissues[i], ]
  tmp_feat <- pathR_feat[pathR_feat$tissue == tissues[i], ]
  # consider only significant differences
  tmp_feat <- tmp_feat[tmp_feat$pval_corr <= 0.05,]
  tmp_info <- tmp_info[tmp_info$new_id %in% tmp_feat$new_id,]
  
  # combine in group, any pathways sharing a gene
  tmp_path <- data.frame(npath = c(), path = c(),  mean_Zstat = c(), sd_Zstat = c(), highest_Zstat = c(), mean_ngenes = c(), sd_ngenes = c(), tissue = c())
  gene_names <- lapply(tmp_info$genes_id, function(x) strsplit(x, split = ',')[[1]])
  merg_cond <- sapply(gene_names, function(x) sapply(gene_names, function(y) length(intersect(x,y))/length(union(x,y)) >= 0.3))
  
  merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
  merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
  merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
  merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
  new_merge_pos <- list()
  all_merg <- F
  it <- 0
  
  if(length(merge_pos)>1){
    while(!all_merg){
      
      it <- it+1
      print(it)
      
      for(l in 1:(length(merge_pos))){
        
        if(!length(which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))) >1){
          id <- sapply(merge_pos, function(x) any(merge_pos[[l]] %in% x))
          new_merge_pos <- list.append(new_merge_pos, sort(unique(unlist(merge_pos[which(id)]))))
        }else{
          id_b <- which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))
          id_b <- id_b[which.max(sapply(merge_pos[id_b], length))]
          new_merge_pos <- list.append(new_merge_pos, sort(merge_pos[[id_b]]))
        }
        
      }
      
      all_merg <- all(!duplicated(unlist(new_merge_pos)))
      merge_pos <- new_merge_pos
      merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
      merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
      merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
      print(length(merge_pos))
      new_merge_pos <- list() 
      
    }
  }
 
  tmp <-  lapply(merge_pos, function(x) data.frame(npath = length(x), path = paste0(tmp_info$path[x], collapse = '-and-'), 
                                                   mean_Zstat = mean(tmp_info$Zstat[x]), sd_Zstat = sd(tmp_info$Zstat[x]), 
                                                   highest_Zstat = tmp_info$Zstat[x][which.max(abs(tmp_info$Zstat[x]))], 
                                                   mean_ngenes = mean(tmp_info$ngenes_tscore[x]), sd_ngenes = sd(tmp_info$ngenes_tscore[x]), 
                                                   tissue = tissues[i]))
  
  tmp_path <-  rbind(tmp_path,do.call(rbind, tmp))
  tissue_spec_pathR[[i]] <- tmp_path
  
  tissue_spec_pathR[[i]]$comp_sign <- NA
  # for each loci, find groups that are significantly different
  for(l in 1:nrow(tissue_spec_pathR[[i]])){
    paths <- strsplit(tissue_spec_pathR[[i]]$path[l], split = '-and-')[[1]]
    tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% paths & tmp_feat$pval_corr <= 0.05]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
    tissue_spec_pathR[[i]]$comp_sign[l] <-  paste0(tmp_gr, collapse = ',') 
  }
}

tissue_spec_pathR <- do.call(rbind, tissue_spec_pathR)
# save
write.table(x = tissue_spec_pathR, file = sprintf('%s%s_%s_cluster%s_summary_path_Reactome_tissueSpec.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)

# ###### all tissues ######### # not reasonable to pull tissues together, genes are regulated differently
# tmp_info <- pathR_info
# tmp_feat <- pathR_feat
# 
# # consider only significant differences
# tmp_feat <- tmp_feat[tmp_feat$pval_corr <= 0.05,]
# tmp_info <- tmp_info[tmp_info$new_id %in% tmp_feat$new_id,]
#   
# # combine in group, any pathways sharing a gene
# gene_names <- lapply(tmp_info$genes_id, function(x) strsplit(x, split = ',')[[1]])
# merg_cond <- sapply(gene_names, function(x) sapply(gene_names, function(y) length(intersect(x,y))/length(union(x,y)) >= 0.3))
# 
# merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
# merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
# merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
# merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
# new_merge_pos <- list()
# all_merg <- F
# it <- 0
# 
# if(length(merge_pos)>1){
#   while(!all_merg){
#     
#     it <- it+1
#     print(it)
#     
#     for(l in 1:(length(merge_pos))){
#       
#       if(!length(which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))) >1){
#         id <- sapply(merge_pos, function(x) any(merge_pos[[l]] %in% x))
#         new_merge_pos <- list.append(new_merge_pos, sort(unique(unlist(merge_pos[which(id)]))))
#       }else{
#         id_b <- which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))
#         id_b <- id_b[which.max(sapply(merge_pos[id_b], length))]
#         new_merge_pos <- list.append(new_merge_pos, sort(merge_pos[[id_b]]))
#       }
#       
#     }
#     
#     all_merg <- all(!duplicated(unlist(new_merge_pos)))
#     merge_pos <- new_merge_pos
#     merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
#     merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
#     merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
#     print(length(merge_pos))
#     new_merge_pos <- list() 
#     
#   }
# }
# 
# tmp <-  lapply(merge_pos, function(x) data.frame(npath_withrep = length(x), npath_unique = length(unique(tmp_info$path[x])), 
#                                                  ntissues = length(unique(tmp_info$tissue[x])),
#                                                  path = paste0(unique(tmp_info$path[x]), collapse = '-and-'), tissues_path = paste0(unique(tmp_info$tissue[x]), collapse = '-and-'), 
#                                                  mean_Zstat = mean(tmp_info$Zstat[x]), sd_Zstat = sd(tmp_info$Zstat[x]), 
#                                                  highest_Zstat = tmp_info$Zstat[x][which.max(abs(tmp_info$Zstat[x]))], 
#                                                  mean_ngenes = mean(tmp_info$ngenes_tscore[x]), sd_ngenes = sd(tmp_info$ngenes_tscore[x]), 
#                                                  tissue = 'AllTissues'))
# 
# tmp_path <-  do.call(rbind, tmp)
# alltissues_pathR <- tmp_path
# 
# alltissues_pathR$comp_sign <- NA
# # for each loci, find groups that are significantly different
# for(l in 1:nrow(alltissues_pathR)){
#   paths <- strsplit(alltissues_pathR$path[l], split = '-and-')[[1]]
#   tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% paths & tmp_feat$pval_corr <= 0.05]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
#   alltissues_pathR$comp_sign[l] <-  paste0(tmp_gr, collapse = ',') 
# }
# 
# # save
# write.table(x = alltissues_pathR, file = sprintf('%s%s_%s_cluster%s_summary_path_Reactome_allTissues.txt',outFold, type_data, type_input, type_cluster), 
#             col.names = T, row.names = F, sep = '\t', quote = F)


### plot distributions ###
# some pathway may be repeated, put together
gr_tot <- unique(pathR_feat$comp)
gr_tot <- sapply(gr_tot, function(x) strsplit(x, split = '_vs_all')[[1]][1])

df <- data.frame(npath = c(), tissue = c(), comp = c(), path_group = c())
for(j in 1:length(tissues)){
  for(i in 1:length(gr_tot)){
    tmp <- tissue_spec_pathR[grepl(gr_tot[i], tissue_spec_pathR$comp) & tissue_spec_pathR$tissue %in% tissues[j],]
    df <- rbind(df, data.frame(npath = tmp$npath, tissue = rep(tissues[j], nrow(tmp)), comp = rep(gr_tot[i], nrow(tmp)), path_group = tmp$path))
  }
}

# too complicated, do not show legend but class
# df$path_group_name <- df$path_group
# id_m <- which(df$npath > 1)
# df$path_group_name <- factor(df$path_group_name, levels = unique(df$path_group_name))

# same plot but n. of loci and not genes
df_pathgroup <- data.frame(tissue = unlist(lapply(tissues, function(x) rep(x, length(gr_tot)))), comp = rep(gr_tot, length(tissues)), 
                      ngroup = as.vector(sapply( tissues, function(x) sapply(gr_tot, function(y) sum(df$comp == y & df$tissue == x)))))

df_pathgroup$comp <- factor(df_pathgroup$comp, levels = unname(gr_tot))
df_pathgroup$tissue <- factor(df_pathgroup$tissue, levels =  tissues)

df_tot <- data.frame(tissue = unlist(lapply(tissues, function(x) rep(x, length(gr_tot)))), comp = rep(gr_tot, length(tissues)),
                     npath = as.vector(sapply(tissues, function(x) sapply(gr_tot, function(y) sum(df$npath[df$tissue == x & df$comp == y])))))
df_tot$comp <- factor(df_tot$comp, levels = unname(gr_tot))
df_tot$tissue <- factor(df_tot$tissue, levels =  tissues)
newcolours <-  color_tissues$color[match(tissues, color_tissues$tissues)]
color_gr <- pal_d3("category10")(length(gr_tot))

pl1 <- ggplot(data = subset(df_tot, tissue != 'All'), aes(x = tissue, y = npath, fill = comp))+
  geom_bar(alpha = 0.9, width = 0.9, stat = 'identity', position = position_dodge())+
  # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. pathways\ncluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'bottom', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 8), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours))+
  guides(fill=guide_legend(nrow=1))+
  scale_fill_manual(values = color_gr)+
  # scale_fill_d3()+
  coord_flip() 

# pl1 <- ggplot_gtable(ggplot_build(pl1))
# stripr <- which(grepl('strip-t', pl1$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl1$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

pl2 <- ggplot(data = subset(df_pathgroup, tissue != 'All'), aes(x = tissue, y = ngroup, fill = comp))+
  geom_bar(alpha = 0.9, width = 0.9, stat = 'identity',  position = position_dodge())+
  # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. pathway macro\ncluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours))+
  # guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values = color_gr)+
  coord_flip() 

# pl2 <- ggplot_gtable(ggplot_build(pl2))
# stripr <- which(grepl('strip-t', pl2$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h', common.legend = T, legend = 'bottom')
# tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=1, nrow=2, align = 'v')
ggsave(filename =  sprintf('%scluster_npath_npathgroup_path_Reactome_pertissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 7, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%scluster_npath_npathgroup_path_Reactome_pertissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 7, height = 4.5, dpi = 500, compress = F)


# pl1 <- ggplot(data = subset(df, tissue == 'All'), aes(x = comp, y = npath, fill = path_group_name))+
#   geom_bar(alpha = 0.9, width = 0.8, stat = 'identity')+
#   ylab('n. pathways cluster relevant')+ 
#   theme_bw()+ 
#   # ggtitle('All tissues combined')+
#   theme(legend.position = 'none', legend.key.size = unit(0.2, "cm"), plot.title = element_text(hjust = 0.5),
#         legend.text = element_text(size = 6), legend.title = element_blank(), 
#         axis.title.x = element_blank(), axis.text.x = element_text(colour = color_gr))+
#   scale_fill_manual(values = coul)+
#   guides(fill=guide_legend(ncol=2))
# # scale_fill_d3()+
# 
# pl2 <- ggplot(data = subset(df_pathgroup, tissue == 'All'), aes(x = comp, y = ngroup))+
#   geom_bar(alpha = 0.9, width = 0.5, stat = 'identity', fill = 'darkgrey')+
#   # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
#   ylab('n. pathway classes cluster relevant')+ 
#   theme_bw()+ 
#   theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
#         legend.text = element_text(size = 6), legend.title = element_blank(), 
#         axis.title.x = element_blank(), axis.text.x = element_text(colour = color_gr))+
#   # guides(fill=guide_legend(ncol=2))+
#   scale_fill_d3()
# 
# tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h')
# ggsave(filename =  sprintf('%scluster_npath_npathgroup_path_Reactome_alltissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 4, height = 5, dpi = 500)
# ggsave(filename = sprintf('%scluster_npath_npathgroup_path_Reactome_alltissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 4, height = 5, dpi = 500, compress = F)



#################################################################################################################################################
## find pathways differences GO
tissue_spec_pathGO <- list()

for(i in 1:length(tissues)){
  
  print(tissues[i])
  tmp_info <- pathGO_info[pathGO_info$tissue == tissues[i], ]
  tmp_feat <- pathGO_feat[pathGO_feat$tissue == tissues[i], ]
  # consider only significant differences
  tmp_feat <- tmp_feat[tmp_feat$pval_corr <= 0.05,]
  tmp_info <- tmp_info[tmp_info$new_id %in% tmp_feat$new_id,]
  
  # combine in group, any pathways sharing a gene
  tmp_path <- data.frame(npath = c(), path = c(),  mean_Zstat = c(), sd_Zstat = c(), highest_Zstat = c(), mean_ngenes = c(), sd_ngenes = c(), tissue = c())
  gene_names <- lapply(tmp_info$genes_id, function(x) strsplit(x, split = ',')[[1]])
  merg_cond <- sapply(gene_names, function(x) sapply(gene_names, function(y) length(intersect(x,y))/length(union(x,y)) >= 0.3))
  
  merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
  merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
  merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
  merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
  new_merge_pos <- list()
  all_merg <- F
  it <- 0
  
  if(length(merge_pos)>1){
    while(!all_merg){
      
      it <- it+1
      print(it)
      
      for(l in 1:(length(merge_pos))){
        
        if(!length(which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))) >1){
          id <- sapply(merge_pos, function(x) any(merge_pos[[l]] %in% x))
          new_merge_pos <- list.append(new_merge_pos, sort(unique(unlist(merge_pos[which(id)]))))
        }else{
          id_b <- which(sapply(merge_pos, function(x) all(merge_pos[[l]] %in% x)))
          id_b <- id_b[which.max(sapply(merge_pos[id_b], length))]
          new_merge_pos <- list.append(new_merge_pos, sort(merge_pos[[id_b]]))
        }
        
      }
      
      all_merg <- all(!duplicated(unlist(new_merge_pos)))
      merge_pos <- new_merge_pos
      merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
      merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
      merge_pos <- lapply( merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
      print(length(merge_pos))
      new_merge_pos <- list() 
      
    }
  }
  
  tmp <-  lapply(merge_pos, function(x) data.frame(npath = length(x), path = paste0(tmp_info$path[x], collapse = '-and-'), 
                                                   mean_Zstat = mean(tmp_info$Zstat[x]), sd_Zstat = sd(tmp_info$Zstat[x]), 
                                                   highest_Zstat = tmp_info$Zstat[x][which.max(abs(tmp_info$Zstat[x]))], 
                                                   mean_ngenes = mean(tmp_info$ngenes_tscore[x]), sd_ngenes = sd(tmp_info$ngenes_tscore[x]), 
                                                   tissue = tissues[i]))
  
  tmp_path <-  rbind(tmp_path,do.call(rbind, tmp))
  tissue_spec_pathGO[[i]] <- tmp_path
  
  tissue_spec_pathGO[[i]]$comp_sign <- NA
  # for each loci, find groups that are significantly different
  for(l in 1:nrow(tissue_spec_pathGO[[i]])){
    paths <- strsplit(tissue_spec_pathGO[[i]]$path[l], split = '-and-')[[1]]
    tmp_gr <- sapply(unique(tmp_feat$comp[tmp_feat$feat %in% paths & tmp_feat$pval_corr <= 0.05]), function(x) strsplit(x, split = '_vs_all')[[1]][1])
    tissue_spec_pathGO[[i]]$comp_sign[l] <-  paste0(tmp_gr, collapse = ',') 
  }
}

tissue_spec_pathGO <- do.call(rbind, tissue_spec_pathGO)
# save
write.table(x = tissue_spec_pathGO, file = sprintf('%s%s_%s_cluster%s_summary_path_GO_tissueSpec.txt',outFold, type_data, type_input, type_cluster), 
            col.names = T, row.names = F, sep = '\t', quote = F)


### plot distributions ###
gr_tot <- unique(pathGO_feat$comp)
gr_tot <- sapply(gr_tot, function(x) strsplit(x, split = '_vs_all')[[1]][1])

df <- data.frame(npath = c(), tissue = c(), comp = c(), path_group = c())
for(j in 1:length(tissues)){
  for(i in 1:length(gr_tot)){
    tmp <- tissue_spec_pathGO[grepl(gr_tot[i], tissue_spec_pathGO$comp) & tissue_spec_pathGO$tissue %in% tissues[j],]
    df <- rbind(df, data.frame(npath = tmp$npath, tissue = rep(tissues[j], nrow(tmp)), comp = rep(gr_tot[i], nrow(tmp)), path_group = tmp$path))
  }
}

# too complicated, do not show legend but class
# df$path_group_name <- df$path_group
# id_m <- which(df$npath > 1)
# df$path_group_name <- factor(df$path_group_name, levels = unique(df$path_group_name))

# same plot but n. of loci and not genes
df_pathgroup <- data.frame(tissue = unlist(lapply(tissues, function(x) rep(x, length(gr_tot)))), comp = rep(gr_tot, length(tissues)), 
                           ngroup = as.vector(sapply( tissues, function(x) sapply(gr_tot, function(y) sum(df$comp == y & df$tissue == x)))))

df_pathgroup$comp <- factor(df_pathgroup$comp, levels = unname(gr_tot))
df_pathgroup$tissue <- factor(df_pathgroup$tissue, levels =  tissues)

df_tot <- data.frame(tissue = unlist(lapply(tissues, function(x) rep(x, length(gr_tot)))), comp = rep(gr_tot, length(tissues)),
                     npath = as.vector(sapply(tissues, function(x) sapply(gr_tot, function(y) sum(df$npath[df$tissue == x & df$comp == y])))))
df_tot$comp <- factor(df_tot$comp, levels = unname(gr_tot))
df_tot$tissue <- factor(df_tot$tissue, levels =  tissues)
newcolours <-  color_tissues$color[match(tissues, color_tissues$tissues)]
color_gr <- pal_d3("category10")(length(gr_tot))

pl1 <- ggplot(data = subset(df_tot, tissue != 'All'), aes(x = tissue, y = npath, fill = comp))+
  geom_bar(alpha = 0.9, width = 0.9, stat = 'identity', position = position_dodge())+
  # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. pathways\ncluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'bottom', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 8), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours))+
  guides(fill=guide_legend(nrow=1))+
  scale_fill_manual(values = color_gr)+
  # scale_fill_d3()+
  coord_flip() 

# pl1 <- ggplot_gtable(ggplot_build(pl1))
# stripr <- which(grepl('strip-t', pl1$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl1$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

pl2 <- ggplot(data = subset(df_pathgroup, tissue != 'All'), aes(x = tissue, y = ngroup, fill = comp))+
  geom_bar(alpha = 0.9, width = 0.9, stat = 'identity',  position = position_dodge())+
  # facet_wrap(.~comp, ncol = length(gr_tot), scales = 'free_x')+
  ylab('n. pathway macro\ncluster relevant')+ 
  theme_bw()+ 
  theme(legend.position = 'none', legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6), legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_text(colour = newcolours))+
  # guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values = color_gr)+
  coord_flip() 

# pl2 <- ggplot_gtable(ggplot_build(pl2))
# stripr <- which(grepl('strip-t', pl2$layout$name))
# fills <- color_gr
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', pl2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   pl2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }

tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=2, nrow=1, align = 'h', common.legend = T, legend = 'bottom')
# tot_pl <- ggarrange(plotlist = list(pl1, pl2), ncol=1, nrow=2, align = 'v')
ggsave(filename =  sprintf('%scluster_npath_npathgroup_path_GO_pertissues_%s_%s.png', outFold, type_data, type_input), plot = tot_pl, width = 7, height = 4.5, dpi = 500)
ggsave(filename = sprintf('%scluster_npath_npathgroup_path_GO_pertissues_%s_%s.pdf', outFold, type_data, type_input), plot = tot_pl, width = 7, height = 4.5, dpi = 500, compress = F)

