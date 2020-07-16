# load enrichment result for each tissue and make plot
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(apcluster))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="find enrichment of CAD associated element with UKBB gene-pheno or path-pheno phenotypes")
parser$add_argument("--phenoFold", type = "character", help = "Folder with results phenotype annotation UKBB")
parser$add_argument("--inputFold_CAD", type = "character", nargs = '*', help = "Folder with results from CAD (1 for each tissue)")
parser$add_argument("--pval_FDR_CAD", type = "double", default = 0.05, help = "pval threshold to filter the genes and pathways (after BH correction) in CAD pheno")
parser$add_argument("--pval_FDR_CADrel", type = "double", default = 0.05,  help = "pval threshold to filter the genes and pathways (after BH correction) in UKBB phenotypes")
parser$add_argument("--pval_FDR_fisher", type = "double", default = 0.05, help = "pval threshold to filter the phenotypes (after BH correction) for enrichment")
parser$add_argument("--tissues_name", type = "character",nargs = '*', help = "tissue considered")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
phenoFold <- args$phenoFold
tissues_name <- args$tissues_name
inputFold_CAD <- args$inputFold_CAD
pval_FDR_fisher <- args$pval_FDR_fisher
pval_FDR_CAD <- args$pval_FDR_CAD
pval_FDR_CADrel <- args$pval_FDR_CADrel
outFold <- args$outFold

# #########################################################################################################################
# phenoFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/'
# tissues_name <- c('Adipose_Subcutaneous' ,'Adipose_Visceral_Omentum' ,'Adrenal_Gland' ,'Artery_Aorta', 'Artery_Coronary' ,'Colon_Sigmoid' ,'Colon_Transverse' ,'Liver')
# inputFold_CAD <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/',tissues_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/'
# pval_FDR_fisher <- 0.05
# pval_FDR_CAD <- 0.05
# pval_FDR_CADrel <- 0.05
# #########################################################################################################################

df_tscore <- list()
df_pathR <- list()
df_pathGO <- list()

# load results
for(i in 1:length(tissues_name)){
  
  t <- tissues_name[i]
  print(t)
  print(sprintf('%s%s_enrichment_significant_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.RData', inputFold_CAD[i], t, pval_FDR_CAD, pval_FDR_CADrel))
  tmp <- get(load(sprintf('%s%s_enrichment_significant_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.RData', inputFold_CAD[i], t, pval_FDR_CAD, pval_FDR_CADrel)))
  
  df_tscore[[i]] <- cbind(tmp$tscore, data.frame(tissue = t, stringsAsFactors = F))
  df_pathR[[i]] <- cbind(tmp$pathScore_reactome, data.frame(tissue = t, stringsAsFactors = F))
  df_pathGO[[i]] <- cbind(tmp$pathScore_GO, data.frame(tissue = t, stringsAsFactors = F))
  
}

df_tscore <- do.call(rbind,df_tscore)
df_pathR <- do.call(rbind,df_pathR)
df_pathGO <- do.call(rbind,df_pathGO)

# create matrices for heatmap
# pval
mat_tscore_pval <- matrix(nrow = length(tissues_name), ncol = nrow(df_tscore)/length(tissues_name), -log10(df_tscore$fisher_pval), byrow = T)
colnames(mat_tscore_pval) <- df_tscore$pheno[1:ncol(mat_tscore_pval)]
rownames(mat_tscore_pval) <- tissues_name
mat_tscore_pvalcorr <- matrix(nrow = length(tissues_name), ncol = nrow(df_tscore)/length(tissues_name), df_tscore$fisher_pval_BHcorr, byrow = T)
id_excl <- mat_tscore_pvalcorr>pval_FDR_fisher
mat_tscore_pval[id_excl] <- NA
mat_tscore_pval <- mat_tscore_pval[, colSums(is.na(mat_tscore_pval))<nrow(mat_tscore_pval)]
mat_tscore_pval <- mat_tscore_pval[rowSums(is.na(mat_tscore_pval))<ncol(mat_tscore_pval),]
# OR
mat_tscore_OR <- matrix(nrow = length(tissues_name), ncol = nrow(df_tscore)/length(tissues_name), df_tscore$fisher_OR, byrow = T)
colnames(mat_tscore_OR) <- df_tscore$pheno[1:ncol(mat_tscore_OR)]
rownames(mat_tscore_OR) <- tissues_name
mat_tscore_OR[id_excl] <- NA
mat_tscore_OR <- mat_tscore_OR[, colSums(is.na(mat_tscore_OR))<nrow(mat_tscore_OR)]
mat_tscore_OR <- mat_tscore_OR[rowSums(is.na(mat_tscore_OR))<ncol(mat_tscore_OR),]

#pval
mat_pathR_pval <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathR)/length(tissues_name), -log10(df_pathR$fisher_pval), byrow = T)
colnames(mat_pathR_pval) <- df_pathR$pheno[1:ncol(mat_pathR_pval)]
rownames(mat_pathR_pval) <- tissues_name
mat_pathR_pvalcorr <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathR)/length(tissues_name), df_pathR$fisher_pval_BHcorr, byrow = T)
id_excl <- mat_pathR_pvalcorr>pval_FDR_fisher
mat_pathR_pval[id_excl] <- NA
mat_pathR_pval <- mat_pathR_pval[, colSums(is.na(mat_pathR_pval))<nrow(mat_pathR_pval)]
mat_pathR_pval <- mat_pathR_pval[rowSums(is.na(mat_pathR_pval))<ncol(mat_pathR_pval),]
# OR
mat_pathR_OR <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathR)/length(tissues_name), df_pathR$fisher_OR, byrow = T)
colnames(mat_pathR_OR) <- df_tscore$pheno[1:ncol(mat_pathR_OR)]
rownames(mat_pathR_OR) <- tissues_name
mat_pathR_OR[id_excl] <- NA
mat_pathR_OR <- mat_pathR_OR[, colSums(is.na(mat_pathR_OR))<nrow(mat_pathR_OR)]
mat_pathR_OR <- mat_pathR_OR[rowSums(is.na(mat_pathR_OR))<ncol(mat_pathR_OR),]

# pval
mat_pathGO_pval <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathGO)/length(tissues_name), -log10(df_pathGO$fisher_pval), byrow = T)
colnames(mat_pathGO_pval) <- df_pathGO$pheno[1:ncol(mat_pathGO_pval)]
rownames(mat_pathGO_pval) <- tissues_name
mat_pathGO_pvalcorr <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathGO)/length(tissues_name), df_pathGO$fisher_pval_BHcorr, byrow = T)
id_excl <- mat_pathGO_pvalcorr > pval_FDR_fisher
mat_pathGO_pval[id_excl] <- NA
mat_pathGO_pval <- mat_pathGO_pval[, colSums(is.na(mat_pathGO_pval))<nrow(mat_pathGO_pval)]
mat_pathGO_pval <- mat_pathGO_pval[rowSums(is.na(mat_pathGO_pval))<ncol(mat_pathGO_pval),]
# OR
mat_pathGO_OR <- matrix(nrow = length(tissues_name), ncol = nrow(df_pathGO)/length(tissues_name), df_pathGO$fisher_OR, byrow = T)
colnames(mat_pathGO_OR) <- df_tscore$pheno[1:ncol(mat_pathGO_OR)]
rownames(mat_pathGO_OR) <- tissues_name
mat_pathGO_OR[id_excl] <- NA
mat_pathGO_OR <- mat_pathGO_OR[, colSums(is.na(mat_pathGO_OR))<nrow(mat_pathGO_OR)]
mat_pathGO_OR <- mat_pathGO_OR[rowSums(is.na(mat_pathGO_OR))<ncol(mat_pathGO_OR),]

pheno_ann <- read.delim(sprintf('%scolor_pheno_type_UKBB.txt', phenoFold), header = T, stringsAsFactors = F)
pheno_info <- df_tscore[1:(nrow(df_pathR)/length(tissues_name)), c('names_field', 'pheno','pheno_type')]

color_tissues <- read.table('/psycl/g/mpsziller/lucia/color_tissues.txt', h=T, stringsAsFactors = F)
color_tissues <- color_tissues[match(tissues_name,color_tissues$tissue),]

# save total df
write.table(x = df_tscore, file = sprintf('%s/tscore_enrichment_alltissues_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.txt', outFold, pval_FDR_CAD, pval_FDR_CADrel), 
            col.names = T, row.names = F, sep = '\t', quote = F)  
write.table(x = df_pathR, file = sprintf('%s/pathScore_Reactome_enrichment_alltissues_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.txt', outFold, pval_FDR_CAD, pval_FDR_CADrel), 
            col.names = T, row.names = F, sep = '\t', quote = F)  
write.table(x = df_pathGO, file = sprintf('%s/pathScore_GO_enrichment_alltissues_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.txt', outFold, pval_FDR_CAD, pval_FDR_CADrel), 
            col.names = T, row.names = F, sep = '\t', quote = F)  


############# make plot ################
save_pheatmap_png <- function(x, filename, width=10, height=10, res = 200) {
  png(filename, width = width, height = height, res = res, units = 'in')
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
# used to make heatmaps 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_90",
  ns = asNamespace("pheatmap")
)

plot_heatmap_split <- function(type_mat, mat_pval, pheno_ann, pheno_info, color_tissues, height_pl = 17, width_pl=11, cap_val = 10, split_pheno = T, type_input = 'pvalue'){
  
  coul <- colorRampPalette(brewer.pal(8, "BuPu"))(25)
  coul[1] <- 'white'
  title_sub <- ifelse(type_input == 'pvalue', '-log10(pval)','OR' )
  
  if(split_pheno){
    
    # split blood, ICD10 and rest
    id_blood <- colnames(mat_pval) %in% pheno_info$pheno[pheno_info$pheno_type == 'Blood_count'] 
    mat_pval_blood <- mat_pval[, which(id_blood)]
    
    id_icd10 <- colnames(mat_pval) %in% pheno_info$pheno[grepl('ICD10', pheno_info$pheno_type)] 
    mat_pval_icd10 <- mat_pval[, which(id_icd10)]
    
    id_blood_bio <- colnames(mat_pval) %in% pheno_info$pheno[pheno_info$pheno_type == 'Blood_biochemistry'] 
    mat_pval_bloodbio <- mat_pval[, which(id_blood_bio)]
    
    mat_pval_rest <- mat_pval[, which(!id_blood & !id_icd10 & !id_blood_bio)]
    
    ########################
    ##### blood count ######
    ########################
    
    if(ncol(mat_pval_blood)>0){
      tmp_mat <- as.matrix(t(mat_pval_blood))
      tmp_mat[is.na(tmp_mat)] <- 0
      val <- min(cap_val,round(max(tmp_mat)))
      mat_breaks <- seq(0, val, length.out = 25)
      pheno_tmp <- pheno_info[match(rownames(tmp_mat), pheno_info$pheno),] 
      
      # Data frame with column annotations.
      mat_row <- data.frame(pheno = pheno_tmp$pheno_type)
      rownames(mat_row) <- pheno_tmp$names_field
      
      mat_colors <- list(pheno = pheno_ann$color)
      names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
      
      tmp_mat_capped <- tmp_mat
      tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
      rownames(tmp_mat_capped) <- pheno_tmp$names_field
      
      mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat_capped), color_tissues$tissue)])
      rownames(mat_col) <- colnames(tmp_mat_capped)
      
      mat_colors_t <- list(tissue_type = unique(color_tissues$color))
      names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
      
      new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
      
      
      hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, 
                        annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                        cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                        annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 9, fontsize_col = 10, fontsize = 12, 
                        main = sprintf('%s enrichment phenotypes and CAD-HARD\n(%s)', title_sub, type_mat),
                        breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 11)
      
      
      save_pheatmap_png(hm_pl, sprintf("%s/%s_enrichment%s_Blood_count_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
      save_pheatmap_pdf(hm_pl, sprintf("%s/%s_enrichment%s_Blood_count_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
    }
    
    ########################
    ##### blood count ######
    ########################
    
    if(ncol(mat_pval_bloodbio)>0){
      tmp_mat <- as.matrix(t(mat_pval_bloodbio))
      tmp_mat[is.na(tmp_mat)] <- 0
      val <- min(cap_val,round(max(tmp_mat)))
      mat_breaks <- seq(0, val, length.out = 25)
      pheno_tmp <- pheno_info[match(rownames(tmp_mat), pheno_info$pheno),] 
      
      # Data frame with column annotations.
      mat_row <- data.frame(pheno = pheno_tmp$pheno_type)
      rownames(mat_row) <- pheno_tmp$names_field
      
      mat_colors <- list(pheno = pheno_ann$color)
      names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
      
      tmp_mat_capped <- tmp_mat
      tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
      rownames(tmp_mat_capped) <- pheno_tmp$names_field
      
      mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat_capped), color_tissues$tissue)])
      rownames(mat_col) <- colnames(tmp_mat_capped)
      
      mat_colors_t <- list(tissue_type = unique(color_tissues$color))
      names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
      
      new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
      
      
      hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, 
                        annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                        cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                        annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 9, fontsize_col = 10, fontsize = 12, 
                        main = sprintf('%s enrichment phenotypes and CAD-HARD\n(%s)', title_sub, type_mat),
                        breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 11)
      
      
      save_pheatmap_png(hm_pl, sprintf("%s/%s_enrichment%s_Blood_biochemistry_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
      save_pheatmap_pdf(hm_pl, sprintf("%s/%s_enrichment%s_Blood_biochemistry_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
    }
    
    ##################
    ##### icd10 ######
    ##################
    
    if(ncol(mat_pval_icd10)>0){
      tmp_mat <- as.matrix(t(mat_pval_icd10))
      tmp_mat[is.na(tmp_mat)] <- 0
      val <- min(cap_val,round(max(tmp_mat)))
      mat_breaks <- seq(0, val, length.out = 25)
      pheno_tmp <- pheno_info[match(rownames(tmp_mat), pheno_info$pheno),] 
      
      # Data frame with column annotations.
      mat_row <- data.frame(pheno = pheno_tmp$pheno_type)
      rownames(mat_row) <- pheno_tmp$names_field
      
      mat_colors <- list(pheno = pheno_ann$color)
      names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
      
      tmp_mat_capped <- tmp_mat
      tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
      rownames(tmp_mat_capped) <- pheno_tmp$names_field
      
      mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat_capped), color_tissues$tissue)])
      rownames(mat_col) <- colnames(tmp_mat_capped)
      
      mat_colors_t <- list(tissue_type = unique(color_tissues$color))
      names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
      
      new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
      
      hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, 
                        annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                        cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                        annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 9, fontsize_col = 10, fontsize = 12, 
                        main = sprintf('%s enrichment phenotypes and CAD-HARD\n(%s)', title_sub, type_mat),
                        breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 11)
      
      
      save_pheatmap_png(hm_pl, sprintf("%s/%s_enrichment%s_icd10_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
      save_pheatmap_pdf(hm_pl, sprintf("%s/%s_enrichment%s_icd10_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl-5, width =width_pl)
    }
  }else{
    mat_pval_rest <- mat_pval
  }
  
  #################
  ##### rest ######
  #################
  
  tmp_mat <- as.matrix(t(mat_pval_rest))
  tmp_mat[is.na(tmp_mat)] <- 0
  val <- min(cap_val,round(max(tmp_mat)))
  mat_breaks <- seq(0, val, length.out = 25)
  pheno_tmp <- pheno_info[match(rownames(tmp_mat), pheno_info$pheno),] 
  
  # Data frame with column annotations.
  mat_row <- data.frame(pheno = pheno_tmp$pheno_type)
  rownames(mat_row) <- pheno_tmp$names_field
  
  mat_colors <- list(pheno = pheno_ann$color)
  names(mat_colors$pheno) <- unique(pheno_ann$pheno_type)
  
  tmp_mat_capped <- tmp_mat
  tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
  rownames(tmp_mat_capped) <- pheno_tmp$names_field
  
  mat_col <- data.frame(tissue_type = color_tissues$type[match(colnames(tmp_mat_capped), color_tissues$tissue)])
  rownames(mat_col) <- colnames(tmp_mat_capped)
  
  mat_colors_t <- list(tissue_type = unique(color_tissues$color))
  names(mat_colors_t$tissue_type) <- unique(color_tissues$type)
  
  new_color = list(pheno = mat_colors[[1]], tissue_type = mat_colors_t[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, 
                    annotation_row = mat_row,  annotation_col = mat_col, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = 'darkgrey', 
                    annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 10, fontsize_col = 10, fontsize = 12, 
                    main = sprintf('%s enrichment phenotypes and CAD-HARD\n(%s)', title_sub, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0, cellwidth = 10)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s/%s_enrichment%s_rest_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl, width =width_pl, res = 200)
  save_pheatmap_pdf(hm_pl, sprintf("%s/%s_enrichment%s_rest_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, type_mat, type_input, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl, width =width_pl)
  
}

# pvalue
plot_heatmap_split(type_mat = 'tscore', mat_pval = mat_tscore_pval, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val=15,
                   color_tissues = color_tissues, height_pl = 15, width_pl = 14)

plot_heatmap_split(type_mat = 'pathScore_Reactome', mat_pval = mat_pathR_pval, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val=15, 
                        color_tissues = color_tissues, height_pl = 10, width_pl = 14, split_pheno = F)

plot_heatmap_split(type_mat = 'pathScore_GO', mat_pval = mat_pathGO_pval, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val=15, 
                        color_tissues = color_tissues, height_pl = 12, width_pl = 14)

# OR
plot_heatmap_split(type_mat = 'tscore', mat_pval = mat_tscore_OR, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val = 10,
                   color_tissues = color_tissues, height_pl = 15, width_pl = 14, type_input = 'OR')

plot_heatmap_split(type_mat = 'pathScore_Reactome', mat_pval = mat_pathR_OR, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val = 10, 
                   color_tissues = color_tissues, height_pl = 10, width_pl = 14, split_pheno = F,  type_input = 'OR')

plot_heatmap_split(type_mat = 'pathScore_GO', mat_pval = mat_pathGO_OR, pheno_ann = pheno_ann, pheno_info = pheno_info, cap_val = 10,
                   color_tissues = color_tissues, height_pl = 12, width_pl = 14,  type_input = 'OR')

