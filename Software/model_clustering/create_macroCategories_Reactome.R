library(igraph)

setwd('/psycl/g/mpsziller/lucia/refData/')

df_R <- read.table('ReactomePathwaysRelation_Human.txt', h=F, stringsAsFactors = F, sep = '\t')
colnames(df_R) <- c('Parent', 'Child')
tot_names <- unique(df_R$Parent, df_R$Child)
class_names <- tot_names[tot_names %in% df_R$Parent & !tot_names %in% df_R$Child]
list_path <- read.delim('ReactomePathways_Human.txt', h=F, stringsAsFactors = F, sep = '\t')
class_names <- list_path[list_path$V1 %in% class_names, ]

gr_R = graph_from_data_frame(df_R, directed = TRUE, vertices = NULL)
distMatrix <- shortest.paths(gr_R, v=class_names$V1, to=V(gr_R))
distMatrix_original <- distMatrix
# some paths can have multiple macro, put Inf to the ones that are higher
id_edit_manually <- which(apply(distMatrix, 2, function(x) sum(x == min(x)))>1)
# distMatrix['R-HSA-69306', 'R-HSA-113510'] <- Inf
# distMatrix['R-HSA-1474165', 'R-HSA-1500620'] <- Inf

for(j in 1:ncol(distMatrix)){
  if(sum(distMatrix[, j] != Inf)>1){
    distMatrix[distMatrix[, j]!=min(distMatrix[, j]), j] <- Inf
  }
}
group_class <- apply(distMatrix, 1, function(x) names(which(x!=Inf)))
Reactome_macro <- list()
for(i in 1:length(group_class)){
  Reactome_macro[[i]] <- list(macro = class_names[i, 1:2], 
                              subclasses = list_path[list_path$V1 %in% group_class[[i]], 1:2])
  Reactome_macro[[i]]$subclasses$shortestPath <- distMatrix_original[i,match(Reactome_macro[[i]]$subclasses$V1,colnames(distMatrix_original))]
}
save(Reactome_macro, file = 'ReactomePathways_macro_2021.RData')
# remove macro with less than 10 elements:
len_macro <- sapply(Reactome_macro, function(x) nrow(x[[2]]))
Reactome_macro_red <- Reactome_macro[-which(len_macro<10)]
save(Reactome_macro_red, file = 'ReactomePathways_macro_filtN10_2021.RData')

# extract class
extract_class <- function(code_macro, df_parent_child, Reactome_macro, graph_reactome){
  
  ### class of 'xx' ###
  class_names_sub <- df_parent_child$Child[df_parent_child$Parent == code_macro]
  id <- which(sapply(Reactome_macro, function(x) x$macro$V1) == code_macro)
  class_names_sub <- list_path[list_path$V1 %in% class_names_sub, ]
  distMatrix_sub <- shortest.paths(graph_reactome, v=class_names_sub$V1, to=Reactome_macro[[id]]$subclasses$V1, mode = 'out')
  distMatrix_sub_original <- distMatrix_sub
  for(j in 1:ncol(distMatrix_sub)){
    if(sum(distMatrix_sub[, j] != Inf)>1){
      distMatrix_sub[distMatrix_sub[, j]!=min(distMatrix_sub[, j]), j] <- Inf
    }
  }
  group_class <- apply(distMatrix_sub, 1, function(x) names(which(x!=Inf)))
  Reactome_macro_sub <- list()
  for(i in 1:length(group_class)){
    Reactome_macro_sub[[i]] <- list(macro = class_names_sub[i, 1:2], 
                                       subclasses = list_path[list_path$V1 %in% group_class[[i]], 1:2])
    Reactome_macro_sub[[i]]$subclasses$shortestPath <- distMatrix_sub_original[i,match(Reactome_macro_sub[[i]]$subclasses$V1,colnames(distMatrix_sub_original))]
  }
  return(Reactome_macro_sub)
}

### class of 'Immune System' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-168256',
                                       graph_reactome = gr_R, 
                                       df_parent_child = df_R, 
                                       Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_Immune_2021.RData')


### class of 'Neuronal System' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-112316',
                                       graph_reactome = gr_R, 
                                       df_parent_child = df_R, 
                                       Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_NeuronalSystem_2021.RData')

### class of 'Cell Cycle System' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-1640170',
                                       graph_reactome = gr_R, 
                                       df_parent_child = df_R, 
                                       Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_CellCycle_2021.RData')

### class of 'Gene Expression Transcription' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-74160',
                                       graph_reactome = gr_R, 
                                       df_parent_child = df_R, 
                                       Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_GeneExpr_2021.RData')

### class of 'Programmed Cell Death' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-5357801',
                                    graph_reactome = gr_R, 
                                    df_parent_child = df_R, 
                                    Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_ProgramCellDeath_2021.RData')


### class of 'Metabolism of Protein' ###
Reactome_macro_sub <- extract_class(code_macro = 'R-HSA-8953854',
                                    graph_reactome = gr_R, 
                                    df_parent_child = df_R, 
                                    Reactome_macro = Reactome_macro)
save(Reactome_macro_sub, file = 'ReactomePathways_macro_MetabolismProtein_2021.RData')



