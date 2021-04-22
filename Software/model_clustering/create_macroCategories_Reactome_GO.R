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

### class of 'Immune system' ###
class_names_immune <- df_R$Child[df_R$Parent == 'R-HSA-168256']
class_names_immune <- list_path[list_path$V1 %in% class_names_immune, ]
distMatrix_immune <- shortest.paths(gr_R, v=class_names_immune$V1, to=Reactome_macro[[15]]$subclasses$V1, mode = 'out')
distMatrix_immune_original <- distMatrix_immune
for(j in 1:ncol(distMatrix_immune)){
  if(sum(distMatrix_immune[, j] != Inf)>1){
    distMatrix_immune[distMatrix_immune[, j]!=min(distMatrix_immune[, j]), j] <- Inf
  }
}
group_class <- apply(distMatrix_immune, 1, function(x) names(which(x!=Inf)))
Reactome_macro_immune <- list()
for(i in 1:length(group_class)){
  Reactome_macro_immune[[i]] <- list(macro = class_names_immune[i, 1:2], 
                              subclasses = list_path[list_path$V1 %in% group_class[[i]], 1:2])
  Reactome_macro_immune[[i]]$subclasses$shortestPath <- distMatrix_immune_original[i,match(Reactome_macro_immune[[i]]$subclasses$V1,colnames(distMatrix_immune_original))]
}
save(Reactome_macro_immune, file = 'ReactomePathways_macro_Immune_2021.RData')


