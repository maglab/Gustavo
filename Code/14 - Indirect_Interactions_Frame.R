# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., GenArdArcFrm)
# All -All (e.g., Grp$All)
# Arc - Ageing-related Communities (e.g., ArdArcFrm)
# Ard - Ageing-related Diseases (e.g., GenArdArcFrm)
# Arr - Array (e.g., GenArr)
# Cod - Code (e.g., dplyr::rename(Cod=Code))
# Cox - Coexpression (e.g., 'COX90')
# Frm - Frame (e.g., GenArd_IndFrm)
# Gen - Gene (e.g., GenArdArcFrm)
# Grp - Graph (e.g., Grp$All)
# Ind - Indirect (e.g., GenArd_IndFrm)
# Len - Length (e.g., GenLen)
# Lng - Long (i.e., Matrix converted into a Two/Three colums dataframe) (e.g., GenArd_IndFrmLng)
# Nod - Node (e.g., Nod$Ard)
# Ntw - Network (e.g., NtwArr)
# Pth - Path (e.g., Pth = paste(...))
# Qry - Query (e.g., ArdQry)
# Var - Variable (e.g., rownames_to_column(var = "rowname"))

### LIBRARIES ##################################################################

library(igraph)
library(rje)
library(dplyr)
library(tidyr)
library(tidyverse)

GenArdArcFrm = readRDS("maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)

ArdArcFrm = GenArdArcFrm %>% select(Ard,Arc) %>% unique()

####################################################################
### COMPUTE INDIRECT INTERACTIONS ##################################
####################################################################

NtwArr = c('PPI','COX90','COX95','KEGG')
NtwLen = length(NtwArr)
for(i in 1:NtwLen){
  
  NtwQry = NtwArr[i]
  
  Pth = paste("maglab/Gustavo/Data/Generated/",NtwQry,"/",sep="")
  
  Nod = paste(Pth,"Lists/","NodesList.rds",sep="") %>% readRDS()
  Grp = paste(Pth,"Lists/","GraphsList.rds",sep="") %>% readRDS()

  QryGrp = Grp$All %>% delete_vertices(Nod$Ard)
  GenAdjFrm = QryGrp %>% as_adjacency_matrix %>% as.matrix() %>% as.data.frame()
  
  GenArr = row.names(GenAdjFrm)
  GenLen = length(GenArr)
  
  ### GEN-PHN ########################################################
  
  ArdArr = GenArdArcFrm$Ard %>% unique()
  ArdLen = length(ArdArr)
  GenArd_IndFrm = matrix(0,GenLen,ArdLen) %>% as.data.frame()
  colnames(GenArd_IndFrm) = ArdArr
  row.names(GenArd_IndFrm) = row.names(GenAdjFrm)
  
  dim(GenArd_IndFrm)
  dim(GenAdjFrm)
  length(GenArr)
  
  # Get the number of association of genes to each phenotype
  for(j in 1:ArdLen){
    paste(i,1,j) %>% print()
    ArdQry = ArdArr[j] 
    GenArd = GenArdArcFrm %>% filter(Ard %in% ArdQry) %>% pull(Gen)
    GenArd_IndFrm[, ArdQry] = GenAdjFrm[, GenArd, drop = FALSE] %>% rowSums()
  }
  
  # LONG FORM OF INTERACTOR FRAME
  GenArd_IndFrmLng <- GenArd_IndFrm %>%
    rownames_to_column(var = "rowname") %>%   # Convert row names to a column
    pivot_longer(
      cols = -rowname,                        # Exclude the 'rowname' column from pivoting
      names_to = "variable",                  # Name for the column with original column names
      values_to = "value"                     # Name for the column with values
    ) %>% 
    filter(value>0) %>%
    rename(Gen=rowname,Ard=variable,Val=value) %>%
    merge(ArdArcFrm) %>% 
    select(Gen, Ard, Arc) %>%
    unique()
    
  GenArd_IndFrmLng %>% write.csv(paste(Pth,"Proximity_and_Distances/Indirect_Interactions_Long.csv",sep=""))
  
}

