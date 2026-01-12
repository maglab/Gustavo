library(igraph)
library(tidyverse)
#install.packages('tidyverse')

# DICTIONARY ##################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Asc - Ascendants
# Cat - Category
# Cl1 - Cluster 1
# Cl2 - Cluster 2
# Cod - Code
# Cnt - Count
# Dis - Disease
# Frm - Frame
# Fcn - Function
# Grp - Graph
# Hrc - Hierarchy
# Len - Length
# Mng - Meaning
# Nod - Node
# Pth - Path
# Qry - Query
# Sel - Selected
# Srt - Sorted
# Top - Top
# Vec - Vector
# Vrt - Vertice

# START ########################################################################

# Load Hierarchy data of UK Bioank self-reported diseases
DisCat = read.csv("Data/Retrieved/Genes_and_diseases/Showcase/Ukb_DiseaseHierarchy.csv")

GrpDis = graph_from_data_frame(DisCat[,c("Node", "Parent")])

NodVec = V(GrpDis)$name

HrcFrm = data.frame()
for(i in 1:length(NodVec)){
  print(i)
  QryNod = NodVec[i]
  if (QryNod != "Top"){
    HrClustert = AgeAscDis_Fcn(QryNod, GrpDis, DisCat) # Cluster based on root ascendant (ARC) of ARDs
    sHrcFrm = HrClustert$NodFrm
    HrcFrm = rbind(HrcFrm, sHrcFrm[nrow(sHrcFrm),])
  }
}

HrcFrm = unique(HrcFrm)
HrcFrm$X = NULL

HrcFrm = HrcFrm %>% rename(Node=Nod,
                  Coding=Cod,
                  Meaning=Mng,
                  Selected=Sel,
                  RootMeaning=RootMng,
                  RootNode=RootNod)

HrcFrm$Node = paste("D",HrcFrm$Node,sep="")

write.csv(HrcFrm, "Data/Retrieved/Genes_and_diseases/Diseases/Ukb_DiseaseRoots.csv")



#################################################################################################
### FUNCTIONS ###################################################################################
#################################################################################################


AgeAscDis_Fcn <- function(QryNod, GrpDis, DisCat){
  
  ShrtPth = shortest_paths(GrpDis,
                                  from = QryNod,
                                  to = "Top",
                                  mode = c("out", "all", "in")
  )
  
  AllAscNod = as_ids(ShrtPth$vpath[[1]])
  
  if(!(setequal(AllAscNod,"Top"))){
    RootNodIdx = length(AllAscNod) - 1
    AscFrm = DisCat %>% filter(Node %in% AllAscNod) %>% dplyr::rename(Cod = Coding, Mng = Meaning, Sel = Selectable, Nod = Node, Parent  = Parent, Cluster = Ageing)
    AscFrm$RootMng = DisCat %>% filter(Node == AllAscNod[RootNodIdx]) %>% pull("Meaning")
    AscFrm$RootNod = DisCat %>% filter(Node == AllAscNod[RootNodIdx]) %>% pull("Node")
  
    AllAscNod = rev(AllAscNod)
    AllAscLen = length(AllAscNod) - 1
    DepthFrm = data.frame("Nod" = AllAscNod, "Depth" = c(0:AllAscLen))
    
    NodFrm = merge(AscFrm, DepthFrm, by = "Nod") %>% arrange(Depth)
    
    AgeAscNod = list()
    AgeAscNod$All = AgeAscNod_Fcn(NodFrm, AllAscNod, Cluster = c(1,2))
    AgeAscNod$Cl1 = AgeAscNod_Fcn(NodFrm, AllAscNod, Cluster = c(1))
    AgeAscNod$Cl2 = AgeAscNod_Fcn(NodFrm, AllAscNod, Cluster = c(2))
    
    AllAscMng = c("Top" , NodFrm$Mng)
    
    AgeAscMng = list()
    AgeAscMng$All = AgeAscMng_Fcn(NodFrm, AllAscMng, Cluster = c(1,2))
    AgeAscMng$Cl1 = AgeAscMng_Fcn(NodFrm, AllAscMng, Cluster = c(1))
    AgeAscMng$Cl2 = AgeAscMng_Fcn(NodFrm, AllAscMng, Cluster = c(2))
    
    IncAgeAscCluster = c("Top",NodFrm$Cluster)
    AgeAscCluster = IncAgeAscCluster[1:(length(IncAgeAscCluster)-1)]
    
    IncPthAgeClusterCnt = list()
    IncPthAgeClusterCnt$All = sum(IncAgeAscCluster %in% c(1,2))
    IncPthAgeClusterCnt$Cl1 = sum(IncAgeAscCluster %in% c(1))
    IncPthAgeClusterCnt$Cl2 = sum(IncAgeAscCluster %in% c(2))
    
    PthAgeClusterCnt = list()
    PthAgeClusterCnt$All = sum(AgeAscCluster %in% c(1,2))
    PthAgeClusterCnt$Cl1 = sum(AgeAscCluster %in% c(1))
    PthAgeClusterCnt$Cl2 = sum(AgeAscCluster %in% c(2))
    
  NodLst = list("NodFrm"=NodFrm, "AgeAscNod"=AgeAscNod, "AllAscNod"=AllAscNod, "MskAllAscNod"=AgeAscNod, 
                "AllAscMng" = AllAscMng, "AgeAscMng" = AgeAscMng, "AgeAscCluster" = AgeAscCluster,   # Remove MskAgeAscMng
                "IncAgeAscCluster" = IncAgeAscCluster, "PthAgeClusterCnt" = PthAgeClusterCnt, "IncPthAgeClusterCnt" = IncPthAgeClusterCnt)
  }
  return(NodLst)
}



AgeAscNod_Fcn = function(NodFrm, AllAscNod, Cluster){  
  AgeAscNod_Lst = list()
  AgeAscNod = NodFrm %>% filter(Cluster %in% Cluster) %>% pull("Nod") %>% as.character() #%>% rev()
  MskAgeAscNod = ifelse(AllAscNod %in% AgeAscNod,AllAscNod,"X")
  MskAgeAscNod[1] = "Top"
  AgeAscNod_Lst$AgeAscNod = AgeAscNod
  AgeAscNod_Lst$MskAgeAscNod = MskAgeAscNod
  return(AgeAscNod_Lst)
}


AgeAscMng_Fcn = function(NodFrm, AllAscMng, Cluster){  
  AgeAscMng_Lst = list()
  AgeAscMng = NodFrm %>% filter(Cluster %in% Cluster) %>% pull("Mng") %>% as.character()
  AgeAscMng = c("Top", AgeAscMng)
  MskAgeAscMng = ifelse(AllAscMng %in% AgeAscMng,AllAscMng,"X")
  MskAgeAscMng[1] = "Top"
  AgeAscMng_Lst$AgeAscMng = AgeAscMng
  AgeAscMng_Lst$MskAgeAscMng = MskAgeAscMng
return(AgeAscMng_Lst)
}



