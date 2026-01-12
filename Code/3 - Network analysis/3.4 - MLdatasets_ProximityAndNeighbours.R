### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Alg - Algorithm (i.e., metod for stablishing the relationship between gene and ARC/ARC: Average.Proximity2ARD, etc.)
# Age (e.g., AgeGenAll)
# Arc - Ageing-related Community (e.g., GenArdArcFrm)
# Ard - Ageing-related Community (e.g., GenArdArcFrm)
# Auc - Area Under the roc Curve (e.g., AucFrm)
# Cod - Coding (e.g., dplyr::rename(Cod=Code))
# Cox - Coexpression (e.g., "COX95")
# Dts - Dataset (e.g., DtsTypArr)
# Frm - Frame (e.g., AucFrm)
# Gen - Gene (e.g., AgeGenAll)
# Grp - Graph (e.g., Grp$GenGenNtw)
# Hum - Human (e.g., AgeGenAll)
# Hrc - Hierarchy (e.g., HrcArd)
# Len - Length
# Lst - List (e.g., AucGenLst)
# Ngb - Neighbours (e.g., AgeNgb)
# Nod - Node (e.g., Nod$GenNtw)
# Ntw - Network (e.g., NtwArr)
# Pth - Path (e.g., Pth = paste(...))
# Qry = Query (e.g., NtwQry)
# Uad - Union of Ageing(Human)- and Disease-related genes (e.g., UadNtwNgb)

# DATASET CREATION ALGOROTHMS
# Neighbours2ARC
# Neighbours2ARD
# Average.Proximity2ARC
# Closest.Proximity2ARC
# Closest.Proximity2ARD
# Average.Proximity2ARD

### LIBRARIES ##################################################################

library(dplyr)
library(igraph)
library(grid)
library(pROC)
library(data.table)



### INITIAL PARAMETERS #########################################################################################

GenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)

AgeGenAll = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
ArdGenAll = GenArdArcFrm$Gen %>% unique()

NtwArr = c("PPI", "COX95", "COX90", "KEGG")
NtwQry = NtwArr[1]
AucFrm = data.frame()
for(NtwQry in NtwArr){

  Pth = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/",sep="")
  
  FrmLst = paste(Pth,"Proximity_and_Distances/Proximities_To_ARC_ARD.rds", sep="") %>% readRDS()
  
  Nod = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/Lists/NodesList.rds", sep='') %>% readRDS()
  Grp = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/Lists/GraphsList.rds", sep='') %>% readRDS()
  
  AgeGen = AgeGenAll %>% intersect(Nod$GenNtw)
  ArdGen = ArdGenAll %>% intersect(Nod$GenNtw)
  
  AgeNgb = Grp$GenGenNtw %>% ego(order = 1, nodes = AgeGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  ArdNgb = Grp$GenGenNtw %>% ego(order = 1, nodes = ArdGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()

  AllNgb = c(ArdNgb, AgeNgb) %>% unique()
  ExcNgb = AllNgb %>% setdiff(AgeGen) %>% setdiff(ArdGen) # Exclusive neighbhours - Neighboours od Ageing and diseases that are not associated with them
  
  AucGenLst = list("AgeGen"= AgeGen, "ArdGen"= ArdGen,
                   "AgeNgb"= AgeNgb, "ArdNgb"= ArdNgb,
                   "AllNgb"= AllNgb, "ExcNgb"= ExcNgb)
  
  GenAll = Nod$GenAll
  GenNtw = Nod$GenNtw
  GenUadAll = c(ArdGenAll, AgeGenAll) %>% unique()
  GenUadNtw = intersect(GenUadAll, GenNtw)
  GenUadNtwNgb = ego(Grp$GenGenNtw, order = 1, nodes = GenUadNtw, mode = "all", mindist = 0) %>% unlist() %>% names() %>% unique()
  
  DtsLst = list()
  
  ### Average.Proximity2ARD #################################################
  
  Average.Proximity2ARD = FrmLst$Average.Proximity2ARD
  RowGenArr = row.names(Average.Proximity2ARD)
  DtsLst$Average.Proximity2ARD = Average.Proximity2ARD %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)

  sAucFrm = AucFcn(DtsLst$Average.Proximity2ARD, AucGenLst, NtwQry, Alg="Average.Proximity2ARD")
  AucFrm = rbind(AucFrm, sAucFrm)
  
  ### Closest.Proximity2ARD #################################################
  
  Closest.Proximity2ARD = FrmLst$Closest.Proximity2ARD
  RowGenArr = row.names(Closest.Proximity2ARD)
  DtsLst$Closest.Proximity2ARD = Closest.Proximity2ARD %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)

  sAucFrm = AucFcn(DtsLst$Closest.Proximity2ARD, AucGenLst, NtwQry, Alg="Closest.Proximity2ARD")
  AucFrm = rbind(AucFrm, sAucFrm)
  
  ### Closest.Proximity2ARC #############################################
  
  Closest.Proximity2ARC = FrmLst$Closest.Proximity2ARC
  RowGenArr = row.names(Closest.Proximity2ARC)
  DtsLst$Closest.Proximity2ARC = Closest.Proximity2ARC %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)

  sAucFrm = AucFcn(DtsLst$Closest.Proximity2ARC, AucGenLst, NtwQry, Alg="Closest.Proximity2ARC")
  AucFrm = rbind(AucFrm, sAucFrm)

  
  
  ### Average.Proximity2ARC #################################################
  
  Average.Proximity2ARC = FrmLst$Average.Proximity2ARC
  RowGenArr = row.names(Average.Proximity2ARC)
  DtsLst$Average.Proximity2ARC = Average.Proximity2ARC %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)

  sAucFrm = AucFcn(DtsLst$Average.Proximity2ARC, AucGenLst, NtwQry, Alg="Average.Proximity2ARC")
    AucFrm = rbind(AucFrm, sAucFrm)

  
  ### Neigbhours2ARD ###########################################################

  
  Neighbours2ARD = FrmLst$Neighbours2ARD
  RowGenArr = row.names(Neighbours2ARD)
  DtsLst$Neighbours2ARD = Neighbours2ARD %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)

  sAucFrm = AucFcn(DtsLst$Neighbours2ARD, AucGenLst, NtwQry, Alg="Neighbours2ARD")
  AucFrm = rbind(AucFrm, sAucFrm)

  
  ### Neighbours2ARC ####################################################
  
  Neighbours2ARC = FrmLst$Neighbours2ARC
  RowGenArr = row.names(Neighbours2ARC)
  DtsLst$Neighbours2ARC = Neighbours2ARC %>% 
    mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "Age", "NotAge")) %>%
    Keep_Rows(GenNtw)
  
  sAucFrm = AucFcn(DtsLst$Neighbours2ARC, AucGenLst, NtwQry, Alg="Neighbours2ARC")
  AucFrm = rbind(AucFrm, sAucFrm)

  
  ### SAVE #####################################################################
  
  DtsNamArr = DtsLst %>% names()
  DtsNamLen = length(DtsNamArr)

  
  for(i in 1:DtsNamLen){
    DtsNamQry = DtsNamArr[i]

    #paste(NtwQry,i,j) %>% print()
    paste(NtwQry,i) %>% print()

    QryFrm = DtsLst[[DtsNamQry]]
    RowNam = row.names(QryFrm)
    
    ClsArr = QryFrm$Class
    QryFrm$Class = NULL
    QryFrm$AgeGen = NULL
    QryFrm$ArdGen = NULL
    MeanFrm = rowMeans(QryFrm) %>% as.data.frame()
    MeanFrm$Gen = row.names(MeanFrm)
    MeanFrm$Class = ClsArr
    colnames(MeanFrm) = c("Scr","Gen","Class")
    MeanFrm = MeanFrm %>% select(Gen,Scr,Class)
    
    QryFrm = QryFrm %>% mutate(AgeGen = ifelse(RowNam %in% AgeGen, "Age", "NotAge"))
    
    OrgAucScr = AucFrm %>% filter(Grp==NtwQry, Gen=="Age", Alg==DtsNamQry) %>% unique() %>% pull(Auc)
    AucScr = (OrgAucScr*100) %>% round(0)
    
    
    QryFrm %>% write.csv( paste(Pth,"Ageing_Prediction/Datasets/",DtsNamQry,".csv",sep="") )
    #MeanFrm %>% write.csv( paste(Pth,"Ageing_Prediction/Predictions/Mean-Based/",DtsNamQry,"_Auc",AucScr,".csv",sep="") )

  }
} 


AucFrm = AucFrm %>% rename(Dataset=Alg,Gene=Gen,Network=Grp)

##############################################################################
# FUNCTIONS
##############################################################################


AucFcn = function(AllFrm, AucGenLst, NtwQry, Alg){
  
  AgeGen = AucGenLst$AgeGen
  ArdGen = AucGenLst$ArdGen
  ExcNgb = AucGenLst$ExcNgb

  AllFrm$Class = NULL
  
  AvgFrm = AllFrm %>% rowMeans() %>% as.data.frame()
  
  colnames(AvgFrm) = "Scr"
  
  AvgFrm$AgeGen = 0
  AvgFrm$ArdGen = 0
  
  AvgFrm[AgeGen,"AgeGen"] = 1
  AvgFrm[ArdGen,"ArdGen"] = 1
  AvgNamQry = "Scr"
  

  AucArr.AgeAll = auc(AvgFrm[,'AgeGen'], AvgFrm[,AvgNamQry]) %>% as.numeric()
  AucFrm.AgeAll = AucArr.AgeAll %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(Gen = "Age", Typ = "All", Grp = NtwQry)
  sAucFrm = AucFrm.AgeAll 
  colnames(sAucFrm) = c("Auc", "Gen", "Typ", "Grp")        
  
  sAucFrm$Alg = Alg
  sAucFrm = sAucFrm %>% select(Alg, Auc, Gen, Grp)
  
  return(sAucFrm)
  
}

Keep_Rows = function(df, rows){
  return(df[rows,])
}
