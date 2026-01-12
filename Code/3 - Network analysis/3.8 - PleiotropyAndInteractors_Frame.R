# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., GenAgeHum_Genes)
# Arc - Ageing-related Communities (e.g., DirGenArdArcFrm)
# Ard - Ageing-related Diseases (e.g., DirGenArdArcFrm)
# Arr - Array (e.g., NtwArr)
# Cod - Code (e.g., dplyr::rename(Cod=Code))
# Dir - Direct (e.g., DirGenArdArcFrm)
# Frm - Frame (e.g., DirIntArdFrm)
# Gen - Gene (e.g., DirGenArdArcFrm)
# Hum - Human (e.g., GenAgeHum_Genes)
# Ind - Indirect (e.g., IndGenArdArcFrm)
# Int - Interactor (e.g., DirIntArdFrm)
# Len - Length (e.g., NtwLen)
# Ntw - Network(e.g., NtwArr)
# Pth - Path (e.g., Pth = paste(...))
# Qry - Query (e.g., QryNtw)

# LIBRARIES ####################################################################

library(dplyr)
library(igraph)
library(rje)

# START ########################################################################

# LOADING FATA 
DirGenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)

GenAgeHum_Genes = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")

# NUMBER OF DIRECT INTERACTIONS WITH ARDs
DirIntArdFrm = DirGenArdArcFrm %>% 
  select(Gen,Ard) %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(Gen='.',DirIntArd='Freq')

# NUMBER OF DIRECT INTERACTIONS WITH ARCs
DirIntArcFrm = DirGenArdArcFrm %>% 
  select(Gen,Arc) %>% 
  unique() %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame %>% 
  rename(Gen='.',DirIntArc='Freq')

# LOOP OVER THE NETWORKS
NtwArr = c('PPI','COX90','COX95','KEGG')
NtwLen = length(NtwArr)
i=1
for(i in 1:NtwLen){

  print(i)
  QryNtw = NtwArr[i]
  
  Pth = paste("Data/Generated/Networks_and_predictions/Networks/",QryNtw,"/",sep="")
  
  IndGenArdArcFrm = paste(Pth,"Proximity_and_Distances/","Indirect_Interactions_Long.csv",sep="") %>% read.csv()
  # NUMBER OF INDIRECT INTERACTIONS WITH ARDs
  IndIntArdFrm = IndGenArdArcFrm %>% 
    select(Gen,Ard) %>% 
    pull(Gen) %>% 
    table() %>% 
    as.data.frame() %>% 
    rename(Gen='.',IndIntArd='Freq')
  
  # NUMBER OF DIRECT INTERACTIONS WITH ARCs
  IndIntArcFrm = IndGenArdArcFrm %>% 
    select(Gen,Arc) %>% 
    unique() %>% 
    pull(Gen) %>% 
    table() %>% 
    as.data.frame() %>% 
    rename(Gen='.',IndIntArc='Freq')

  # MERGING DIRECT AND INDIRECT INTRACTIONS
  DirIndIntFrm = DirIntArdFrm %>% 
    merge(DirIntArcFrm,by='Gen',all=TRUE) %>%
    merge(IndIntArdFrm,by='Gen',all=TRUE) %>%
    merge(IndIntArcFrm,by='Gen',all=TRUE) %>% 
    mutate(GenAge = ifelse(Gen %in% GenAgeHum_Genes, "Age", "NotAge"))
  
  DirIndIntFrm[is.na(DirIndIntFrm)] = 0
  
  # RENAMING FOR SAVING
  DirIndIntFrm = DirIndIntFrm %>% 
    rename(ARD_Pleiotropy=DirIntArd,
           ARC_Pleiotropy=DirIntArc,
           iARD_Interactor=IndIntArd,
           iARC_Interactor=IndIntArc)
  
  DirIndIntFrm %>% write.csv(paste(Pth,"Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv",sep=""))
  
}
