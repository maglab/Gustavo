# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Arr - Array (e.g., NtwArr)
# Dst - Distance (e.g., GenGenAll_DstFrm)
# Frm - Frame (e.g., GenGenAll_DstFrm, GenGenAll_PrxFrm)
# Grp - Graph (e.g., Grp)
# Len - Length (e.g., NtwLen)
# Ntw - Network (e.g., NtwArr)
# Prx - Proximity (e.g., GenGenAll_PrxFrm)
# Qry - Query (e.g., NtwQry)

# LIBRARIES ####################################################################

library(dplyr)
library(igraph)
library(rje)

# START ########################################################################

NtwArr = c("PPI","COX90","COX95","KEGG")
NtwLen = length(NtwArr)
for(i in 1:NtwLen){
  print(i)
  NtwQry = NtwArr[i]
  
  ### LOAD NETWROK DATA ########################################################
  
  Grp = paste("maglab/Gustavo/Data/Generated/",NtwQry,"/Lists/GraphsList.rds", sep='') %>% readRDS()
  
  ### CREATE AND SAVE GEN_GEN DISTANCE AND PROXIMITY FRAMES ####################
  
  GenGenAll_DstFrm = Grp$GenGenAll %>% distances() %>% as.matrix() %>% as.data.frame()
  saveRDS(GenGenAll_DstFrm, paste("maglab/Gustavo/Data/Generated/",NtwQry,"/Proximity_and_Distances/GenGen_DistanceFrame.rds", sep='') )
  
  GenGenAll_PrxFrm = GenGenAll_DstFrm %>% ProximityFunction()
  saveRDS(GenGenAll_PrxFrm, paste("maglab/Gustavo/Data/Generated/",NtwQry,"/Proximity_and_Distances/GenGen_ProximityFrame.rds", sep='') )
  
}

################################################################################
# FUNCTIONS
################################################################################

ProximityFunction = function(x){
  return(1/(1+x))
}

