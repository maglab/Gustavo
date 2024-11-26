### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# ABBREVIATIONS
# Arc - Arcmunity (e.g., ArcArr, ArcQry)
# Ard - Phenotype (e.g., GenArdArcFrm, ArdArr)
# Arr - Array (e.g., NtwArr)
# Avg - Average (e.g., sAvgPrx2Ard)
# Cls - Closest (e.g., ClsDst2Ard)
# Dst - Distance (e.g., GenGen_DstFrm, QryGenGen_DstFrm)
# Frm - Frame (e.g., GenGen_DstFrm, GenGen_PrxFrm)
# Gen - Gene (e.g., ArdQryGen)
# Len - Length (e.g., NtwLen, ArdLen, ArcLen)
# Max - Maximum (e.g., Compute_column(type="Max"))
# Min - Minimum (e.g., Compute_column(type="Min"))
# Ngb - Neighbors (e.g., sGenArdNgb_PrxFrm)
# Prx - Proximity (e.g., GenGen_PrxFrm, QryGenGen_PrxFrm)
# Qry - Query (e.g., ArdQry, ArcQry)

## GENE RELATIONSHIPS with ARC and ARD
# ClsDst2Ard - Closest.Distance2ARD
# ClsPrx2Ard - Closest.Proximity2ARD
# ClsPrx2Arc - Closest.Proximity2ARC
# AvgPrx2Ard - Average.Proximity2ARD
# AvgPrx2Arc - Average.Proximity2ARC
# Ngb2Ard - Neighbours2ARD
# Ngb2Arc - Neighbours2ARC

## LIBRARIES ###################################################################

library(dplyr)
library(rje)

################################################################################

NtwArr = c("PPI","COX90","COX95","KEGG")
NtwLen = length(NtwArr)

for(x in 2:NtwLen){
  print(x)
  QryNtw = NtwArr[x]

  GenGen_DstFrm = paste("maglab/Gustavo/Data/Generated/",QryNtw,"/Proximity_and_Distances/","GenGen_DistanceFrame.rds",sep="") %>% readRDS()
  GenGen_PrxFrm = paste("maglab/Gustavo/Data/Generated/",QryNtw,"/Proximity_and_Distances/","GenGen_ProximityFrame.rds",sep="") %>% readRDS()
  
  
  
  GenArdArcFrm = readRDS("maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds") %>% 
    dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)
  
  # ARDs #########################################################################
  ArdArr = GenArdArcFrm$Ard %>% unique()
  ArdLen = length(ArdArr)
  for(i in 1:ArdLen){
    paste(x,1,i) %>% print()
    ArdQry = ArdArr[i]
    
    # RETRIEVE GENES ASSOCIATED WITH THE QUERIES ARD
    ArdQryGen = GenArdArcFrm %>% filter(Ard %in% ArdQry) %>% pull(Gen)
    
    # RETRIEVE DISTANCES/PROXIMITY FOR THE QUERIED ARD-RELATED GENES
    QryGenGen_DstFrm = GenGen_DstFrm[,ArdQryGen,drop=FALSE]
    QryGenGen_PrxFrm = GenGen_PrxFrm[,ArdQryGen,drop=FALSE]
    
    #COMPUTE PROXIMITY AND DISTANCES TO THE QUERIED ARD
    sClsDst2Ard = QryGenGen_DstFrm %>% Compute_column(type="Min",colname=ArdQry)
    sClsPrx2Ard = QryGenGen_PrxFrm %>% Compute_column(type="Max",colname=ArdQry)
    sAvgPrx2Ard = QryGenGen_PrxFrm %>% Compute_column(type="Mean",colname=ArdQry) 
    sNgb2Ard = QryGenGen_PrxFrm %>% Compute_column(type="Ngb",colname=ArdQry) 
    if(i==1){
      ClsDst2Ard = sClsDst2Ard
      ClsPrx2Ard = sClsPrx2Ard
      AvgPrx2Ard = sAvgPrx2Ard
      Ngb2Ard = sNgb2Ard
    } else{
      ClsDst2Ard = cbind(ClsDst2Ard,sClsDst2Ard)
      ClsPrx2Ard = cbind(ClsPrx2Ard,sClsPrx2Ard)
      AvgPrx2Ard = cbind(AvgPrx2Ard,sClsPrx2Ard)
      Ngb2Ard = cbind(Ngb2Ard,sNgb2Ard)
    }
    
  }
  

  # ARCs #########################################################################
  
  ArcArr = GenArdArcFrm$Arc %>% unique()
  ArcLen = length(ArcArr)
  for(i in 1:ArcLen){
    paste(x,2,i) %>% print()
    ArcQry = ArcArr[i]
    
    # RETRIEVE GENES ASSOCIATED WITH THE QUERIES ARC
    ArcQryGen = GenArdArcFrm %>% filter(Arc %in% ArcQry) %>% pull(Gen) %>% unique()
    
    # RETRIEVE DISTANCES/PROXIMITY FOR THE QUERIED ARC-RELATED GENES
    QryGenGen_DstFrm = GenGen_DstFrm[,ArcQryGen,drop=FALSE]
    QryGenGen_PrxFrm = GenGen_PrxFrm[,ArcQryGen,drop=FALSE]
    
    #COMPUTE PROXIMITY AND DISTANCES TO THE QUERIED ARC
    sClsDst2Arc = QryGenGen_DstFrm %>% Compute_column(type="Min",colname=ArcQry)
    sClsPrx2Arc = QryGenGen_PrxFrm %>% Compute_column(type="Max",colname=ArcQry)
    sAvgPrx2Arc = QryGenGen_PrxFrm %>% Compute_column(type="Mean",colname=ArcQry) 
    sNgb2Arc = QryGenGen_PrxFrm %>% Compute_column(type="Ngb",colname=ArcQry) 
    if(i==1){
      ClsDst2Arc = sClsDst2Arc
      ClsPrx2Arc = sClsPrx2Arc
      AvgPrx2Arc = sAvgPrx2Arc
      Ngb2Arc = sNgb2Arc
    } else{
      ClsDst2Arc = cbind(ClsDst2Arc,sClsDst2Arc)
      ClsPrx2Arc = cbind(ClsPrx2Arc,sClsPrx2Arc)
      AvgPrx2Arc = cbind(AvgPrx2Arc,sAvgPrx2Arc)
      Ngb2Arc = cbind(Ngb2Arc,sNgb2Arc)
    }
    
  }
  

  # SAVE #######################################################################
  
  GenArd_PrxDstLst = list(Closest.Distance2ARD=ClsDst2Ard, Closest.Proximity2ARD=ClsPrx2Ard, 
                          Average.Proximity2ARD=AvgPrx2Ard, Neighbours2ARD=Ngb2Ard,
                          Closest.Distance2ARC=ClsDst2Arc, Closest.Proximity2ARC=ClsPrx2Arc, 
                          Average.Proximity2ARC=AvgPrx2Arc, Neighbours2ARC=Ngb2Arc)
  
  saveRDS(GenArd_PrxDstLst, paste("maglab/Gustavo/Data/Generated/",QryNtw,"/Proximity_and_Distances/","Proximities_To_ARC_ARD.rds", sep="") )

}



################################################################################
# FUNCTIONS
################################################################################

ProximityFunction = function(x){
  return(1/(1+x))
}

Compute_column <- function(data, type="Max", colname = "Col") {
  # Arcpute row sums
  
  if(type=="Max"){
    row_values <- rowMaxs(data)
  }
  if(type=="Min"){
    row_values <- rowMins(data)
  }
  if(type=="Mean"){
    row_values <- rowMeans(data)
  }
  if(type=="Ngb"){
    row_values <- count_row_elements_equal_to_0.5(data)
  }
  
  
  # Create a data frame with the row sums as a single row
  result <- data.frame(row_values)
  
  # Set the column names to the original row names of the input data
  row.names(result) <- rownames(data)
  
  # Set the row name for the result
  colnames(result) <- colname
  
  return(result)
}


count_row_elements_equal_to_0.5 <- function(data) {
  # Apply a function to each row that counts the number of elements equal to 0.5
  row_counts <- apply(data, 1, function(row) sum(row == 0.5))
  return(row_counts)
}
