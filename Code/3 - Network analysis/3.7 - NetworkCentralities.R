# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgeModGen)
# Arc - Ageing-relaed Communities (e.g., GenArdArcFrm)
# Ard - Ageing-related Diseases (e.g., GenArdArcFrm)
# Arr - Array (e.g., NtwArr)
# Btw - Betweenness Centrality (e.g., BtwHumFrm)
# Cls - Closeness Centrality (e.g., ClsHumFrm)
# Cod - Coding (e.g., dplyr::rename(Cod=Code))
# Cnt - Centrality (e.g., CntFrm)
# Deg - Degree Centrality (e.g., colnames(IntFrm) = c("Gen", "Deg", "Ovl", "Ratio"))
# Eig - Eigenvalue Centrality (e.g., EigHumFrm)
# Frm - Frame (e.g., GenPhnFrm)
# Gen - Gene (e.g., GenArdArcFrm)
# Grp - Graph (e.g., QryGrp = Grp$GenGenNtw)
# Hum - Human (e.g., GenAgeHum_Genes)
# Imm - Immunological-Disorder Diseases (e.g., ImmGen)
# Int - Intersection (e.g., IntFrm)
# Mat - Matrix (e.g., MatCntFrm)
# Mod - Model Orgnisms (e.g., AgeModFrm)
# Ntw - Network (e.g., NtwArr)
# Oth - Others (e.g., OthGenGrp)
# Ovl - Overlap (e.g., colnames(IntFrm) = c("Gen", "Deg", "Ovl", "Ratio"))
# Phn - Phenotype (eihter ARCs or GenAgeHim/GenAgeMod) (e.g., GenPhnFrm)
# Pth - Path (e.g., Pth = paste(...))
# Qry - Quert (e.g., NtwQry)
# Rat - Ratio (e.g., Ratio)
# Trn - Transtivity (e.g., TrnHumFrm)

# LIBRARIES ####################################################################

library(igraph)

### LOADING DATA ###############################################################

NtwArr = c("PPI", "COX90", "COX95", "KEGG")

#GenPhnPhn = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnPhnFrm.rds")
GenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)


GenAge = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
GenPhnAge = data.frame("Gen"=GenAge,"Phn"="GenAge_Hum")
GenPhnFrm = data.frame(GenArdArcFrm$Gen, GenArdArcFrm$Arc) %>% unique() 
colnames(GenPhnFrm) = c("Gen", "Phn")
GenPhnFrm = GenPhnFrm%>% rbind(GenPhnAge)


AgeModGen = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")
AgeModFrm = data.frame(Gen=AgeModGen, Phn="GenAge_Mod")
GenPhnFrm = GenPhnFrm %>% rbind(AgeModFrm)

HumGen = GenPhnFrm %>% filter(Phn == "GenAge_Hum") %>% pull(Gen)
ModGen = GenPhnFrm %>% filter(Phn == "GenAge_Mod") %>% pull(Gen)
ImmGen = GenPhnFrm %>% filter(Phn == "immunological/systemic disorders") %>% pull(Gen)
ArcGen = GenPhnFrm %>% filter(!(Phn %in% c("GenAge_Hum", "GenAge_Mod", "cancer"))) %>% pull(Gen)
AgeGen = intersect(HumGen, ModGen)
ArcAgeGen = GenPhnFrm$Gen %>% unique()


## NEXT #################################################

NtwLen = length(NtwArr)
CntFrm = data.frame()
for(x in 1:NtwLen){
  
  NtwQry = NtwArr[x]
  Pth = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/",sep="")
  
  ### RETRIEVE DATA #################################################
  
  Grp = readRDS(paste(Pth,"Lists/","GraphsList.rds",sep=""))
  QryGrp = Grp$GenGenNtw
  IntFrm = calculate_overlaps(QryGrp, ArcGen)
  colnames(IntFrm) = c("Gen", "Deg", "Ovl", "Ratio")
  IntGrp = IntFrm$Ratio
  names(IntGrp) = IntFrm$Gen
  
  ### COMPUTE CENTALITIES ###########################################
  paste(x,".Deg",sep="") %>% print()
  DegGrp = igraph::degree(QryGrp) %>% sort() %>% rev()
  
  paste(x,".Btw",sep="") %>% print()
  BtwGrp = igraph::betweenness(QryGrp) %>% sort() %>% rev()
  
  paste(x,".Cls",sep="") %>% print()
  ClsGrp = igraph::closeness(QryGrp) %>% sort() %>% rev()
  
  paste(x,".Trn",sep="") %>% print()
  TrnGrp = full_transitivity(QryGrp) %>% sort() %>% rev()
  
  paste(x,".Eig",sep="") %>% print()
  EigGrp = igraph::eigen_centrality(QryGrp)$vector %>% sort() %>% rev() 


  ### MAP CENTRALITIES ##############################################
  
  GrpGen = V(QryGrp)$name
  HumGenGrp = intersect(HumGen, GrpGen)
  ModGenGrp = intersect(ModGen, GrpGen)
  ImmGenGrp = intersect(ImmGen, GrpGen)
  ArcGenGrp = intersect(ArcGen, GrpGen)
  AgeGenGrp = intersect(AgeGen, ModGen)
  OthGenGrp = setdiff(GrpGen, ArcAgeGen)

  # CREATE CENTRALITY FRAMES ---------------------------------------------------
  
  # DEGREE CENTRALITY
  DegHumFrm = data.frame(Phn="Hum", Scr=DegGrp[HumGenGrp], Cnt="Deg", Ntw=NtwQry)
  DegModFrm = data.frame(Phn="Mod", Scr=DegGrp[ModGenGrp], Cnt="Deg", Ntw=NtwQry)
  DegImmFrm = data.frame(Phn="Imm", Scr=DegGrp[ImmGenGrp], Cnt="Deg", Ntw=NtwQry)
  DegArcFrm = data.frame(Phn="Arc", Scr=DegGrp[ArcGenGrp], Cnt="Deg", Ntw=NtwQry)
  DegAgeFrm = data.frame(Phn="Age", Scr=DegGrp[AgeGenGrp], Cnt="Deg", Ntw=NtwQry)
  DegOthFrm = data.frame(Phn="Oth", Scr=DegGrp[OthGenGrp], Cnt="Deg", Ntw=NtwQry)
  
  # BETWEENESS CENTRALITY
  BtwHumFrm = data.frame(Phn="Hum", Scr=BtwGrp[HumGenGrp], Cnt="Btw", Ntw=NtwQry)
  BtwModFrm = data.frame(Phn="Mod", Scr=BtwGrp[ModGenGrp], Cnt="Btw", Ntw=NtwQry)
  BtwImmFrm = data.frame(Phn="Imm", Scr=BtwGrp[ImmGenGrp], Cnt="Btw", Ntw=NtwQry)
  BtwArcFrm = data.frame(Phn="Arc", Scr=BtwGrp[ArcGenGrp], Cnt="Btw", Ntw=NtwQry)
  BtwAgeFrm = data.frame(Phn="Age", Scr=BtwGrp[AgeGenGrp], Cnt="Btw", Ntw=NtwQry)
  BtwOthFrm = data.frame(Phn="Oth", Scr=BtwGrp[OthGenGrp], Cnt="Btw", Ntw=NtwQry)
  
  # CLOSENESS CENTRALITY
  ClsHumFrm = data.frame(Phn="Hum", Scr=ClsGrp[HumGenGrp], Cnt="Cls", Ntw=NtwQry)
  ClsModFrm = data.frame(Phn="Mod", Scr=ClsGrp[ModGenGrp], Cnt="Cls", Ntw=NtwQry)
  ClsImmFrm = data.frame(Phn="Imm", Scr=ClsGrp[ImmGenGrp], Cnt="Cls", Ntw=NtwQry)
  ClsArcFrm = data.frame(Phn="Arc", Scr=ClsGrp[ArcGenGrp], Cnt="Cls", Ntw=NtwQry)
  ClsAgeFrm = data.frame(Phn="Age", Scr=ClsGrp[AgeGenGrp], Cnt="Cls", Ntw=NtwQry)
  ClsOthFrm = data.frame(Phn="Oth", Scr=ClsGrp[OthGenGrp], Cnt="Cls", Ntw=NtwQry)
  
  # TRANSITIVITY (CLUSTERING COEFFICIENT)
  TrnHumFrm = data.frame(Phn="Hum", Scr=TrnGrp[HumGenGrp], Cnt="Trn", Ntw=NtwQry)
  TrnModFrm = data.frame(Phn="Mod", Scr=TrnGrp[ModGenGrp], Cnt="Trn", Ntw=NtwQry)
  TrnImmFrm = data.frame(Phn="Imm", Scr=TrnGrp[ImmGenGrp], Cnt="Trn", Ntw=NtwQry)
  TrnArcFrm = data.frame(Phn="Arc", Scr=TrnGrp[ArcGenGrp], Cnt="Trn", Ntw=NtwQry)
  TrnAgeFrm = data.frame(Phn="Age", Scr=TrnGrp[AgeGenGrp], Cnt="Trn", Ntw=NtwQry)
  TrnOthFrm = data.frame(Phn="Oth", Scr=TrnGrp[OthGenGrp], Cnt="Trn", Ntw=NtwQry) 
  
  # EIGENVALUE CENTRALITY
  EigHumFrm = data.frame(Phn="Hum", Scr=EigGrp[HumGenGrp], Cnt="Eig", Ntw=NtwQry)
  EigModFrm = data.frame(Phn="Mod", Scr=EigGrp[ModGenGrp], Cnt="Eig", Ntw=NtwQry)
  EigImmFrm = data.frame(Phn="Imm", Scr=EigGrp[ImmGenGrp], Cnt="Eig", Ntw=NtwQry)
  EigArcFrm = data.frame(Phn="Arc", Scr=EigGrp[ArcGenGrp], Cnt="Eig", Ntw=NtwQry)
  EigAgeFrm = data.frame(Phn="Age", Scr=EigGrp[AgeGenGrp], Cnt="Eig", Ntw=NtwQry)
  EigOthFrm = data.frame(Phn="Oth", Scr=EigGrp[OthGenGrp], Cnt="Eig", Ntw=NtwQry) 
  
  # INTERSECT
  IntHumFrm = data.frame(Phn="Hum", Scr=IntGrp[HumGenGrp], Cnt="Int", Ntw=NtwQry)
  IntModFrm = data.frame(Phn="Mod", Scr=IntGrp[ModGenGrp], Cnt="Int", Ntw=NtwQry)
  IntImmFrm = data.frame(Phn="Imm", Scr=IntGrp[ImmGenGrp], Cnt="Int", Ntw=NtwQry)
  IntArcFrm = data.frame(Phn="Arc", Scr=IntGrp[ArcGenGrp], Cnt="Int", Ntw=NtwQry)
  IntAgeFrm = data.frame(Phn="Age", Scr=IntGrp[AgeGenGrp], Cnt="Int", Ntw=NtwQry)
  IntOthFrm = data.frame(Phn="Oth", Scr=IntGrp[OthGenGrp], Cnt="Int", Ntw=NtwQry)
  
  # BINDING
  paste(x,".BndDeg",sep="") %>% print()
  DegFrm = DegHumFrm %>% rbind(DegModFrm) %>% rbind(DegImmFrm) %>% rbind(DegArcFrm) %>% rbind(DegAgeFrm) %>% rbind(DegOthFrm)
  
  paste(x,".BndBtw",sep="") %>% print()
  BtwFrm = BtwHumFrm %>% rbind(BtwModFrm) %>% rbind(BtwImmFrm) %>% rbind(BtwArcFrm) %>% rbind(BtwAgeFrm) %>% rbind(BtwOthFrm)
  
  paste(x,".BndCls",sep="") %>% print()
  ClsFrm = ClsHumFrm %>% rbind(ClsModFrm) %>% rbind(ClsImmFrm) %>% rbind(ClsArcFrm) %>% rbind(ClsAgeFrm) %>% rbind(ClsOthFrm)
  
  paste(x,".BndTrn",sep="") %>% print()
  TrnFrm = TrnHumFrm %>% rbind(TrnModFrm) %>% rbind(TrnImmFrm) %>% rbind(TrnArcFrm) %>% rbind(TrnAgeFrm) %>% rbind(TrnOthFrm)
  
  paste(x,".BndEig",sep="") %>% print()
  EigFrm = EigHumFrm %>% rbind(EigModFrm) %>% rbind(EigImmFrm) %>% rbind(EigArcFrm) %>% rbind(EigAgeFrm) %>% rbind(EigOthFrm)
  
  paste(x,".BndInt",sep="") %>% print()
  IntFrm = IntHumFrm %>% rbind(IntModFrm) %>% rbind(IntImmFrm) %>% rbind(IntArcFrm) %>% rbind(IntAgeFrm) %>% rbind(IntOthFrm)

  paste(x,".sBndAll",sep="") %>% print()
  sCntFrm = DegFrm %>% rbind(BtwFrm) %>% rbind(ClsFrm) %>% rbind(TrnFrm) %>% rbind(EigFrm) %>% rbind(IntFrm)
  
  paste(x,".BndAll",sep="") %>% print()
  CntFrm = rbind(CntFrm, sCntFrm)
  
}


### TABLE ##########################################################


MatCntFrm = transform_data(CntFrm)
MatCntFrm$Btw = round(MatCntFrm$Btw)
MatCntFrm$Deg = round(MatCntFrm$Deg)

MatCntFrm$Cls = paste(round(MatCntFrm$Cls*100), "%", sep="")
MatCntFrm$Eig = paste(round(MatCntFrm$Eig*100), "%", sep="")
MatCntFrm$Trn = paste(round(MatCntFrm$Trn*100), "%", sep="")
MatCntFrm$Int = paste(round(MatCntFrm$Int*100), "%", sep="")

MatCntFrm = MatCntFrm %>% select(Ntw,Phn,Deg,Btw,Cls,Eig,Trn,Int)
MatCntFrm = MatCntFrm %>% filter(Phn!="Age")
MatCntFrm$Phn[MatCntFrm$Phn=="Hum"] = "GenAge.Hum"
MatCntFrm$Phn[MatCntFrm$Phn=="Mod"] = "GenAge.Mod"
MatCntFrm$Phn[MatCntFrm$Phn=="Arc"] = "Diseases"
MatCntFrm$Phn[MatCntFrm$Phn=="Imm"] = "Immune Disorders"
MatCntFrm$Phn[MatCntFrm$Phn=="Oth"] = "Others"
MatCntFrm$Eig = NULL
colnames(MatCntFrm) = c("Network", "Genes", "Degree Centrality", "Betweenness Centrality", "Closeness Centrality",
                        "Clustering Coefficient", "Percentage of Neighbour Genes related with Diseases") 

NtwArr = MatCntFrm$Network %>% unique()
NtwNum = c(3,4,2,1) 
NtwFrm = data.frame(Network=NtwArr,NtwNum)

MatCntFrm = MatCntFrm %>% 
  merge(NtwFrm) %>%
  group_by(Network) %>%
  arrange(NtwNum, Genes) %>%
  ungroup() 
MatCntFrm$NtwNum = NULL

write.csv(MatCntFrm,"Data/Generated/Networks_and_predictions/Topological_data/Distance_and_centrality/Centrality_Table.csv")


################################################################################
### FUNCTIONS
################################################################################

calculate_overlaps <- function(graph, gene_set) {
  # Initialize a data frame to store the results
  result_df <- data.frame(Gene = V(graph)$name, Degree = igraph::degree(graph,mode="out"), Number_of_overlaps = NA)
  #result_df <- data.frame(Gene = V(graph)$name, Number_of_overlaps = NA)
  
  # Loop through each gene (node) in the graph
  NumNod = vcount(graph)
  i=3
  for (i in 1:NumNod) {
    
    paste(i,NumNod) %>% print()
    
    # Get the first-order neighbors of the gene
    neighbors <- neighbors(graph, i, mode="out")
    
    # Calculate the overlap with the gene set (excluding the gene itself)
    #overlap <- length(intersect(neighbors$name, gene_set)) - (V(graph)$name[i] %in% gene_set)
    overlap <- length(intersect(neighbors$name, gene_set)) 
    
    # Store the result in the data frame
    result_df$Number_of_overlaps[i] <- overlap
  }
  
  result_df$Ratio = result_df$Number_of_overlaps/result_df$Degree
  
  return(result_df)
}


full_transitivity <- function(graph) {
  # Check if the input is a graph object
  if (!inherits(graph, "igraph")) {
    stop("The input is not an igraph graph object.")
  }
  
  # Calculate local transitivity for non-isolated vertices
  transitivity_scores <- transitivity(graph, type = "local")
  
  # Get all vertex names
  all_vertex_names <- V(graph)$name
  
  # Initialize a named vector with NA values for all vertices
  full_scores <- setNames(rep(NA, length(all_vertex_names)), all_vertex_names)
  
  # Identify non-isolated vertices (vertices with at least one edge)
  non_isolated <- which(igraph::degree(graph) > 0)
  
  # Update the full_scores with the transitivity scores for non-isolated vertices
  full_scores[names(full_scores) %in% all_vertex_names[non_isolated]] <- transitivity_scores
  
  return(full_scores)
}
  

transform_data <- function(df) {
  # Reshape the data and compute the mean score for each combination of Ntw, Phn, and Cnt
  result_df <- df %>%
    group_by(Ntw, Phn, Cnt) %>%
    summarize(Mean_Scr = mean(Scr, na.rm = TRUE), .groups = 'drop') %>%
    spread(key = Cnt, value = Mean_Scr)
  
  return(result_df)
}
  