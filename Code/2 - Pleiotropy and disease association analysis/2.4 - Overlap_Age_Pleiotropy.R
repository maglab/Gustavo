library(dplyr)
library(UpSetR)
library(grid)
library(data.table)
library(tidyr) 

### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Arc - Ageing-related community (e.g., GenArdArcFrm)
# Ard - Ageing-related disease (e.g., GenArdArcFrm)
# Arr - Array (e.g., GenAgeHumArr)
# Ass - Association (e.g., ArcAssFrm)
# Col - Column (e.g., colnames(UpsFrm))
# Cnt - Count (e.g., ArcArr)
# Dis - Disease (e.g., DisGenArr)
# Dst - Distance (e.g., HumArcDstColArr)
# Fcn - Function (e.g., UpsPltFcn)
# Frm - Frame (e.g., GenArdArcFrm)
# Frq - Frequency
# Gen - Gene (e.g., GenArdArcFrm)
# Hum - Human (e.g., GenAge.Hum)
# Hrc - Hierarchy (e.g., ArdHrcFrm)
# Idx - Index (e.g., DifIdxArr)
# Inf - Infinite (e.g., GenArcInf)
# Len - Length (e.g., DisGenLen)
# Log - Logical (e.g., UpsNam = UpsNam[!(LogRow)])
# Mod - Model organism (e.g., GenAge.Mod)
# Nam - Name (e.g., UpsNam)
# Non - Non-Human (e.g., NonHumArc)
# Num - Number (e.g., ArcNumFrm)
# Org - Original (e.g., OrgUpsFrm)
# Ovl - Overlap (e.g., AllOvlLst)
# Plt - Pleiotropy (e.g., ArcPltFrm)
# Qry - Query (e.g., GenQry)
# Row - Row (e.g., RowGen)
# Tit - Title (e.g., TitTxt)
# Txt - Text (e.g., TitTxt)
# Typ - Type (e.g., DisTyp)
# Ups - Upset plot (e.g., UpsFrm)
# Yes - Yes (e.g., YesHumArc)

### START ######################################################################

# LOADING DATA

GenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)

GenAgeHumArr = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
GenAgeModArr = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")

ArcPltFrm = GenArdArcFrm %>% 
  dplyr::select(Arc,Gen) %>% 
  unique() %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
colnames(ArcPltFrm) = c('Gen','Frq')


# Outer pleitropy
GenPlt6 = ArcPltFrm %>% filter(Frq==6) %>% pull(Gen) %>% as.character()
GenPlt5 = ArcPltFrm %>% filter(Frq==5) %>% pull(Gen) %>% as.character()
GenPlt4 = ArcPltFrm %>% filter(Frq==4) %>% pull(Gen) %>% as.character()
GenPlt3 = ArcPltFrm %>% filter(Frq==3) %>% pull(Gen) %>% as.character()
GenPlt2 = ArcPltFrm %>% filter(Frq==2) %>% pull(Gen) %>% as.character()
GenPlt1 = ArcPltFrm %>% filter(Frq==1) %>% pull(Gen) %>% as.character()

RowGen = c(GenAgeHumArr, GenAgeModArr, GenPlt1, GenPlt2, GenPlt3, GenPlt4, GenPlt5, GenPlt6) %>% unique()


UpsFrm = 
  data.frame(
             'GenAge.Hum' = RowGen %in% GenAgeHumArr,
             'GenAge.Mod' = RowGen %in% GenAgeModArr,
             'ARC-Pleitropy_1' = RowGen %in% GenPlt1,
             'ARC-Pleitropy_2' = RowGen %in% GenPlt2,
             'ARC-Pleitropy_3' = RowGen %in% GenPlt3,
             'ARC-Pleitropy_4' = RowGen %in% GenPlt4,
             'ARC-Pleitropy_5' = RowGen %in% GenPlt5,
             'ARC-Pleitropy_6' = RowGen %in% GenPlt6
  ) * 1

row.names(UpsFrm) = RowGen

colnames(UpsFrm) = c("GenAge.Hum","GenAge.Mod","ARC-Pleitropy_1","ARC-Pleitropy_2","ARC-Pleitropy_3","ARC-Pleitropy_4","ARC-Pleitropy_5","ARC-Pleitropy_6")

# PRINT UPSET
ap = UpsPltFcn(UpsFrm)
ap


################################################################################
# FUNCTIONS
################################################################################

UpsPltFcn = function(OrgUpsFrm){

  OrgUpsFrm[,'Gen'] = row.names(OrgUpsFrm)
  UpsNam = colnames(OrgUpsFrm)
  UpsNam = UpsNam[!(UpsNam %in% "Gen")]
  
  
  QryColArr = colnames(OrgUpsFrm)
  QryColArr = QryColArr[!(QryColArr %in% "Gen")]
  QryColArr
  
  OrgUpsFrm$Sum = OrgUpsFrm[QryColArr] %>% rowSums()
  OrgUpsFrm = OrgUpsFrm %>% filter(Sum>0)
  
  MaxSum = OrgUpsFrm$Sum %>% max()
  SumArr = OrgUpsFrm$Sum %>% unique() %>% sort()
  SumLen = length(SumArr)
  SumTbl = OrgUpsFrm$Sum %>% table()
  
  
  PltCol = QryColArr[grepl("Pleitropy",QryColArr)]
  AgeCol = QryColArr[grepl("Age",QryColArr)]
  SumLen = length(AgeCol) + 1 # AQUI LE MODIFIQUE!!!!!!!!!!!!!!!
  MaxSum = SumLen
  
  # RANGE BETWEEN TWO COLORS
  ColLen = SumLen
  FunColRng <- colorRampPalette(c("#48AAAD", "grey"))   # Apply colorRampPalette
  ColArrTwo <- FunColRng(ColLen) 
  #pie(rep(1, ColLen), col = ColArrTwo, main = "Pie plot with rainbow")  # TEST COLOR
  
  # DIVERSITY THROGHOUT MANY COLORS
  ColLen = SumLen
  ColArrDiv = c(hcl.colors(ColLen+1, palette = "Reds 2")[1:(ColLen-1)],  "#333333")
  #pie(rep(1, ColLen), col = ColArrDiv, main = "Pie plot with rainbow")  # TEST COLOR

  # RETRIEVE LIST OF NAMES
  OvlDisLst = OrgUpsFrm %>% 
    unite(col='ColPst', QryColArr, sep='-') %>% 
    pull(ColPst) %>% 
    unique() %>% 
    strsplit("-") %>% 
    lapply(function(x){x==TRUE | x==1}) %>%
    lapply(function(x){QryColArr[x]})
  
  NumUseSet = OvlDisLst %>% unlist() %>% unique() %>% length()
  OvlLstLenArr = OvlDisLst %>% lapply(function(x){length(x)}) %>% unlist()
  OvlLstLenMax = max(OvlLstLenArr)
  TopIdxArr = which(OvlLstLenArr %in% OvlLstLenMax)
  TopIdxLen = length(TopIdxArr)
  
  OvlDisNamArr = c()
  for(i in 1:TopIdxLen){
    Idx = TopIdxArr[i]
    OvlDisNamArr = c(OvlDisNamArr, OvlDisLst[[Idx]])
  }
  OvlDisNamArr = OvlDisNamArr %>% unique()
  TopDisOvlLen = length(OvlDisNamArr)
  
  # GENERATE QUERIES FOR UPSET, BARCOLORS
  AllOvlLst = list()
  MaxOvlLst = list()
  OvlDisLen = length(OvlDisLst)
  i = 1
  for(i in 1:OvlDisLen){
    print(i)
    DisArr = OvlDisLst[[i]]
    DisLen = DisArr %>% length()
    
    if("GenAge.Hum" %in% DisArr){
      DisTyp = 1
    } else if("GenAge.Mod" %in% DisArr){
      DisTyp = 2
    } else{
      DisTyp = 3
    }

    j=1
    AllOvlLst[[i]] = list(query = intersects, params = list(DisArr), color = ColArrDiv[DisTyp], active = T)
    if(DisLen == MaxSum){
      MaxOvlLst[[j]] = AllOvlLst[[i]]
      j = j+1
    }
  }
  
  # CREATE ZEROS AND ONES FROM TRUE-FALSE
  MaxOrgUpsFrm = OrgUpsFrm %>% filter(Sum == MaxSum)
  MaxOrgUpsFrm[QryColArr] = 1*MaxOrgUpsFrm[QryColArr]
  
  # ALL
  OrgUpsFrm[QryColArr] = 1*OrgUpsFrm[QryColArr] 
  
  
  # META FRAME  
  MetaFrm = data.frame(sets = QryColArr, Met = QryColArr)
  
  ### BY CLUSTER ###############################################################
  
  PltLen = PltCol %>% length()
  AgeLen = AgeCol %>% length()
  QryColArr
  
  
  Col1 <- hcl.colors(PltLen+2, palette = "Blues 2")[1:PltLen] %>% rev()
  Col2 <- hcl.colors(AgeLen+3, palette = "Purples 2")[1:(AgeLen)] %>% rev()
  

  
  PhnArr =  QryColArr[QryColArr!="Sum"]
  A = PhnArr[1]
  B = PhnArr[2]
  PhnArr[1] = B
  PhnArr[2] = A
  
  
  ColArr = c(Col2, Col1)
  names(ColArr) = PhnArr

  
  UpsFrm = OrgUpsFrm %>% dplyr::select(c("Gen",PhnArr,"Sum"))
  UpsMeta = MetaFrm %>% filter(Met %in% PhnArr)
  
  DifArr = c("Gen", "Sum")
  
  DifIdxArr = AllOvlLst %>% lapply(function(x){sum(unlist(x$params) %in% DifArr) == 0}) %>% unlist()
  DifIdxLen = length(DifIdxArr)

  
  UpsQry = list()
  for(i in 1:DifIdxLen){
    if(DifIdxArr[i]){
      UpsQry[[i]] = AllOvlLst[[i]]
    }
  }
  

  # PRINT UPSET
  TitTxt = "ARCs' pleitropies and GenAge's PPI-based neighbourhood"
  upset(UpsFrm ,
        keep.order = T,
        number.angles = 45,
        sets = PhnArr, 
        point.size = 3,
        line.size = 1.5, 
        text.scale = c(1.2,1,1.2,1,1.2,1.1),
        nsets = NumUseSet,
        mainbar.y.label = "Gene Count", 
        sets.x.label = "Gene Count", 
        queries = UpsQry,
        set.metadata = list(data = UpsMeta, 
                            plots = list(list(type = "matrix_rows", column = "Met", colors = ColArr,#colors = c(Amigo = "green", Amiga = "red"), 
                                              alpha = 0.5)))
  )

  
}

