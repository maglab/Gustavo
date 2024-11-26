### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgeHumGen, AgeModGen)
# Arc - Ageing-related community (e.g., GenArdArcFrm)
# Ard - Ageing-related disease (e.g., GenArdArcFrm)
# Arr - Array (e.g., NamArr, ColArr)
# Bin - Binary (e.g., NodOvlBinFrm)
# Col - Column (e.g., QryColArr)
# Cod - Code (e.g., KegCod)
# Dif - Difference (e.g., DifArr)
# Dis - Disease 
# Flt - Filtered (e.g., MrgFlt)
# Frm - Frame (e.g., GenArdArcFrm, AgeHumFrm)
# Gen - Gene (e.g., GenPhnFrm, GenAgeHum_Genes)
# Hum - Human (e.g., AgeHumGen, GenAge_Hum)
# Idx - Index (e.g., TopIdxArr)
# Len - Length (e.g., KegLen)
# Lst - List (e.g., OvlDisLst)
# Met - Method (e.g., NewUpsMet)
# Mod - Model_Organisms (e.g., AgeModGen, GenAge_Mod)
# Mrg - Merged (e.g., MrgTbl)
# Nam - Name (e.g., NewNamArr)
# Ovl - Overlap (e.g., NodOvlBinFrm)
# Phn - Phenotype (e.g., GenPhnFrm)
# Qry - Query (e.g., QryColArr)
# Ren - Rename (e.g., RenNamFcn)
# Sel - Select (e.g., dplyr::select())
# Set - Set (e.g., SetFrm)
# Sum - Summary (e.g., NodOvlBinFrm$Sum)
# Tbl - Table (e.g., MrgTbl)
# Top - Top (e.g., TopIdxArr)
# Typ - Type (e.g., DisTyp)
# Ups - Upset (e.g., NewUpsMet)

### LIBRARIES ##################################################################

library(UpSetR)
library(tidyr)
library(dplyr)

### START ######################################################################

GenArdArcFrm = readRDS('maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds') %>%
  dplyr::select(ARC_Meaning,Gene) %>%
  dplyr::rename(Phn=ARC_Meaning, Gen=Gene)

AgeHumGen = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeHum_Genes.rds")
AgeModGen = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeMod_Genes.rds")

AgeHumFrm = data.frame(Gen=AgeHumGen, Phn="GenAge_Hum")
AgeModFrm = data.frame(Gen=AgeModGen, Phn="GenAge_Mod")

#GenPhnFrm = GenArdArcFrm[,c('Gen','Phn')] %>% rbind(AgeHumFrm) %>% rbind(AgeModFrm) %>% unique()
GenPhnFrm = GenArdArcFrm %>% rbind(AgeHumFrm) %>% rbind(AgeModFrm) %>% unique()

# Convert the data frame to wider format
GenPhnFrm$value <- 1
wider_df <- GenPhnFrm %>% pivot_wider(names_from = Phn, values_from = value, values_fill = 0)
wider_dt = wider_df %>% as.data.frame()
row.names(wider_dt) = wider_dt$Gen
wider_dt$Gen = NULL


NamArr = GenPhnFrm$Phn %>% table() %>% sort() %>% rev() %>% names()
NewNamArr = c(NamArr[!(NamArr %in% c("GenAge_Hum", "GenAge_Mod"))], c("GenAge_Hum","GenAge_Mod")) %>% rev()

NewUpsMet = data.frame(sets = colnames(wider_dt), Met = colnames(wider_dt))

ColArr = c(rep( c("#023FA5","#4D60A9"),4),"#006027", "#41784E")# COOL!


names(ColArr) = c("cardiovascular", "endocrine/diabetes", "musculoskeletal/trauma",
                  "gastrointestinal/abdominal", "neurology/eye/psychiatry", 
                  "immunological/systemic disorders", "haematology/dermatology", "renal/urology",
                  "GenAge_Hum","GenAge_Mod")


### BINARY DATASET OF COMBINATIONS WITH SUM INFO ###############################

wider_dt$Merged <- apply(wider_dt, 1, function(row) paste(row, collapse = " "))
SetFrm = data.frame(Gen = row.names(wider_dt), Set = wider_dt$Merged)
wider_dt$Merged = NULL
YesAgeGen = c(AgeHumGen,AgeModGen)
NonAgeGen = row.names(wider_dt) %>% setdiff(YesAgeGen)

NodOvlBinFrm = wider_dt

QryColArr = colnames(wider_dt)
QryColArr = QryColArr[!(QryColArr %in% "Gen")]

NodOvlBinFrm$Sum = NodOvlBinFrm[QryColArr] %>% rowSums()
NodOvlBinFrm$Mrg = apply(NodOvlBinFrm, 1, function(row) paste(row, collapse = " "))
MrgTbl = NodOvlBinFrm$Mrg %>% table()
MrgFlt = MrgTbl[MrgTbl>1] %>% names() 
NodOvlBinFrmOrg = NodOvlBinFrm
NodOvlBinFrm = NodOvlBinFrm %>% dplyr::filter(Sum>0 & Mrg%in%MrgFlt)  
NodOvlBinFrm$Mrg = NULL

MaxSum = NodOvlBinFrm$Sum %>% max()
SumArr = NodOvlBinFrm$Sum %>% unique() %>% sort()
SumLen = length(SumArr)
SumTbl = NodOvlBinFrm$Sum %>% table()

AgeCol = QryColArr[grepl("Age",QryColArr)]
SumLen = length(AgeCol) + 1 
MaxSum = SumLen

### COLORS OF SETS AND ROWS ####################################################

# DIVERSITY THROGHOUT COLORS
ColLen = SumLen
ColArrDiv = c(hcl.colors(ColLen+1, palette = "Reds 2")[1:(ColLen-1)],  "#333333", "#104a8b")

# RETRIEVE LIST OF SET NAMES
OvlDisLst = NodOvlBinFrm %>% 
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
j=1 # INVENTADO POR MI
for(i in 1:OvlDisLen){
  print(i)
  DisArr = OvlDisLst[[i]]
  DisLen = DisArr %>% length()
  
  if("GenAge_Hum" %in% DisArr){
    DisTyp = 1
  } else if("GenAge_Mod" %in% DisArr){
    DisTyp = 2
  } else if("immunological/systemic disorders" %in% DisArr){
    DisTyp = 4 }else{
    DisTyp = 3
  }
  
  AllOvlLst[[i]] = list(query = intersects, params = list(DisArr), color = ColArrDiv[DisTyp], active = T)
  if(DisLen == MaxSum){
    MaxOvlLst[[j]] = AllOvlLst[[i]]
    j = j+1
  }
}

# CREATE ZEROS AND ONES FROM TRUE-FALSE
# MAX
MaxNodOvlBinFrm = NodOvlBinFrm %>% dplyr::filter(Sum == MaxSum)
MaxNodOvlBinFrm[QryColArr] = 1*MaxNodOvlBinFrm[QryColArr]

# ALL
NodOvlBinFrm[QryColArr] = 1*NodOvlBinFrm[QryColArr] 



### APPLY QUERY FOR SET-LINES COLORS ###########################################

PhnArr =  QryColArr
NodOvlBinFrm$Gen = row.names(NodOvlBinFrm)

UpsFrm = NodOvlBinFrm %>% dplyr::select(c("Gen",PhnArr,"Sum"))

DifArr = c("Gen", "Sum")

DifIdxArr = AllOvlLst %>% lapply(function(x){sum(unlist(x$params) %in% DifArr) == 0}) %>% unlist()
DifIdxLen = length(DifIdxArr)

j = 1
UpsQry = list()
for(i in 1:DifIdxLen){
  if(DifIdxArr[i]){
    UpsQry[[i]] = AllOvlLst[[i]]
  }
}



### PLOT UPSET #################################################################

NodOvlBinFrm$Sum = NULL

ColArr['GenAge_Hum'] = "#3C2692" 
ColArr['GenAge_Mod'] = "#6C64A2" 
ColArr[c(3,4,5,6,7,8)] = c("#023FA5","#4D60A9","#023FA5","#4D60A9","#023FA5","#4D60A9")

upset(NodOvlBinFrm,# wider_dt, 
      nsets = 10,
      nintersects = 90,
      keep.order = T, 
      number.angles = -45,
      sets = NewNamArr,
      line.size = 1.2, 
      text.scale = c(1.2,1,1,1,1.2,0.8),
      mainbar.y.label = "Ageing- and ARCs-related gene count", 
      sets.x.label = "Ageing- and ARCs-related gene count",
      queries = UpsQry,
      set.metadata = list(data = NewUpsMet, 
                          plots = list(list(type = "matrix_rows", column = "Met", colors = ColArr,
                                            alpha = 0.5)))
      
)


##################################################################################
### FUNCTIONS ####################################################################
##################################################################################

CmpSrtFcn = function(EndArr, PhnArr){
  DifNam = setdiff(PhnArr,names(EndArr)) 
  if(length(DifNam)>0){
    DifArr = rep(0,length(DifNam))
    names(DifArr) = DifNam
    FinArr = c(EndArr,DifArr) %>% sort() %>% rev()
  } else{
    FinArr = c(EndArr) %>% sort() %>% rev()
  }
  return(FinArr)
}


Arr2FrmFcn = function(CrdArr, PhnArr, Nam){
  tCrdFrm = CrdArr[PhnArr] %>% data.frame() 
  colnames(tCrdFrm) = Nam
  CrdFrm = t(tCrdFrm)
  return(CrdFrm)
}


normalize_by_max <- function(df) {
  normalized_df <- as.data.frame(t(apply(df, 1, function(row) row / max(row))))
  # Return the normalized data frame
  return(normalized_df)
}
