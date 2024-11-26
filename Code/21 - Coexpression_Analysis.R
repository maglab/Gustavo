# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., GenAgeHum)
# Arc - Ageing-related Communities (e.g., ArcArr)
# Ard - Ageing-related Diseases (e.g., GenAgeArd)
# Arr - Array (e.g., PhnArr)
# Bnf - Bonferroni_Correction (e.g., BnfTst)
# Bth - Both (e.g., BthScr)
# Col - Column (e.g., ColPhnSet)
# Cod - Code (e.g., Cod=Code)
# Cox - Coexpression (e.g., MatrixCox_AgeARD_Hgnc.rds)
# Dir - Direct (e.g., DirIntFrm)
# Frm - Frame (e.g., GenArcFrm)
# Gen - Gene (e.g., GenAgeHum)
# Hum - Human (e.g., GenAgeHum)
# Lng - Long (e.g., LngPhnGenFrm)
# Lst - List (e.g., LngPhnPhnLst)
# Mat - Matrix (e.g., GenMat)
# Mix - Mixed (Col and Row)
# Mod - Model_Organisms (e.g., GenAgeMod)
# Ntw - Network (e.g., NtwGen)
# Oth - Others (e.g., OthGen)
# Phn - Phenotype (e.g., PhnGen)
# Scr - Score (e,g., PhnPhnScrFrm)
# Set - Set (e.g., PhnSetArr)
# Tst - Test (pvalue) (e.g., TstTxt)
# Txt - Text (e.g., TstTxt)
# Typ - Type (e.g., TypInt)
# Val - Value (e.g., LngPhnGenFrm$Val)

### LIBRARIES ##################################################################

library(dplyr)
library(tidyr)
library(rstatix)
library(reshape2)

## START #######################################################################

# LOADING GEN-ARC FRAME 
GenArcFrm = readRDS("maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds") %>%
  dplyr::rename(Phn=ARC_Meaning, Gen=Gene) %>%
  select(Phn, Gen)

# LOADING GenAge
GenAgeHum = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeMod = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeMod_Genes.rds")

GenAgeFrm = data.frame(Phn="GenAge", Gen=GenAgeHum) %>% 
  rbind( data.frame(Phn="ModAge", Gen=GenAgeMod) )

GenMat =  readRDS("maglab/Gustavo/Data/Generated/Networks_Creation/Coexpression/MatrixCox_AgeARD_Hgnc.rds")

# LAST GENE-COMMUNITY DATA
NtwGen = colnames(GenMat)
PhnGen = GenArcFrm$Gen %>% unique()
OthGen = setdiff(NtwGen, PhnGen)
PhnArr = GenArcFrm$Phn %>% unique()
PhnLen = length(PhnArr)


################################################################################
### PLEIOTROPIC ASSOCIATIONS ###################################################
################################################################################

### GENERATE DATAFRAMES WITH SUM DATA ##########################################

# GENERATE THREE COLUMN DATAFRAME
LngPhnGenFrm = GenArcFrm %>% select(Phn, Gen) %>% unique()
LngPhnGenFrm$Val = 1

# CONVERT THREE COLUMNS FRAME TO RECTANGULAR
MatPhnGenFrm <- LngPhnGenFrm %>%
  spread(Phn, Val)

# Set row names using the 'row' column and remove it from the data frame
rownames(MatPhnGenFrm) <- MatPhnGenFrm$Gen
MatPhnGenFrm <- MatPhnGenFrm %>% select(-Gen)
MatPhnGenFrm[is.na(MatPhnGenFrm)] = 0

ArcArr = MatPhnGenFrm %>% colnames() 


# CREATE COLUMN WITH MERGED COMMUNITIES SEPPARATED BY " - "
PhnSetArr = c()
for(i in 1:nrow(MatPhnGenFrm)){
  print(i)
  PhnSetArr[i] = colnames(MatPhnGenFrm[i,])[MatPhnGenFrm[i,] == 1] %>% paste(collapse=" - ")
}


# INCORPORATE SUM AND MERGED COMMUNITIES
MatPhnGenFrm$PhnNum = rowSums(MatPhnGenFrm)
MatPhnGenFrm$OrgPhnSet = PhnSetArr
MatPhnGenFrm$Gen = row.names(MatPhnGenFrm)


# MAKE NAME OF COMMUNITIES SMALLER
GenSumPhn = MatPhnGenFrm %>% select(Gen, PhnNum, OrgPhnSet)
OrgPhnSetArr = GenSumPhn$OrgPhnSet

AbvPhnSetArr = full_shorten_sub_elements(GenSumPhn$OrgPhnSet)

PhnPhnFrm = data.frame(OrgPhnSet=OrgPhnSetArr, AbvPhnSet=AbvPhnSetArr)
PhnPhnFrm = unique(PhnPhnFrm)

GenSumPhn$AbbPhnSet = AbvPhnSetArr

GenSumPhn = GenSumPhn %>% arrange(desc(PhnNum))
GenSumPhn = GenSumPhn %>% filter(PhnNum>0)

### PLEIOTROPY DATAFRAMES ######################################################

NumPhnPhnFrm = SumPhnFcn(GenSumPhn, GenMat)  # Get coexpression by Group genes
NumPhnPhnFrm$ColPhnNum = full_count_dashes(NumPhnPhnFrm$ColPhnSet)
NumPhnPhnFrm$RowPhnNum = full_count_dashes(NumPhnPhnFrm$RowPhnSet)

IntPhnPhnFrm = NumPhnPhnFrm %>% filter(ColPhnSet == RowPhnSet) %>% select(Scr, ColPhnSet, ColPhnNum, Dbs) %>% rename(Phn=ColPhnSet, Num=ColPhnNum)

################################################################################
### COMPUTE DATA ###############################################################
################################################################################

DirIntFrm = GenSumPhn %>% mutate(Phn=paste("ARC-Pleiotropy_",PhnNum,sep="")) %>% 
  select(Phn, Gen)


GenPhnFrm = rbind(GenArcFrm,DirIntFrm,GenAgeFrm)

NewPhnArr = GenPhnFrm$Phn %>% unique()
NewPhnLen = length(NewPhnArr)

PhnPhnFrm = data.frame()

for(x in 1:NewPhnLen){
  print(x)
  QryPhnCol = NewPhnArr[x]
  
  ColPhnSetGen = GenPhnFrm %>% filter(Phn %in% QryPhnCol) %>% pull(Gen) %>% unique()
  
  ### CROSSED COMMUNITIES ######################################################
  sPhnPhnFrm = data.frame()
  for(y in 1:NewPhnLen){
    
    paste(x, y) %>% print()
    QryPhnRow = NewPhnArr[y]
    
    RowPhnSetGen = GenPhnFrm %>% filter(Phn %in% QryPhnRow) %>% pull(Gen) %>% unique()
    
    RowPhnSetGenNtw = intersect(RowPhnSetGen, NtwGen)
    ColPhnSetGenNtw = intersect(ColPhnSetGen, NtwGen)
    
    Row.NonCol.Gen = setdiff(RowPhnSetGenNtw,ColPhnSetGenNtw)
    Col.NonRow.Gen = setdiff(ColPhnSetGenNtw,RowPhnSetGenNtw)
    Bth.Gen = intersect(RowPhnSetGenNtw,ColPhnSetGenNtw)
    
    BthGenMat = GenMat[Bth.Gen,Bth.Gen] 
    
    BthScr <- BthGenMat[upper.tri(BthGenMat, diag = TRUE)]
    Row.NonCol.Bth.Scr = GenMat[Bth.Gen,Row.NonCol.Gen] %>% unlist() %>% as.numeric()
    Col.NonRow.Bth.Scr = GenMat[Bth.Gen,Col.NonRow.Gen] %>% unlist() %>% as.numeric()
    DifScr = GenMat[Col.NonRow.Gen,Row.NonCol.Gen] %>% unlist() %>% as.numeric()
    
    PhnPhn.Arr = c(BthScr,Row.NonCol.Bth.Scr,Col.NonRow.Bth.Scr,DifScr)
    PhnPhn.Arr = PhnPhn.Arr[PhnPhn.Arr!=1] # Remove elements of the main diagonal (a gene coexpresing with itself equals coexpression 1)
    
    ssPhnPhnFrm = data.frame(Scr = PhnPhn.Arr, ColPhnSet = QryPhnCol, RowPhnSet = QryPhnRow)
    
    sPhnPhnFrm = rbind(sPhnPhnFrm, ssPhnPhnFrm)
    
  }
  
  ### PUTTING ALL TOGETHER #####################################################
  
  PhnPhnFrm = rbind(PhnPhnFrm, sPhnPhnFrm)
  
}

PhnPhnScrFrm = PhnPhnFrm %>% rename(RowPhn=RowPhnSet,ColPhn=ColPhnSet)

################################################################################
### STATISTICS #################################################################
################################################################################

ColQryPhn = 'haematology/dermatology'

PhnPhnTstFrm = data.frame()
for(i in 1:NewPhnLen){
  
  ### PREPARE FRAME ############################################################
  
  #print(i)
  ColQryPhn = NewPhnArr[[i]]
  
  ColQryPhnTxt = ColQryPhn
  
  if(ColQryPhn == "immunological/systemic disorders"){
    ColQryPhnTxt = "immunological/systemic\n disorders"
  }
  

  ### COMPUTE COM-COM DYNAMICS #################################################

  sPhnPhnTstFrm = data.frame()
  for(j in 1:NewPhnLen){
    
    paste(i,j) %>% print()
    
    RowQryPhn = NewPhnArr[[j]]
    
    # COL.COL
    ColScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% ColQryPhn, RowPhn %in% ColQryPhn) %>% pull(Scr)
    ColMeanPhnPhn = ColScrPhnPhn %>% mean()
    ColLenPhnPhn = ColScrPhnPhn %>% length()
    
    # ROW.ROW
    RowScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% RowQryPhn, RowPhn %in% RowQryPhn) %>% pull(Scr)
    RowMeanPhnPhn = RowScrPhnPhn %>% mean()
    RowLenPhnPhn = RowScrPhnPhn %>% length()
    
    # MIX - COL.ROW
    MixScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% ColQryPhn, RowPhn %in% RowQryPhn) %>% pull(Scr)
    MixMeanPhnPhn = MixScrPhnPhn %>% mean()
    MixLenPhnPhn = MixScrPhnPhn %>% length()

    
    # STATISTICAL TESTS
    ColRow.Tst = tryCatch({t.test(ColScrPhnPhn, RowScrPhnPhn)$p.value}, error = function(e){NA})

    # CREATE FRAME
    ssPhnPhnTstFrm =
      data.frame(ColPhn=ColQryPhn,RowPhn=RowQryPhn,
                 ColMean=ColMeanPhnPhn,RowMean=RowMeanPhnPhn, MixMean=MixMeanPhnPhn,
                 ColRow.Dif = ColMeanPhnPhn-RowMeanPhnPhn, 
                 ColLen = ColLenPhnPhn, RowLen=RowLenPhnPhn, MixLen=MixLenPhnPhn, 
                 ColRow.Tst
      )
    
    sPhnPhnTstFrm = rbind(sPhnPhnTstFrm, ssPhnPhnTstFrm) 
    
  }
  
  PhnPhnTstFrm = rbind(sPhnPhnTstFrm, PhnPhnTstFrm)
  
}

### TESTING ADJUSTMENTS #######################################################

# BONFERRININ CORRECTION
PhnPhnTstFrm$ColRow.TstBnf = p.adjust(PhnPhnTstFrm$ColRow.Tst, method="bonferroni")

# ORIGINAL TEST LABEL
PhnPhnTstFrm$ColRow.TstTxt = TstTxtFcn(PhnPhnTstFrm$ColRow.Tst)


# CORRECTED TEST LABEL
PhnPhnTstFrm$ColRow.TstBnfTxt = TstTxtFcn(PhnPhnTstFrm$ColRow.TstBnf)

### MEAN ARRAYS ################################################################

PhnMeanFrm = PhnPhnTstFrm %>% select(ColPhn,ColMean) %>% unique() %>% rename(Phn=ColPhn,Mean=ColMean)
row.names(PhnMeanFrm) = PhnMeanFrm$Phn
PhnMeanFrm$Phn=NULL
PhnMeanFrm$Mean = round(PhnMeanFrm$Mean,2)

################################################################################
### PLOTTING FRAMES ############################################################
################################################################################

LngPhnPhnLst = list()
MatPhnPhnLst = list()

### THRE COLUMNS CASE ##########################################################

# ORIGINAL TEST
LngPhnPhnLst$Tst = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.Tst) %>% rename(Scr=ColRow.Tst)

# CORRECTD TEST
LngPhnPhnLst$TstBnf = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstBnf) %>% rename(Scr=ColRow.TstBnf)

# ORIGINAL TEST LABEL
LngPhnPhnLst$TstTxt = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstTxt) %>% rename(Scr=ColRow.TstTxt)

# CORRECTD TEST LABEL
LngPhnPhnLst$TstBnfTxt = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstBnfTxt) %>% rename(Scr=ColRow.TstBnfTxt)

# MEAN DIFFERENCE
LngPhnPhnLst$Dif = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.Dif) %>% rename(Scr=ColRow.Dif)

# MIX
LngPhnPhnLst$Mix = PhnPhnTstFrm %>% select(ColPhn,RowPhn,MixMean) %>% rename(Scr=MixMean)

### RECTANGULAR CASE ###########################################################

# ORIGINAL TEST 
MatPhnPhnLst$Tst = spread(LngPhnPhnLst$Tst, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# CORRECTED TEST 
MatPhnPhnLst$TstBnf = spread(LngPhnPhnLst$TstBnf, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# ORIGINAL TEST LABEL
MatPhnPhnLst$TstTxt = spread(LngPhnPhnLst$TstTxt, key=ColPhn, value=Scr) %>% rename(Phenotype=RowPhn)

# CORRECTED TEST LABEL
MatPhnPhnLst$TstBnfTxt = spread(LngPhnPhnLst$TstBnfTxt, key=ColPhn, value=Scr) %>% rename(Phenotype=RowPhn)

# MEAN DIFFERENCE
MatPhnPhnLst$Dif = spread(LngPhnPhnLst$Dif, key=ColPhn, value=Scr) %>% rename(Phenotype=RowPhn)

# MIX
MatPhnPhnLst$Mix = spread(LngPhnPhnLst$Mix, key=ColPhn, value=Scr) %>% rename(Phenotype=RowPhn)


################################################################################
### SAVING #####################################################################
################################################################################

PhnPhnScrFrm = PhnPhnScrFrm %>% rename(Column_Phenotype=ColPhn,Row_Phenotype=RowPhn)

saveRDS(MatPhnPhnLst, "maglab/Gustavo/Data/Generated/Coexpression/Coexpression_Matrices_List.rds")  
saveRDS(PhnMeanFrm, "maglab/Gustavo/Data/Generated/Coexpression/Self_Coexpression_Mean_Frame.rds") 
saveRDS(PhnPhnScrFrm, "maglab/Gustavo/Data/Generated/Coexpression/Coexpression_Pairwise_Frame.rds")   


################################################################################
# FUNCTIONS
################################################################################

full_shorten_sub_elements = function(input_vector){
  output_vector <- sapply(input_vector, shorten_sub_elements)
  return(as.character(output_vector))
}

# Function to shorten each sub-element
shorten_sub_elements <- function(x) {
  sub_elements <- unlist(strsplit(x, " - "))
  shortened <- substr(sub_elements, 1, 4)
  paste(shortened, collapse = "-",sep="")
}


SumPhnFcn = function(GenSumPhn, GenMat){ 
  
  RowNtwGen = row.names(GenMat)
  
  SetArr = GenSumPhn$OrgPhnSet %>% unique()
  SetLen = length(SetArr)
   
  PhnArrFrm = data.frame()
  PhnPhnFrm = data.frame()
  for(x in 1:SetLen){
    ColQry = SetArr[x]
    
    ColPhnSetGen = GenSumPhn %>% filter(OrgPhnSet %in% ColQry) %>% pull(Gen) %>% unique()
    sPhnPhnFrm = data.frame()
    
    for(y in 1:SetLen){
      
      paste(x, y) %>% print()
      RowQry = SetArr[y]
      
      RowPhnSetGen = GenSumPhn %>% filter(OrgPhnSet %in% RowQry) %>% pull(Gen) %>% unique()
      
      RowPhnSetGenNtw = intersect(RowPhnSetGen, RowNtwGen)
      ColPhnSetGenNtw = intersect(ColPhnSetGen, RowNtwGen)
      
      UppPhnPhn = GenMat[ColPhnSetGenNtw,RowPhnSetGenNtw]
      LowPhnPhn = GenMat[RowPhnSetGenNtw,ColPhnSetGenNtw]
      
      UppPhnPhn.Arr = UppPhnPhn %>% unlist() %>% as.numeric()
      LowPhnPhn.Arr = LowPhnPhn %>% unlist() %>% as.numeric()
      
      PhnPhn.Arr = c(UppPhnPhn.Arr, LowPhnPhn.Arr)
      
      PhnPhn.Len = length(PhnPhn.Arr)
      if(PhnPhn.Len==0){
        PhnPhn.Arr = NA
      }
      
      ssPhnPhnFrm = data.frame(Scr = PhnPhn.Arr, ColPhnSet = ColQry, RowPhnSet = RowQry, Dbs=TypInt)
      
      sPhnPhnFrm = rbind(sPhnPhnFrm, ssPhnPhnFrm)
      
    }
   
    PhnPhnFrm = rbind(PhnPhnFrm, sPhnPhnFrm)

  }
  
  return(PhnPhnFrm)
}

# Apply the function to each element of the input vector
full_count_dashes <- function(input_vector) {
  dash_counts <- sapply(input_vector, count_dashes)
  return(as.numeric(dash_counts)+1)
}

# Function to count the number of "-" symbols in a string
count_dashes <- function(x) {
  sum(strsplit(x, "")[[1]] == "-")
}

# Pvalues to astherisks
TstTxtFcn = function(Tst){
  TstTxt = rep("ns",length(Tst))
  TstTxt = ifelse(Tst  <= 5e-2, "*", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-2, "**", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-3, "***", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-4, "****", TstTxt)
  TstTxt[is.na(TstTxt)] = "ns"
  return(TstTxt)
}

reorder_ids <- function(df) {
  df <- df %>%
    rowwise() %>%
    mutate(
      FirstID = ifelse(RowID <= ColumnID, RowID, ColumnID),
      SecondID = ifelse(RowID > ColumnID, RowID, ColumnID)
    ) %>%
    select(-RowID, -ColumnID)
  return(df)
}
