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
# Tau - Tau
# Tst - Test (pvalue) (e.g., TstTxt)
# Txt - Text (e.g., TstTxt)
# Typ - Type (e.g., TypInt)
# Val - Value (e.g., LngPhnGenFrm$Val)

### LIBRARIES ##################################################################

library(dplyr)
library(tidyr)
library(rstatix)
library(reshape2)

### START ######################################################################

# LOADING GENE-ARC FRAME
GenArcFrm = readRDS("maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds") %>%
  dplyr::rename(Phn=ARC_Meaning, Gen=Gene) %>%
  select(Phn, Gen)

# LOADING GenAge
GenAgeHum = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeMod = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeMod_Genes.rds")

# MAERING GenAge and Gene-ARC
GenAgeFrm = data.frame(Phn="GenAge", Gen=GenAgeHum) %>% 
  rbind( data.frame(Phn="ModAge", Gen=GenAgeMod) )

# Retreieve gene specificity frame
OrgTauFrm = read.csv("maglab/Gustavo/Data/Retrieved/Specificity/Tau_gene_V8.csv")
EnsHgnFrm = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/HgncEnsembl.rds")

row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
OrgTauFrm$Gen = EnsHgnFrm[OrgTauFrm$gene_id,"external_gene_name"]
TauFrm = OrgTauFrm %>% select(Gen, tau) %>% filter(!is.na(tau))
rm(OrgTauFrm)

# LAST GENE-COMMUNITY DATA
NtwGen = TauFrm$Gen
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

################################################################################
### COMPUTE DATA ###############################################################
################################################################################

DirIntFrm = GenSumPhn %>% mutate(Phn=paste("ARC-Pleiotropy_",PhnNum,sep="")) %>% 
  select(Phn, Gen)

GenPhnFrm = rbind(GenArcFrm,DirIntFrm,GenAgeFrm)

GenPhnTauFrm = GenPhnFrm %>% merge(TauFrm,by='Gen') %>% filter(!is.na(tau))

#dim(PhnTauFrm )
PhnTauFrm = GenPhnTauFrm %>% select(Gen, Phn, tau) %>% unique() %>% rename(Scr=tau)
PhnTauFrm$color_group = 'Other_ARCs'
PhnTauFrm$Lbl = 'Tissue specificity in genes associated to ARCs and GenAge'
PhnTauFrm$color_group[PhnTauFrm$Phn %in% "immunological/systemic disorders"] = "immunological/systemic disorders"
PhnTauFrm$color_group[PhnTauFrm$Phn %in% unique(GenAgeFrm$Phn)] = "immunological/systemic disorders"
PhnTauFrm$color_group[PhnTauFrm$Phn %in% unique(DirIntFrm$Phn)] = "Pleiotropies"

NewPhnArr = PhnTauFrm$Phn %>% unique()
NewPhnLen = length(NewPhnArr)

PhnPhnFrm = data.frame(Scr  =  PhnTauFrm$Scr,
                       ColPhn = PhnTauFrm$Phn,
                       RowPhn = PhnTauFrm$Phn
)

PhnPhnScrFrm = PhnPhnFrm


################################################################################
### STATISTICS #################################################################
################################################################################

PhnPhnTstFrm = data.frame()
for(i in 1:NewPhnLen){
  
  ### PREPARE FRAME ############################################################

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
    
    
    # EXT.EXT
    ColScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% ColQryPhn, RowPhn %in% ColQryPhn) %>% pull(Scr)
    ColMeanPhnPhn = ColScrPhnPhn %>% mean()
    ColLenPhnPhn = ColScrPhnPhn %>% length()
    
    # INN.INN
    RowScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% RowQryPhn, RowPhn %in% RowQryPhn) %>% pull(Scr)
    RowMeanPhnPhn = RowScrPhnPhn %>% mean()
    RowLenPhnPhn = RowScrPhnPhn %>% length()
    
    # MIX - EXT.INN
    MixScrPhnPhn = PhnPhnScrFrm %>% filter(ColPhn %in% ColQryPhn, RowPhn %in% RowQryPhn) %>% pull(Scr)
    MixMeanPhnPhn = MixScrPhnPhn %>% mean()
    MixLenPhnPhn = MixScrPhnPhn %>% length()
    
    
    # STATISTICAL TESTS
   ColRow.Tst = tryCatch({wilcox.test(ColScrPhnPhn, RowScrPhnPhn)$p.value}, error = function(e){NA})

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

PhnPhnTstFrm %>% filter(ColPhn=="ModAge",RowPhn=="ModAge")

dim(PhnPhnTstFrm)

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
LngPhnPhnLst$Tst = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.Tst) %>% rename(Scr=ColRow.Tst) %>% unique()

# CORRECTD TEST
LngPhnPhnLst$TstBnf = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstBnf) %>% rename(Scr=ColRow.TstBnf) %>% unique()

# ORIGINAL TEST LABEL
LngPhnPhnLst$TstTxt = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstTxt) %>% rename(Scr=ColRow.TstTxt) %>% unique()

# CORRECTD TEST LABEL
LngPhnPhnLst$TstBnfTxt = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.TstBnfTxt) %>% rename(Scr=ColRow.TstBnfTxt) %>% unique()

# MEAN DIFFERENCE
LngPhnPhnLst$Dif = PhnPhnTstFrm %>% select(ColPhn,RowPhn,ColRow.Dif) %>% rename(Scr=ColRow.Dif) %>% unique()

# MIX
LngPhnPhnLst$Mix = PhnPhnTstFrm %>% select(ColPhn,RowPhn,MixMean) %>% rename(Scr=MixMean) %>% unique()

### RECTANGULAR CASE ###########################################################

LngPhnPhnLst$Tst %>%
  count(RowPhn, ColPhn) %>%
  filter(n > 1)

# ORIGINAL TEST 
MatPhnPhnLst$Tst = spread(LngPhnPhnLst$Tst, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# CORRECTED TEST 
MatPhnPhnLst$TstBnf = spread(LngPhnPhnLst$TstBnf, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# ORIGINAL TEST LABEL
MatPhnPhnLst$TstTxt = spread(LngPhnPhnLst$TstTxt, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# CORRECTED TEST LABEL
MatPhnPhnLst$TstBnfTxt = spread(LngPhnPhnLst$TstBnfTxt, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# MEAN DIFFERENCE
MatPhnPhnLst$Dif = spread(LngPhnPhnLst$Dif, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)

# MIX
MatPhnPhnLst$Mix = spread(LngPhnPhnLst$Mix, key=ColPhn, value=Scr) %>% dplyr::rename(Phenotype=RowPhn)


################################################################################
### SAVING #####################################################################
################################################################################

PhnPhnScrFrm = PhnPhnScrFrm %>% rename(Column_Phenotype=ColPhn,Row_Phenotype=RowPhn)

saveRDS(MatPhnPhnLst, "maglab/Gustavo/Data/Generated/Specificity/Specificity_Matrices_List.rds") 
saveRDS(PhnMeanFrm, "maglab/Gustavo/Data/Generated/Specificity/Specificity_Mean_Frame.rds") 
saveRDS(PhnPhnScrFrm, "maglab/Gustavo/Data/Generated/Specificity/Specificity_Pairwise_Frame.rds") 


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
