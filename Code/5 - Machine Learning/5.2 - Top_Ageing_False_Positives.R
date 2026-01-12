### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# ARC -  (e.g., ARC_Pleiotropy)
# COX - Coexpression (e.g., IndFrm_COX90)
# Dir - Direct (e.g., KEGG.DirIndFrm)
# Dts - Dataset (e.g., Dts)
# Ens - Ensembl (e.g., EnsHgnFrm)
# Frm - Frame (e.g., TauFrm)
# Gen - Gene (e.g., Gen)
# Hgn - HGNC (e.g., EnsHgnFrm)
# Ind - Indirect (e.g., IndFrm_PPI)
# Int - Interactor (e.g., TypInt)
# Num - Number (e.g., NumThr)
# Shr - Share (e.g., ShrGen)
# Tau- Tau (e.g., Tau)
# Thr - Threshold (e.g., NumThr)
# Typ - Type (e.g., TypInt)

### LIBRARIES ##################################################################

library(ggpubr)
library(dplyr)

# KEGG PLEIOTROPY AND INTERACTORS
KEGG.DirIndFrm = read.csv('Data/Generated/Networks_and_predictions/Networks/KEGG/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') 
KEGG.DirIndFrm = KEGG.DirIndFrm %>% select(Gen, iARC_Interactor) %>% rename(ARC_Interactor_KEGG=iARC_Interactor)

PPI.DirIndFrm = read.csv('Data/Generated/Networks_and_predictions/Networks/PPI/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') 
PPI.DirIndFrm = PPI.DirIndFrm %>% select(Gen, ARC_Pleiotropy, iARC_Interactor) %>% rename(ARC_Interactor_PPI=iARC_Interactor)

KEGG.DirIndFrm$X=NULL
PPI.DirIndFrm$X=NULL


# PPI PREDICTIONS AND PLEIOTROPY/INTERACTIONS
#OrgMUX = read.table("Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions/MachineLearning-Based/All_Int_Mux_Auc88.csv", header = TRUE, sep=",")
OrgMUX = read.table("Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions/Predictions/Ard_Auc88.csv", header = TRUE, sep=",")
OrgMUX$X = NULL
OrgMUX = OrgMUX %>% rename(Gen=Label)
MUX = merge(PPI.DirIndFrm,OrgMUX,by="Gen") %>% merge(KEGG.DirIndFrm,by="Gen") %>% 
  filter(Class == "NotAge") %>%
  select(Gen, Prob, ARC_Pleiotropy, ARC_Interactor_PPI, ARC_Interactor_KEGG) %>%
  arrange(desc(Prob))
  #select(Gen, Prob, Class, Test, ARC_Pleiotropy, iARC_Interactor, ARD_Pleiotropy, iARD_Interactor)
MUX$Rank = 1:nrow(MUX)
MUX <- MUX %>% select(Rank, everything())

write.csv(MUX,'Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions/Top_Genes_Best_Algorithm/TopGenes_All_Int_Mux_Auc88.csv', row.names=FALSE)
