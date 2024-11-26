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

library(scatterpie)
library(ggpubr)
library(cowplot)
library(ggpubr)

### LOADING DATA ###############################################################

# KEGG PLEIOTROPY AND INTERACTORS
KEGG.DirIndFrm = read.csv('maglab/Gustavo/Data/Generated/KEGG/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') 
PPI.DirIndFrm = read.csv('maglab/Gustavo/Data/Generated/PPI/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') 
KEGG.DirIndFrm$X=NULL
PPI.DirIndFrm$X=NULL

# PPI PREDICTIONS AND PLEIOTROPY/INTERACTIONS
OrgPPI = read.table("maglab/Gustavo/Data/Generated/PPI/Ageing_Prediction/Predictions/MachineLearning-Based/Closest.Proximity2ARD_Auc81.csv", header = TRUE, sep=",")
OrgPPI$X = NULL
OrgPPI = OrgPPI %>% rename(Gen=Label)
PPI = merge(PPI.DirIndFrm,OrgPPI) %>% 
  select(Gen, Prob, Class, Test, ARC_Pleiotropy, iARC_Interactor, ARD_Pleiotropy, iARD_Interactor)


# KEGG PREDICTIONS AND PLEIOTROPY/INTERACTIONS
OrgKEGG = read.table("maglab/Gustavo/Data/Generated/KEGG/Ageing_Prediction/Predictions/MachineLearning-Based/Closest.Proximity2ARD_Auc83.csv", header = TRUE, sep=",")
OrgKEGG$X = NULL
OrgKEGG = OrgKEGG %>% rename(Gen=Label)
KEGG = merge(KEGG.DirIndFrm,OrgKEGG) %>% 
  select(Gen, Prob, Class, Test, ARC_Pleiotropy, iARC_Interactor, ARD_Pleiotropy, iARD_Interactor)


# SAVING AGEING-PROB-DESCENDING FRAMES
# SAVING PPI
PPI %>% arrange(desc(Prob)) %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredPPI_All.csv")
PPI %>% arrange(desc(Prob)) %>% filter(Class=="Age") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredPPI_Age.csv")
PPI %>% arrange(desc(Prob)) %>% filter(Class=="NotAge") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredPPI_NotAge.csv")
# SAVING KEGG
KEGG %>% arrange(desc(Prob)) %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredKEGG_All.csv")
KEGG %>% arrange(desc(Prob)) %>% filter(Class=="Age") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredKEGG_Age.csv")
KEGG %>% arrange(desc(Prob)) %>% filter(Class=="NotAge") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PredKEGG_NotAge.csv")


# INTEGRATING PPI AND KEGG PREDICTIONS

KEGG = KEGG %>% select(Gen, Prob, Class, Test, ARC_Pleiotropy, iARC_Interactor, ARD_Pleiotropy, iARD_Interactor) %>%
  rename(Prob.KEGG = Prob, iARC_Interactor.KEGG=iARC_Interactor, iARD_Interactor.KEGG = iARD_Interactor)

PPI = PPI %>% select(Gen, Prob, iARC_Interactor, iARD_Interactor) %>%
  rename(Prob.PPI = Prob, iARC_Interactor.PPI=iARC_Interactor, iARD_Interactor.PPI = iARD_Interactor)

# SGARED GENES BETWEEN PPI AND KEGG
ShrGen = intersect(KEGG$Gen, PPI$Gen)

# RANKING COMMONG GENES IN PPI AND KEBHB BASED ON AGEING-PROB
ShrKEGG = KEGG %>% filter(Gen %in% ShrGen)
ShrKEGG$Rank.KEGG = 1:nrow(ShrKEGG)

ShrPPI = PPI %>% filter(Gen %in% ShrGen)
ShrPPI$Rank.PPI = 1:nrow(ShrPPI)

# KEGG-PPI INTEGRATION
KEGG_PPI = merge(ShrKEGG, ShrPPI, by="Gen")
KEGG_PPI = KEGG_PPI %>% mutate(Rank.Mean = (Rank.KEGG+Rank.PPI)/2)

KEGG_PPI = KEGG_PPI %>% mutate(Prob.Mean = (Prob.KEGG+Prob.PPI)/2)
KEGG_PPI = KEGG_PPI %>% mutate(iARC_Interactor.Mean = (iARC_Interactor.KEGG+iARC_Interactor.PPI)/2 )
KEGG_PPI = KEGG_PPI %>% mutate(iARD_Interactor.Mean = (iARD_Interactor.KEGG+iARD_Interactor.PPI)/2 )
KEGG_PPI = KEGG_PPI %>% arrange(desc(Prob.Mean))


KEGG_PPI=
  KEGG_PPI %>% select(Gen, Class, Prob.Mean, Prob.KEGG, Prob.PPI, ARC_Pleiotropy, ARD_Pleiotropy, 
                  iARC_Interactor.KEGG, iARC_Interactor.PPI, iARC_Interactor.Mean, 
                  iARD_Interactor.KEGG, iARD_Interactor.PPI, iARD_Interactor.Mean, 
                  Rank.KEGG, Rank.PPI, Rank.Mean, 
                  Test)

# SAVING KEGG_PPI
KEGG_PPI %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PairwisePred_All.csv")
KEGG_PPI %>% filter(Class=="Age") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PairwisePred_Age.csv")
KEGG_PPI %>% filter(Class=="NotAge") %>% write.csv("maglab/Gustavo/Data/Generated/Pairwise_Predictions/PairwisePred_NotAge.csv")

# CORRELATION
cor(KEGG_PPI$iARC_Interactor.KEGG, KEGG_PPI$iARC_Interactor.PPI)
cor(KEGG_PPI$Pred.x, KEGG_PPI$Pred.y)

