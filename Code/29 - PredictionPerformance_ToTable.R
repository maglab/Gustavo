### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Alg - Algorithm (ML or Mean) (e.g., AlgFrm_Avg)
# Arc - Ageing-related Communities (e.g., Proximity2ARC)
# Arr - Array (e.g., ColArr_Avg)
# Auc -Area Under the Roc Curve (e.g., Auc=round(Auc,2))
# Avg - Average (e.g., Typ="Avg")
# Col - Column (e.g., Col=paste(Network,"-",Typ,sep=""))
# Cox - Coexpression (e.g., COX90-Avg)
# Frm - Frame (e.g., AlgFrm_Avg)
# Mat - Matrix (e.g., MatAlgFrm)
# ML - Machine-Learning (e.g., Typ="ML")
# PPI - Protein-Protein Interactions
# Typ - Type (e.g., Typ="Avg")

### LIBRARIES ##################################################################

library(dplyr)
library(tidyr)

# START ########################################################################

# LOAD MEAN-BASED PERFORMANCES
AlgFrm_Avg = read.csv("maglab/Gustavo/Data/Generated/Ageing_Prediction_Performances/Mean_Based_AUC.csv") %>%
  select(Dataset, Auc, Network) %>% 
  mutate(Typ="Avg") %>%
  mutate(Col=paste(Network,"-",Typ,sep="")) %>%
  select(Dataset,Auc,Col) %>%
  mutate(Auc=round(Auc,2))
  
# LOAD ML-BASED PERFORMANCES
AlgFrm_ML = read.csv("maglab/Gustavo/Data/Generated/Ageing_Prediction_Performances/ML_Based_AUC.csv") %>%
  select(Dataset, Auc, Network) %>% 
  mutate(Typ="ML") %>%
  mutate(Col=paste(Network,"-",Typ,sep="")) %>%
  select(Dataset,Auc,Col) %>%
  mutate(Auc=Auc/100)  %>%
  mutate(Auc=round(Auc,2))

# - INTEGRATION OF ML AND MEAN ALGORITHMS --------------------------------------

ColArr_Avg = AlgFrm_Avg$Col %>% unique()
ColArr_ML = AlgFrm_ML$Col %>% unique()

AlgFrm = rbind(AlgFrm_Avg,AlgFrm_ML)

MatAlgFrm <- AlgFrm %>%
  pivot_wider(names_from = Col, values_from = Auc) %>% 
  as.data.frame()
row.names(MatAlgFrm) = MatAlgFrm$Dataset
MatAlgFrm$Dataset = NULL
MatAlgFrm = MatAlgFrm[c('Closest.Proximity2ARD','Closest.Proximity2ARC','Average.Proximity2ARD','Average.Proximity2ARC','Neighbours2ARD','Neighbours2ARC'),
                      c('PPI-Avg','PPI-ML','COX90-Avg','COX90-ML','COX95-Avg','COX95-ML','KEGG-Avg','KEGG-ML')] 
MatAlgFrm[['Mean-Avg']] = MatAlgFrm %>% select(ColArr_Avg) %>% rowMeans() %>% round(2)
MatAlgFrm[['Mean-ML']] = MatAlgFrm %>% select(ColArr_ML) %>% rowMeans() %>% round(2)
MatAlgFrm[['Mean-All']] = MatAlgFrm %>% select(ColArr_Avg,ColArr_ML) %>% rowMeans() %>% round(2)
MatAlgFrm['Mean',] = MatAlgFrm %>% colMeans() %>% round(2)

write.csv(MatAlgFrm,'maglab/Gustavo/Data/Generated/Ageing_Prediction_Performances/All_Based_AUC.csv')

