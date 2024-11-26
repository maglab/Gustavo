### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Alp - Alpha (Transparency) (e.g. Alp)
# Ang - Angle (e.g., TckAng)
# Ann - Annotation (e.g., AnnTxtSiz)
# ARC - Ageing-related Community (e.g., ARC_Pleiotropy)
# Arr - Array (e.g., DtsArr)
# COX - Coexpresion (e.g., IndFrm_COX90)
# Dot - Dot (e.g., DotArr)
# Dts - Dataset (e.g., DtsArr)
# Ens - Ensembl (e.g., EnsHgnFrm)
# Frm - Frame (e.g., TauFrm)
# Gen - Gene (e.g., Gen)
# Hgn - HGNC (e.g., EnsHgnFrm)
# Jmp - Jump (e.g., YtstJmp)
# Ind - Indirect (e.g., IndFrm_PPI)
# Int - Interactor (e.g., TypInt)
# Len - Length (e.g., LenAnn)
# Lst - List (e.g., TauDirIndLst)
# Num - Number (e.g., NumThr)
# Pnt - Point (e.g., MinPnt)
# Qry - Query (e.g., DtsQry)
# Sep - Separison (e.g., TstSep)
# Siz = Size (e.g., TstTxtSiz)
# Tau - Tau (e.g., Tau)
# Tck - Tick (e.g., TckAng)
# Thr - Threshold (e.g., NumThr)
# Tst - Test (e.g., TstTxtSiz)
# Txt - Text (e.g.,TstTxtSiz )
# Typ - Type (e.g., TypInt)

### LIBRARIES ##################################################################

# Dir - Direct
# Ind - Indirect

library(dplyr)
library(cowplot)
library(ggplot2)
library(rstatix)
library(ggpubr)

### START ######################################################################

TauFrm = read.csv("maglab/Gustavo/Data/Retrieved/Specificity/Tau_gene_V8.csv") 
EnsHgnFrm = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/HgncEnsembl.rds")
row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
TauFrm$Gen = EnsHgnFrm[TauFrm$gene_id,"external_gene_name"]
TauFrm = TauFrm %>% filter(!is.na(tau)) %>% select(Gen, tau) %>% rename(Tau=tau)

DirIndIntFrm_PPI = read.csv('maglab/Gustavo/Data/Generated/PPI/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') 

# Pleiotropy
DirFrm = DirIndIntFrm_PPI %>% filter(ARC_Pleiotropy>0) %>% select(Gen,ARC_Pleiotropy) %>% rename(Num=ARC_Pleiotropy) %>%
  mutate(TypInt = "Dir") %>% 
  mutate(NumThr = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "DISEASE") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, NumThr, Dts)

# PPI iARC-Interactions
IndFrm_PPI = DirIndIntFrm_PPI %>% 
  filter(iARC_Interactor>0) %>% select(Gen,iARC_Interactor) %>% rename(Num=iARC_Interactor) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(NumThr = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "PPI") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, NumThr, Dts)

# COX90 iARC-Interactions
IndFrm_COX90 = read.csv('maglab/Gustavo/Data/Generated/COX90/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') %>%
  filter(iARC_Interactor>0) %>% select(Gen,iARC_Interactor) %>% rename(Num=iARC_Interactor) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(NumThr = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "COX90") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, NumThr, Dts)

# COX95 iARC-Interactions
IndFrm_COX95 = read.csv('maglab/Gustavo/Data/Generated/COX95/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv')  %>%
  filter(iARC_Interactor>0) %>% select(Gen,iARC_Interactor) %>% rename(Num=iARC_Interactor) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(NumThr = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "COX95") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, NumThr, Dts)

# KEGG iARC-Interactions
IndFrm_KEGG = read.csv('maglab/Gustavo/Data/Generated/KEGG/Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv') %>%
  filter(iARC_Interactor>0) %>% select(Gen,iARC_Interactor) %>% rename(Num=iARC_Interactor) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(NumThr = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "KEGG") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, NumThr, Dts)

# BINDING TAU FRAMES
TauDirIndFrm = rbind(DirFrm,IndFrm_PPI,IndFrm_COX90,IndFrm_COX95,IndFrm_KEGG)

# ARRAR AND LENGHT OF DATASETS
DtsArr = TauDirIndFrm$Dts %>% unique()
DtsLen = length(DtsArr)

# Plot parameters
Alp = 0.2
TckAng = 45
TstTxtSiz = 3
AnnTxtSiz = 3
YtstJmp = 0.1
TstSep = 0.07
Ylab = "Tau"
MinPnt = 1.05
DotArr = c(1,3.5,1,3,3.5,3.5)
LenAnn = 4


# Plot list and frame
TauDirIndLst = list()
MeanMedianFrm = data.frame()
for(z in 1:DtsLen){
  
  print(z)
  
  DtsQry = DtsArr[z] 
  
  TypIntQry = TauDirIndFrm %>%
    filter(Dts == DtsQry) %>%
    pull(TypInt) %>%
    unique()
  
  if(DtsQry == "DISEASE"){
    DtsTxt = "Pleiotropy"
  }
  if(DtsQry == "PPI"){
    DtsTxt = "PPI"
  }
  if(DtsQry == "KEGG"){
    DtsTxt = "KEGG"
  }
  if(DtsQry == "COX90"){
    DtsTxt = "COX.90"
  }
  if(DtsQry == "COX95"){
    DtsTxt = "COX.95"
  }
  
  TstFrm = TauDirIndFrm %>% 
    filter(Dts == DtsQry) %>%
    select(NumThr,Tau,Dts)
  
  pwc = TstFrm %>% pairwise_t_test( Tau ~ NumThr , p.adjust.method = "bonferroni")

  if(TypIntQry == "Dir"){
    TstFrm$TypInt = paste("Tissue specificity\nARC-Pleiotropy", sep="")
    Xlab = "Number of directly\nconnected ARCs"
  } else{
    TstFrm$TypInt = paste("Tissue specificity\niARC-Interactor (",DtsTxt,")", sep="")
    Xlab = "Number of indirectly\nconnected ARCs"
  }
  
  # Statistics
  grouped_means <-  TstFrm %>%
    group_by(NumThr) %>%
    summarise(MeanTau = mean(Tau, na.rm = TRUE))
  grouped_means
  
  
  grouped_medians <- TstFrm %>%
    group_by(NumThr) %>%
    summarise(MedianTau = median(Tau, na.rm = TRUE))
  
  # Mean Tau Frame
  sMeanMedianFrm = grouped_means %>% merge(grouped_medians, by='NumThr')
  sMeanMedianFrm$Dts = DtsQry
  MeanMedianFrm = rbind(MeanMedianFrm,sMeanMedianFrm)
  
  # Convert the `NumThr` column to a factor with the desired order
  TstFrm$NumThr <- factor(TstFrm$NumThr, levels = c("1", "2", "3", "4+"))
  
  # BOX PLOT 
    p = ggboxplot(TstFrm, x = "NumThr", y = "Tau",
                  fill="#F8776D",
                  size = 0.9,
                  facet.by = "TypInt",
                  combine = TRUE,
                  outlier.shape = NA,
    ) 
  
  
  p = p+geom_dotplot(binaxis='y', stackdir='center',
                     position=position_dodge(1), dotsize=0.004*(DotArr[z]), color="#202020")

  p = p + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) 
  
  p = ggpar(p, xlab =Xlab, ylab = Ylab)+
    font("xlab", size = 8, color = "black")+
    font("ylab", size = 8, color = "black")+
    font("xy.text", size = 8, color = "black") +
    theme(legend.position = "none")
  
  p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL) #x, y, xy
  
  pwc$y.position = c(MinPnt, MinPnt+(1*YtstJmp), MinPnt, MinPnt+(2*YtstJmp), MinPnt+(3*YtstJmp),MinPnt)
  pwc[[".y."]] = pwc$y.position
  
  # X Separison - Tests
  pwc = pwc %>% mutate(xmin = c(1,1,2,1,2,3)+TstSep )
  pwc = pwc %>% mutate(xmax = c(2,3,3,4,4,4)-TstSep )
  
  p = p + stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0, step.increase = 0.00, label.size = TstTxtSiz) 
  
  p = p +
    annotate("text",
             x = 1:LenAnn,
             y = -0.1,
             label = TstFrm$NumThr %>% table() %>% format_with_commas(),
             col = "black",
             size = AnnTxtSiz,
             vjust = - 1)
  
  p = p + coord_cartesian(ylim = c(-0,1.36))
  
  if(!(z %in% c(2,5))){
    p = p + ylab(NULL)
  }
  
  p = p + scale_x_discrete(labels = c("1" = "1", "2" = "2","3" = "3", "4+" = "4+") )
  
  TauDirIndLst[[DtsQry]] = p

}

# PLOTTING GRID
plot_grid(TauDirIndLst$PPI, TauDirIndLst$COX95, TauDirIndLst$COX90, 
          TauDirIndLst$KEGG, TauDirIndLst$DISEASE, NULL,
          labels = c("a","b","c","d","e","","","",""), label_size = 12, ncol=3) 


# ADAPT TAUS AND MEAN-MEDIAN FRAMES FOR SAVING
TauDirIndFrm = TauDirIndFrm %>% 
  rename(Gene=Gen,
         InteractionType=TypInt,
         ARC_Connected = Num,
         ARC_Connected_Threshold = NumThr,
         Dataset = Dts)

MeanMedianFrm = MeanMedianFrm %>%
  rename(ARC_Connected_Threshold=NumThr,
         Mean_Tau=MeanTau,
         Median_Tau=MedianTau,
         Dataset = Dts)


write.csv(TauDirIndFrm,"maglab/Gustavo/Data/Generated/Specificity/Pleiotropies_Interactors_Scores.csv",row.names=FALSE)
write.csv(MeanMedianFrm,"maglab/Gustavo/Data/Generated/Specificity/Pleiotropies_Interactors_MeanMedian.csv",row.names=FALSE)


################################################################################
# FUNCTIONS
################################################################################

format_with_commas <- function(x) {
  # Check if x is a numeric vector or a matrix
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector or matrix.")
  }
  
  # Format each element in x with commas
  formatted <- format(x, big.mark = ",", scientific = FALSE, justify = "none")
  
  # Remove leading spaces by trimming the formatted strings
  formatted <- trimws(formatted)
  
  # Return the formatted vector
  return(formatted)
}
